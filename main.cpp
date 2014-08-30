/*  This file is part of the Vc library.

    Copyright (C) 2014 Matthias Kretz <kretz@kde.org>

    Permission to use, copy, modify, and distribute this software
    and its documentation for any purpose and without fee is hereby
    granted, provided that the above copyright notice appear in all
    copies and that both that the copyright notice and this
    permission notice and warranty disclaimer appear in supporting
    documentation, and that the name of the author not be used in
    advertising or publicity pertaining to distribution of the
    software without specific, written prior permission.

    The author disclaim all warranties with regard to this
    software, including all implied warranties of merchantability
    and fitness.  In no event shall the author be liable for any
    special, indirect or consequential damages or any damages
    whatsoever resulting from loss of use, data or profits, whether
    in an action of contract, negligence or other tortious action,
    arising out of or in connection with the use or performance of
    this software.

*/

#include <array>
#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <stdexcept>

#include "../tsc.h"
#include "simdize.h"

// make_unique {{{1
template <typename T, typename... Args> std::unique_ptr<T> make_unique(Args &&... args)
{
    return std::unique_ptr<T>{new T{std::forward<Args>(args)...}};
}

// struct Point<T, N> {{{1
template <typename T, std::size_t N> class Point
{
    std::array<T, N> coordinate;

public:
    template <typename... Us> Point(Us &&... init) : coordinate{{std::forward<Us>(init)...}}
    {
    }

    Point() = default;

    T &operator[](std::size_t i) noexcept { return coordinate[i]; }
    const T &operator[](std::size_t i) const noexcept { return coordinate[i]; }
};

template <typename T, std::size_t N> std::ostream &operator<<(std::ostream &out, const Point<T, N> &p)
{
    out << '[' << p[0];
    Vc::Common::unrolled_loop<std::size_t, 1, N>([&](std::size_t i) {
        out << ' ' << p[i];
    });
    return out << ']';
}

// tuple interface to Point<T> {{{1
namespace std
{
template <typename T, std::size_t N>
struct tuple_size<Point<T, N>> : public std::integral_constant<std::size_t, N>
{
};
template <std::size_t I, typename T, std::size_t N> struct tuple_element<I, Point<T, N>>
{
    typedef T type;
};
}  // namespace std
template <std::size_t I, typename T, std::size_t N> T &get(Point<T, N> &x) noexcept
{
    return x[I];
}
template <std::size_t I, typename T, std::size_t N>
const T &get(const Point<T, N> &x) noexcept
{
    return x[I];
}

// get_kdtree_distance {{{1
template <typename T, std::size_t N>
T get_kdtree_distance(const Point<T, N> &p0, const Point<T, N> &p1)
{
    const auto dx = p0[0] - p1[0];
    T r = dx * dx;
    Vc::Common::unrolled_loop<std::size_t, 1, N>([&](std::size_t i) {
        const auto d_ = p0[i] - p1[i];
        r += d_ * d_;
    });
    return r;
}

// get_kdtree_value {{{1
template <std::size_t Plane, typename T, std::size_t N>
T get_kdtree_value(const Point<T, N> &p)
{
    return p[Plane];
}

// get_kdtree_1dim_distance {{{1
template <std::size_t Plane, typename T, std::size_t N>
T get_kdtree_1dim_distance(const Point<T, N> &p0, const Point<T, N> &p1)
{
    const auto dx = get_kdtree_value<Plane>(p0) - get_kdtree_value<Plane>(p1);
    return dx * dx;
}

// class KdTreeV {{{1
template <typename T, std::size_t Dimensions = std::tuple_size<T>::value> class KdTreeV
{
    using V = simdize<T>;
    using OneDimV = typename V::FirstVectorType;

    template <std::size_t SplittingPlane> struct Node;
    template <std::size_t SplittingPlane>
    using NodePtr = std::unique_ptr<Node<SplittingPlane>>;

    template <std::size_t SplittingPlane> struct Node : public V {
        // Node ChildSplittingPlane {{{2
        /// the template parameter to child Node objects
        static constexpr std::size_t ChildSplittingPlane =
            (SplittingPlane + 1) % Dimensions;
        // Node ChildPtr {{{2
        /// unique_ptr type for the children
        using ChildPtr = NodePtr<ChildSplittingPlane>;

        // Node data members {{{2
        std::array<ChildPtr, 2> m_child;
        unsigned int m_entries;

        // Node(T) {{{2
        /// create a new leaf node with payload \p x
        Node(const T &x)
            : V(x)     // broadcast to all entries
            , m_entries(1)  // but mark that only the first entry is valid
        {}

        // Node::insert {{{2
        /// check whether \p x needs to go left or right and hand it off or create a new
        /// leaf node
        void insert(const T &x)
        {
            using namespace std;
            const auto less = get_kdtree_value<SplittingPlane>(*this) <
                              get_kdtree_value<SplittingPlane>(V(x));
            //cerr << SplittingPlane << '|' << setw(10) << x[SplittingPlane] << ' ' << less << ' ' << *this[SplittingPlane];
            if (m_entries < V::Size) {
                if (any_of(less)) {
                    if (all_of(less)) {
                        //cerr << " (1)";
                        for (std::size_t i = 0; i < Dimensions; ++i) {
                            (*this)[i][m_entries] = x[i];
                        }
                    } else {
                        int pos = (!less).firstOne();
                        //cerr << " (2) " << pos;
                        for (std::size_t i = 0; i < Dimensions; ++i) {
                            where(OneDimV::IndexesFromZero() > pos) | (*this)[i] =
                                (*this)[i].shifted(-1);
                            (*this)[i][pos] = x[i];
                        }
                    }
                } else {
                    //cerr << " (3)";
                    for (std::size_t i = 0; i < Dimensions; ++i) {
                        (*this)[i] = (*this)[i].shifted(-1);
                        (*this)[i][0] = x[i];
                    }
                }
                ++m_entries;
                //cerr << ' ' << (*this)[SplittingPlane] << ' ' << m_entries << '\n';
            } else {
                if (all_of(less)) {
                    //cerr << " (4)\n";
                    // needs to go to the left child node
                    if (m_child[0]) {
                        m_child[0]->insert(x);
                    } else {
                        m_child[0] = make_unique<Node<ChildSplittingPlane>>(x);
                    }
                } else if (none_of(less)) {
                    //cerr << " (5)\n";
                    // needs to go to the right child node
                    if (m_child[1]) {
                        m_child[1]->insert(x);
                    } else {
                        m_child[1] = make_unique<Node<ChildSplittingPlane>>(x);
                    }
                } else {
                    // The new value must go into this node. The value it pushes out of (*this) could
                    // go to either the left or the right child - depending on what side we want to
                    // push out. Determine the bias from the position.
                    std::size_t pos = (!less).firstOne();
                    T new_x;
                    if (pos <= V::Size / 2) {
                        pos -= 1;
                        //cerr << " (6) " << pos;
                        for (std::size_t i = 0; i < Dimensions; ++i) {
                            // store largest value (for this SplittingPlane)
                            new_x[i] = (*this)[i][0];
                            where(less) | (*this)[i] = (*this)[i].shifted(1);
                            (*this)[i][pos] = x[i];
                        }
                        //cerr << ' ' << (*this)[SplittingPlane] << '\n';
                        if (m_child[1]) {
                            m_child[1]->insert(new_x);
                        } else {
                            m_child[1] = make_unique<Node<ChildSplittingPlane>>(new_x);
                        }
                    } else {
                        //cerr << " (7) " << pos;
                        for (std::size_t i = 0; i < Dimensions; ++i) {
                            // store smallest value (for this SplittingPlane)
                            new_x[i] = (*this)[i][V::Size - 1];
                            where(!less) | (*this)[i] = (*this)[i].shifted(-1);
                            (*this)[i][pos] = x[i];
                        }
                        //cerr << ' ' << (*this)[SplittingPlane] << '\n';
                        if (m_child[0]) {
                            m_child[0]->insert(new_x);
                        } else {
                            m_child[0] = make_unique<Node<ChildSplittingPlane>>(new_x);
                        }
                    }
                }
            }
        }

        // Node::operator<< {{{2
        friend std::ostream &operator<<(std::ostream &out, const Node &node)
        {
            out << SplittingPlane << ' ' << static_cast<const V &>(node);
            if (node.m_child[0]) {
                out << "\nl" << *node.m_child[0];
            }
            if (node.m_child[1]) {
                out << "\nr" << *node.m_child[1];
            }
            return out << " up";
        }

        // Node::findNearest {{{2
        T findNearest(const T &x) const
        {
            using namespace std;
            const V xv(x);

            // if we have a child on the "wrong" side it could still be/find the
            // nearest neighbor, but only if the shortest distance of the search point
            // x to the splitting plane is less than the distance to the current
            // candidate. We only need to look at dx[0] and dx[V::Size - 1]
            const auto dx = get_kdtree_1dim_distance<SplittingPlane>(xv, *this);

            const auto less = get_kdtree_value<SplittingPlane>(*this) < get_kdtree_value<SplittingPlane>(xv);
            const auto d = get_kdtree_distance(xv, *this);
            const auto pos = (d.min() == d).firstOne();
            T candidate = simdize_get(*this, pos);
            auto bestDistance = d[pos];
            if (all_of(less) && m_child[0]) {
                const auto tmp = m_child[0]->findNearest(x);
                const auto tmpDistance = get_kdtree_distance(x, tmp);
                if (tmpDistance < bestDistance) {
                    candidate = tmp;
                    bestDistance = tmpDistance;
                }
            } else if (none_of(less) && m_child[1]) {
                const auto tmp = m_child[1]->findNearest(x);
                const auto tmpDistance = get_kdtree_distance(x, tmp);
                if (tmpDistance < bestDistance) {
                    candidate = tmp;
                    bestDistance = tmpDistance;
                }
            }
            if (!all_of(less) &&  // no need to go down the same path again
                dx[V::Size - 1] < bestDistance && m_child[0]) {
                const auto tmp = m_child[0]->findNearest(x);
                const auto tmpDistance = get_kdtree_distance(x, tmp);
                if (tmpDistance < bestDistance) {
                    candidate = tmp;
                    bestDistance = tmpDistance;
                }
            }
            if (!none_of(less) &&  // no need to go down the same path again
                dx[0] < bestDistance && m_child[1]) {
                const auto tmp = m_child[1]->findNearest(x);
                const auto tmpDistance = get_kdtree_distance(x, tmp);
                if (tmpDistance < bestDistance) {
                    candidate = tmp;
                    // bestDistance = tmpDistance;
                }
            }

            return candidate;
        }
    };

    NodePtr<0> m_root; //{{{2

public: //{{{2
    KdTreeV() = default; //{{{2

    // insert {{{2
    template <typename U> void insert(U &&x)
    {
        if (m_root) {
            m_root->insert(std::forward<U>(x));
        } else {
            m_root = make_unique<Node<0>>(std::forward<U>(x));
        }
    }

    // findNearest {{{2
    T findNearest(const T &x) const
    {
        if (!m_root) {
            throw std::runtime_error(
                "No values in the KdTree, which is required for findNearest.");
        }
        return m_root->findNearest(x);
    }

    // operator<< {{{2
    friend std::ostream &operator<<(std::ostream &out, const KdTreeV &tree)
    {
        return out << *tree.m_root;
    }
    // }}}2
};

// class KdTree {{{1
template <typename T, std::size_t Dimensions = std::tuple_size<T>::value> class KdTree
{
    template <std::size_t SplittingPlane> struct Node;
    template <std::size_t SplittingPlane>
    using NodePtr = std::unique_ptr<Node<SplittingPlane>>;

    template <std::size_t SplittingPlane> struct Node {
        /// the template parameter to child Node objects
        static constexpr std::size_t ChildSplittingPlane =
            (SplittingPlane + 1) % Dimensions;
        /// unique_ptr type for the children
        using ChildPtr = NodePtr<ChildSplittingPlane>;

        T m_data;
        std::array<ChildPtr, 2> m_child;

        /// create a new leaf node with payload \p x
        template <typename U> Node(U &&x) : m_data(std::forward<U>(x)) {}

        /// check whether \p x needs to go left or right and hand it off or create a new
        /// leaf node
        template <typename U> void insert(U &&x)
        {
            ChildPtr &child = (get_kdtree_value<SplittingPlane>(x) <
                               get_kdtree_value<SplittingPlane>(m_data))
                                  ? m_child[0]
                                  : m_child[1];
            if (child) {
                child->insert(std::forward<U>(x));
            } else {
                child = make_unique<Node<ChildSplittingPlane>>(std::forward<U>(x));
            }
        }

        friend std::ostream &operator<<(std::ostream &out, const Node &node)
        {
            out << SplittingPlane << ' ' << node.m_data;
            if (node.m_child[0]) {
                out << "\nl" << *node.m_child[0];
            }
            if (node.m_child[1]) {
                out << "\nr" << *node.m_child[1];
            }
            return out << " up";
        }

        static const T &selectCloser(const T &search, const T &candidate0,
                                     const T &candidate1)
        {
            return get_kdtree_distance(search, candidate0) <
                           get_kdtree_distance(search, candidate1)
                       ? candidate0
                       : candidate1;
        }

        const T &findNearest(const T &x) const
        {
            const std::size_t index = (get_kdtree_value<SplittingPlane>(x) <
                                       get_kdtree_value<SplittingPlane>(m_data))
                                          ? 0
                                          : 1;
            const T *candidate = &m_data;  // first candidate to return is our point
            if (m_child[index]) {
                // if we have a child node on the side where the search point would get
                // stored, that node is an obvious candidate. Compare whether the existing
                // candidate is still better, though.
                candidate = &selectCloser(x, m_child[index]->findNearest(x), *candidate);
            }
            if (m_child[index ^ 1]) {
                // if we have a child on the "wrong" side it could still be/find the
                // nearest neighbor, but only if the shortest distance of the search point
                // x to the splitting plane is less than the distance to the current
                // candidate.
                const auto dx = get_kdtree_1dim_distance<SplittingPlane>(x, m_data);
                if (dx < get_kdtree_distance(x, *candidate)) {
                    candidate =
                        &selectCloser(x, m_child[index ^ 1]->findNearest(x), *candidate);
                }
            }
            return *candidate;
        }
    };

    NodePtr<0> m_root;

public:
    KdTree() = default;

    template <typename U> void insert(U &&x)
    {
        if (m_root) {
            m_root->insert(std::forward<U>(x));
        } else {
            m_root = make_unique<Node<0>>(std::forward<U>(x));
        }
    }

    const T &findNearest(const T &x) const
    {
        if (!m_root) {
            throw std::runtime_error(
                "No values in the KdTree, which is required for findNearest.");
        }
        return m_root->findNearest(x);
    }

    friend std::ostream &operator<<(std::ostream &out, const KdTree &tree)
    {
        return out << *tree.m_root;
    }
};

// class LinearNeighborSearch {{{1
template <typename T> class LinearNeighborSearch
{
public:
    LinearNeighborSearch() = default;
    LinearNeighborSearch(std::size_t reserve) { m_data.reserve(reserve); }

    template <typename U> void insert(U &&x) { m_data.push_back(std::forward<U>(x)); }
    T findNearest(const T &x) const
    {
        /*
        auto bestDistance = get_kdtree_distance(m_data[0], x);
        std::size_t best = 0;
        for (std::size_t i = 1; i < m_data.size(); ++i) {
            const auto d = get_kdtree_distance(m_data[i], x);
            if (d < bestDistance) {
                bestDistance = d;
                best = i;
            }
        }
        return m_data[best];
        */
        using TT = decltype(get_kdtree_distance(x, x));
        TT bestDistance = std::numeric_limits<TT>::max();
        T p2;
        for (const auto &p : m_data) {
            const auto d = get_kdtree_distance(p, x);
            if (d < bestDistance) {
                bestDistance = d;
                p2 = p;
            }
        }
        return p2;
    }

private:
    std::vector<T> m_data;
};

int main() //{{{1
{
    constexpr int SetSize = 20000;
    constexpr int NumberOfSearches = 50000;

    using T = float;
    using Point = ::Point<T, 3>;

    std::default_random_engine randomEngine(1);
    std::uniform_real_distribution<T> uniform(-99, 99);

    std::vector<Point> randomPoints;
    for (int i = 0; i < SetSize; ++i) {
        randomPoints.emplace_back(uniform(randomEngine), uniform(randomEngine), uniform(randomEngine));
    }

    TimeStampCounter tsc;

    tsc.start();
    KdTree<Point> pointsTree;
    for (const Point &p : randomPoints) {
        pointsTree.insert(p);
    }
    tsc.stop();
    const auto kdtree_inserts = tsc.cycles();

    tsc.start();
    KdTreeV<Point> pointsTreeV;
    for (const Point &p : randomPoints) {
        pointsTreeV.insert(p);
    }
    tsc.stop();
    const auto kdtreev_inserts = tsc.cycles();

    tsc.start();
    LinearNeighborSearch<Point> pointsVector(SetSize);
    for (const Point &p : randomPoints) {
        pointsVector.insert(p);
    }
    tsc.stop();
    const auto linear_inserts = tsc.cycles();
    std::cout << "-- inserts --"
              << "\n   kd-tree: " << std::setw(11) << kdtree_inserts
              << "\nVc kd-tree: " << std::setw(11) << kdtreev_inserts
              << "\n    linear: " << std::setw(11) << linear_inserts
              << "\n  Vc ratio: " << std::setw(11)
              << double(kdtree_inserts) / double(kdtreev_inserts) << '\n';

    //std::cout << pointsTreeV << '\n';

    std::vector<Point> searchPoints;
    searchPoints.reserve(NumberOfSearches);
    for (int i = 0; i < NumberOfSearches; ++i) {
        searchPoints.emplace_back(uniform(randomEngine), uniform(randomEngine), uniform(randomEngine));
    }

    tsc.start();
    for (int i = 0; i < NumberOfSearches; ++i) {
        const auto &p = searchPoints[i];
        const auto &p2 = pointsTree.findNearest(p);
        asm(""::"m"(p2));
        //std::cout << "looking near " << p << ", found " << p2 << ", distance " << get_kdtree_distance(p, p2) << '\n';
    }
    tsc.stop();
    const auto time_kdtree = tsc.cycles();

    tsc.start();
    for (int i = 0; i < NumberOfSearches; ++i) {
        const auto &p = searchPoints[i];
        const auto &p2 = pointsTreeV.findNearest(p);
        asm(""::"m"(p2));
        //std::cout << "looking near " << p << ", found " << p2 << ", distance " << get_kdtree_distance(p, p2) << '\n';
    }
    tsc.stop();
    const auto time_kdtreev = tsc.cycles();

    tsc.start();
    for (int i = 0; i < NumberOfSearches; ++i) {
        const auto &p = searchPoints[i];
        const auto &p2 = pointsVector.findNearest(p);
        asm(""::"m"(p2));
//#define COMPARE_KDTREE_LINEAR
#ifdef COMPARE_KDTREE_LINEAR
        const auto &p3 = pointsTree.findNearest(p);
        const auto &p4 = pointsTreeV.findNearest(p);
        if (get_kdtree_distance(p, p2) != get_kdtree_distance(p, p3)) {
            std::cerr << p << " failed " << p2 << " (" << get_kdtree_distance(p, p2)
                      << ") vs. KdTree " << p3 << " (" << get_kdtree_distance(p, p3) << ")\n";
        }
        if (get_kdtree_distance(p, p2) != get_kdtree_distance(p, p4)) {
            std::cerr << p << " failed " << p2 << " (" << get_kdtree_distance(p, p2)
                      << ") vs. KdTreeV " << p4 << " (" << get_kdtree_distance(p, p4) << ")\n";
        }
#endif
    }
    tsc.stop();
    const auto time_linear = tsc.cycles();

    std::cout << "-- searches --"
              << "\n   kd-tree: " << std::setw(11) << time_kdtree
              << "\nVc kd-tree: " << std::setw(11) << time_kdtreev
              << "\n    linear: " << std::setw(11) << time_linear
              << "\n     ratio: " << std::setw(11)
              << double(time_linear) / double(time_kdtree)
              << "\n  Vc ratio: " << std::setw(11)
              << double(time_kdtree) / double(time_kdtreev) << '\n';

    return 0;
} //}}}1

// vim: foldmethod=marker
