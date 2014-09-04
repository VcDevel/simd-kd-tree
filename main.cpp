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
#include "point.h"
#include "linearsearch.h"

// make_unique {{{1
template <typename T, typename... Args> std::unique_ptr<T> make_unique(Args &&... args)
{
    return std::unique_ptr<T>{new T{std::forward<Args>(args)...}};
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
const T &get_kdtree_value(const Point<T, N> &p) noexcept
{
    return p[Plane];
}
template <std::size_t Plane, typename T, std::size_t N>
T &get_kdtree_value(Point<T, N> &p) noexcept
{
    return p[Plane];
}

template <typename T, std::size_t N>
const T &get_kdtree_value(const Point<T, N> &p, std::size_t Plane) noexcept
{
    return p[Plane];
}
template <typename T, std::size_t N>
T &get_kdtree_value(Point<T, N> &p, std::size_t Plane) noexcept
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
    using DistanceType = decltype(get_kdtree_distance(std::declval<T>(), std::declval<T>()));

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
        Node(T x)
            : V(x)     // broadcast to all entries
            , m_entries(1)  // but mark that only the first entry is valid
        {}

        // Node::childInsert {{{2
        void childInsert(ChildPtr &child, const T new_x)
        {
            if (child) {
                child->insert(new_x);
            } else {
                // the Node constructor could be optimized to not have to redo the
                // V(x) broadcast
                child = make_unique<Node<ChildSplittingPlane>>(new_x);
            }
        };

        // Node::insert {{{2
        /// check whether \p x needs to go left or right and hand it off or create a new
        /// leaf node
        void insert(const T x)
        {
            using namespace std;
            //cerr << SplittingPlane << '|' << setw(10) << get_kdtree_value<SplittingPlane>(x)
            if (m_entries < V::Size) {
                //cerr << " (1)";
                simdize_assign(*this, m_entries, x);
                ++m_entries;
                //cerr << ' ' << get_kdtree_value<SplittingPlane>(*this) << ' ' << m_entries << '\n';
            } else {
                const auto less = get_kdtree_value<SplittingPlane>(*this) <
                                  get_kdtree_value<SplittingPlane>(V(x));
                //cerr << ' ' << less << ' ' << get_kdtree_value<SplittingPlane>(*this);
                if (all_of(less)) {  // needs to go to the left child node
                    //cerr << " (2)\n";
                    childInsert(m_child[0], x);
                } else if (none_of(less)) {  // needs to go to the right child node
                    //cerr << " (3)\n";
                    childInsert(m_child[1], x);
                } else {
                    // The new value must go into this node. The value it pushes out of (*this) could
                    // go to either the left or the right child - depending on which value we want to
                    // push out. Determine the bias from the number of values that are "less".
                    // If most values are "less" then the value to insert is larger than the median
                    // and thus we push out the max value. Otherwise push out the min value.
                    const bool majorityIsLess = less.count() >= int(V::Size / 2);
                    const auto &data = get_kdtree_value<SplittingPlane>(*this);
                    const auto replaceValue = majorityIsLess ? data.max() : data.min();
                    // replaceValue is constructed such that:
                    // all_of(get_kdtree_value<SplittingPlane>(*this) < get_kdtree_value<SplittingPlane>(V(new_x)))
                    // if the majority of values in this node are less than x (for the current
                    // SplittingPlane)
                    const std::size_t pos = (replaceValue == data).firstOne();
                    //cerr << " (4) " << pos << ' ' << replaceValue;
                    // swap largest/smallest value (of this SplittingPlane):
                    const T new_x = simdize_get(*this, pos);
                    simdize_assign(*this, pos, x);
                    //cerr << ' ' << get_kdtree_value<SplittingPlane>(*this) << '\n';
                    childInsert(m_child[majorityIsLess ? 0 : 1], new_x);
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
        std::pair<T, DistanceType> findNearest(const V &x) const
        {
            using namespace std;

            // if we have a child on the "wrong" side it could still contain the
            // nearest neighbor, but only if the shortest distance of the search point
            // x to the splitting plane is less than the distance to the current
            // candidate. We only need to look at dx[0] and dx[V::Size - 1]
            const auto dx = get_kdtree_1dim_distance<SplittingPlane>(x, *this);

            const auto less = get_kdtree_value<SplittingPlane>(*this) < get_kdtree_value<SplittingPlane>(x);
            const auto distance = get_kdtree_distance(x, *this);
            const auto pos = (distance.min() == distance).firstOne();
            T candidate = simdize_get(*this, pos);
            auto bestDistance = distance[pos];
            auto recurseChild = [&](const ChildPtr &child) {
                if (child) {
                    const auto candidateFromChild = child->findNearest(x);
                    if (candidateFromChild.second < bestDistance) {
                        candidate = candidateFromChild.first;
                        bestDistance = candidateFromChild.second;
                    }
                }
            };
            if (all_of(less)) {
                recurseChild(m_child[0]);
                if (dx.max() <= bestDistance) {
                    recurseChild(m_child[1]);
                }
            } else if (none_of(less)) {
                recurseChild(m_child[1]);
                if (dx.min() <= bestDistance) {
                    recurseChild(m_child[0]);
                }
            } else {
                const auto bestDx = dx[pos];
                if (bestDx <= bestDistance) {
                    recurseChild(m_child[0]);
                }
                if (bestDx <= bestDistance) {
                    recurseChild(m_child[1]);
                }
            }

            return std::make_pair(candidate, bestDistance);
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
    T findNearest(T x) const
    {
        if (!m_root) {
            throw std::runtime_error(
                "No values in the KdTree, which is required for findNearest.");
        }
        return m_root->findNearest(V(x)).first;
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

        const T &findNearest(const T x) const
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

    const T &findNearest(const T x) const
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

int main() //{{{1
{
    constexpr int SetSize = 20000;  // required memory ~ SetSize * sizeof(Node)
                                     // ~ SetSize * (sizeof(Point<T>) + 16)
                                     // = SetSize * 32
    constexpr int NumberOfSearches = 50000;

    using T = float;
    using Point = ::Point<T, 3>;

    std::default_random_engine randomEngine(1);
    typename std::conditional<std::is_floating_point<T>::value,
                              std::uniform_real_distribution<T>,
                              std::uniform_int_distribution<T>>::type uniform(-99, 99);

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
    //std::cout << pointsTreeV << '\n';

    tsc.start();
    LinearNeighborSearch<Point> pointsVector(SetSize);
    for (const Point &p : randomPoints) {
        pointsVector.insert(p);
    }
    tsc.stop();
    const auto linear_inserts = tsc.cycles();

    tsc.start();
    LinearNeighborSearchV<Point> linearSearchV(SetSize);
    for (const Point &p : randomPoints) {
        linearSearchV.insert(p);
    }
    tsc.stop();
    const auto linearv_inserts = tsc.cycles();

    std::cout << "             KdTree    KdTreeV  LinearNeighborSearch  LinearNeighborSearchV  KdTree/KdTreeV  "
                 "Linear/KdTree  Linear/LinearV\n";
    std::cout << "inserts "
              << std::setw(11) << kdtree_inserts
              << std::setw(11) << kdtreev_inserts
              << std::setw(22) << linear_inserts
              << std::setw(23) << linearv_inserts
              << std::setw(16) << double(kdtree_inserts) / double(kdtreev_inserts) << std::endl;

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

    tsc.start();
    for (int i = 0; i < NumberOfSearches; ++i) {
        const auto &p = searchPoints[i];
        const auto &p2 = linearSearchV.findNearest(p);
        asm(""::"m"(p2));
//#define COMPARE_KDTREE_LINEAR
#ifdef COMPARE_KDTREE_LINEAR
        const auto &p4 = pointsTreeV.findNearest(p);
        if (get_kdtree_distance(p, p2) != get_kdtree_distance(p, p4)) {
            std::cerr << p << " failed " << p2 << " (" << get_kdtree_distance(p, p2)
                      << ") vs. KdTreeV " << p4 << " (" << get_kdtree_distance(p, p4) << ")\n";
        }
#endif
    }
    tsc.stop();
    const auto time_linearv = tsc.cycles();

    std::cout << "searches"
              << std::setw(11) << time_kdtree
              << std::setw(11) << time_kdtreev
              << std::setw(22) << time_linear
              << std::setw(23) << time_linearv
              << std::setw(16) << double(time_kdtree) / double(time_kdtreev)
              << std::setw(15) << double(time_linear) / double(time_kdtree)
              << std::setw(16) << double(time_linear) / double(time_linearv) << '\n';

    return 0;
} //}}}1

// vim: foldmethod=marker
