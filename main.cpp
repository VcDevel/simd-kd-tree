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

//#include <Vc/Vc>
#include <array>
#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <stdexcept>

#include "../tsc.h"

template <typename T, typename... Args> std::unique_ptr<T> make_unique(Args &&... args)
{
    return std::unique_ptr<T>{new T{std::forward<Args>(args)...}};
}

struct Point {
    //Point() = default;
    //Point(const Point &) = default;
    //Point(Point &&) = default;
    std::array<int, 3> coordinate;

    friend std::ostream &operator<<(std::ostream &out, const Point &p)
    {
        return out << '[' << p.coordinate[0] << ' ' << p.coordinate[1] << ' '
                   << p.coordinate[2] << ']';
    }

    friend int get_kdtree_distance(const Point &p0, const Point &p1)
    {
        const auto dx = p0.coordinate[0] - p1.coordinate[0];
        const auto dy = p0.coordinate[1] - p1.coordinate[1];
        const auto dz = p0.coordinate[2] - p1.coordinate[2];
        return dx * dx + dy * dy + dz * dz;
    }
};

template <std::size_t Plane> int get_kdtree_value(const Point &p)
{
    return p.coordinate[Plane];
}

template <std::size_t Dimensions, typename T> class KdTree
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
                const auto dx = get_kdtree_value<SplittingPlane>(x) -
                                get_kdtree_value<SplittingPlane>(m_data);
                if (dx * dx < get_kdtree_distance(x, *candidate)) {
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

int main()
{
    constexpr int SetSize = 20000;
    constexpr int NumberOfSearches = 50000;

    std::default_random_engine randomEngine(1);
    std::uniform_int_distribution<int> uniform(0, 99);

    KdTree<3, Point> pointsTree;
    LinearNeighborSearch<Point> pointsVector(SetSize);
    for (int i = 0; i < SetSize; ++i) {
        const Point p{{{uniform(randomEngine), uniform(randomEngine), uniform(randomEngine)}}};
        pointsTree.insert(p);
        pointsVector.insert(p);
    }
    //std::cout << pointsTree << '\n';

    std::vector<Point<float>> searchPoints;
    searchPoints.reserve(NumberOfSearches);
    for (int i = 0; i < NumberOfSearches; ++i) {
        searchPoints.push_back({{{uniform(randomEngine), uniform(randomEngine), uniform(randomEngine)}}});
    }

    TimeStampCounter tsc;
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
        const auto &p2 = pointsVector.findNearest(p);
        asm(""::"m"(p2));
    }
    tsc.stop();
    const auto time_linear = tsc.cycles();

    std::cout << "kd-tree: " << std::setw(11) << time_kdtree
              << "\n linear: " << std::setw(11) << time_linear
              << "\n  ratio: " << std::setw(11)
              << double(time_linear) / double(time_kdtree) << '\n';

    return 0;
}
