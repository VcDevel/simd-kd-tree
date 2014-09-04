/*  This file is part of the Vc library. {{{
Copyright © 2014 Matthias Kretz <kretz@kde.org>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the names of contributing organizations nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

}}}*/

#ifndef KDTREE_H_
#define KDTREE_H_

#include "simdize.h"

// make_unique {{{1
template <typename T, typename... Args> std::unique_ptr<T> make_unique(Args &&... args)
{
    return std::unique_ptr<T>{new T{std::forward<Args>(args)...}};
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
            // if we have a child on the "wrong" side it could still contain the
            // nearest neighbor, but only if the shortest distance of the search point
            // x to the splitting plane is less than the distance to the current
            // candidate. We only need to look at dx[0] and dx[V::Size - 1]
            const auto dx = get_kdtree_1dim_distance<SplittingPlane>(x, *this);

            const auto less = get_kdtree_value<SplittingPlane>(*this) < get_kdtree_value<SplittingPlane>(x);
            if (all_of(less)) {
                if (m_child[0]) {
                    auto candidate = m_child[0]->findNearest(x);
                    if (dx.min() < candidate.second) {
                        const auto distance = get_kdtree_distance(x, *this);
                        if (distance.min() < candidate.second) {
                            candidate.second = distance.min();
                            candidate.first = simdize_get(
                                *this, (distance.min() == distance).firstOne());
                        }
                        if (m_child[1] && dx.max() < candidate.second) {
                            const auto candidate2 = m_child[1]->findNearest(x);
                            if (candidate2.second < candidate.second) {
                                return candidate2;
                            }
                        }
                    }
                    return candidate;
                }
            } else if (none_of(less)) {
                if (m_child[1]) {
                    auto candidate = m_child[1]->findNearest(x);
                    if (dx.min() < candidate.second) {
                        const auto distance = get_kdtree_distance(x, *this);
                        if (distance.min() < candidate.second) {
                            candidate.second = distance.min();
                            candidate.first = simdize_get(
                                *this, (distance.min() == distance).firstOne());
                        }
                        if (m_child[0] && dx.max() < candidate.second) {
                            const auto candidate2 = m_child[0]->findNearest(x);
                            if (candidate2.second < candidate.second) {
                                return candidate2;
                            }
                        }
                    }
                    return candidate;
                }
            }
            const auto distance = get_kdtree_distance(x, *this);
            const int pos = (distance.min() == distance).firstOne();
            auto candidate = std::make_pair(simdize_get(*this, pos), distance.min());
            const auto bestDx = dx[pos];
            if (bestDx <= candidate.second && m_child[0]) {
                const auto candidate2 = m_child[0]->findNearest(x);
                if (candidate2.second < candidate.second) {
                    candidate = candidate2;
                }
            }
            if (bestDx <= candidate.second && m_child[1]) {
                const auto candidate2 = m_child[1]->findNearest(x);
                if (candidate2.second < candidate.second) {
                    return candidate2;
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
    using DistanceType = decltype(get_kdtree_distance(std::declval<T>(), std::declval<T>()));

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

        std::pair<T, DistanceType> findNearest(const T x) const
        {
            const auto dx = get_kdtree_1dim_distance<SplittingPlane>(x, m_data);

            const std::size_t index = (get_kdtree_value<SplittingPlane>(x) <
                                       get_kdtree_value<SplittingPlane>(m_data))
                                          ? 0
                                          : 1;
            if (m_child[index]) {
                // if we have a child node on the side where the search point would get
                // stored, that node is an obvious candidate
                auto candidate = m_child[index]->findNearest(x);
                if (dx < candidate.second) {
                    // if the resulting distance yields a sphere that intersects the current
                    // splitting plane, then we have to look further
                    const auto distance = get_kdtree_distance(x, m_data);
                    if (distance < candidate.second) {
                        candidate = std::make_pair(m_data, distance);
                    }
                    if (dx < candidate.second && m_child[index ^ 1]) {
                        const auto candidate2 = m_child[index ^ 1]->findNearest(x);
                        if (candidate2.second < candidate.second) {
                            return candidate2;
                        }
                    }
                }
                return candidate;
            }
            auto candidate = std::make_pair(m_data, get_kdtree_distance(x, m_data));
            if (m_child[index ^ 1] && dx < candidate.second) {
                // if we have a child on the "wrong" side it could still be/find the
                // nearest neighbor, but only if the shortest distance of the search point
                // x to the splitting plane is less than the distance to the current
                // candidate.
                const auto candidate2 = m_child[index ^ 1]->findNearest(x);
                if (candidate2.second < candidate.second) {
                    return candidate2;
                }
            }
            return candidate;
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

    T findNearest(const T x) const
    {
        if (!m_root) {
            throw std::runtime_error(
                "No values in the KdTree, which is required for findNearest.");
        }
        return m_root->findNearest(x).first;
    }

    friend std::ostream &operator<<(std::ostream &out, const KdTree &tree)
    {
        return out << *tree.m_root;
    }
};  //}}}1

#endif  // KDTREE_H_

// vim: foldmethod=marker
