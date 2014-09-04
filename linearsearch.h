/*{{{
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

#include <vector>
#include "simdize.h"

// class LinearNeighborSearch {{{1
template <typename T> class LinearNeighborSearch
{
public:
    // insert {{{2
    LinearNeighborSearch() = default;
    LinearNeighborSearch(std::size_t reserve) { m_data.reserve(reserve); }

    // insert {{{2
    template <typename U> void insert(U &&x) { m_data.push_back(std::forward<U>(x)); }

    // findNearest {{{2
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

private:  // {{{2
    std::vector<T> m_data;
    // }}}2
};

// class LinearNeighborSearchV {{{1
template <typename T> class LinearNeighborSearchV
{
    using V = simdize<T>;

public:
    // insert {{{2
    LinearNeighborSearchV() = default;
    LinearNeighborSearchV(std::size_t reserve) { m_data.reserve(reserve); }

    // insert {{{2
    void insert(const T &x) {
        if (m_lastFill == 0) {
            m_data.push_back(V(x));
        } else {
            V &last = m_data.back();
            simdize_assign(last, m_lastFill, x);
        }
        ++m_lastFill;
        m_lastFill %= V::Size;
    }

    // findNearest {{{2
    T findNearest(const T &x_) const
    {
        const V x(x_);
        auto it = m_data.begin();
        const auto end = m_data.end();
        V p2 = *it;
        auto bestDistance = get_kdtree_distance(x, *it);
        for (++it; it != end; ++it) {
            const auto tmp = get_kdtree_distance(x, *it);
            where(tmp < bestDistance) | p2 = *it;
            bestDistance = min(tmp, bestDistance);
        }
        return simdize_get(p2, (bestDistance.min() == bestDistance).firstOne());
    }

private:  // {{{2
    std::vector<simdize<T>> m_data;
    int m_lastFill = 0;
    // }}}2
};

// vim: foldmethod=marker
