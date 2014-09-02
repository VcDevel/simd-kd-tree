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

#include <array>
#include <type_traits>
#include <iostream>
#include <Vc/Vc>

// struct Point<T, N> {{{1
template <typename T, std::size_t N> class Point
{
    std::array<T, N> coordinate;

public:
    template <typename U0, typename... Us,
              typename =
                  typename std::enable_if<(!std::is_convertible<U0, Point>::value)>::type>
    Point(U0 &&init0, Us &&... init)
        : coordinate{{std::forward<U0>(init0), std::forward<Us>(init)...}}
    {
    }

    Point() = default;
    Point(const Point &) = default;
    Point(Point &&) = default;
    Point &operator=(const Point &) = default;
    Point &operator=(Point &&) = default;

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

