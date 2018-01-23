/*{{{
Copyright © 2014-2018 Matthias Kretz <kretz@kde.org>

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
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

}}}*/

#include <array>
#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <stdexcept>

#include "../tsc.h"
#include "point.h"
#include "linearsearch.h"
#include "kdtree.h"

int main()  // {{{1
{
    // settings {{{2
    constexpr int SetSize = 20000;  // required memory ~ SetSize * sizeof(Node)
                                     // ~ SetSize * (sizeof(Point<T>) + 16)
                                     // = SetSize * 32
    constexpr int NumberOfSearches = 50000;

    using T = float;
    using Point = ::Point<T, 3>;

    // random points {{{2
    std::default_random_engine randomEngine(1);
    typename std::conditional<std::is_floating_point<T>::value,
                              std::uniform_real_distribution<T>,
                              std::uniform_int_distribution<T>>::type uniform(-99, 99);

    std::vector<Point> randomPoints;
    for (int i = 0; i < SetSize; ++i) {
        randomPoints.emplace_back(uniform(randomEngine), uniform(randomEngine), uniform(randomEngine));
    }

    TimeStampCounter tsc;

    tsc.start();  // create KdTree {{{2
    KdTree<Point> pointsTree;
    for (const Point &p : randomPoints) {
        pointsTree.insert(p);
    }
    tsc.stop();
    const auto kdtree_inserts = tsc.cycles();

    tsc.start();  // create KdTreeV {{{2
    KdTreeV<Point> pointsTreeV;
    for (const Point &p : randomPoints) {
        pointsTreeV.insert(p);
    }
    tsc.stop();
    const auto kdtreev_inserts = tsc.cycles();
    //std::cout << pointsTreeV << '\n';

    tsc.start();  // create LinearNeighborSearch {{{2
    LinearNeighborSearch<Point> pointsVector(SetSize);
    for (const Point &p : randomPoints) {
        pointsVector.insert(p);
    }
    tsc.stop();
    const auto linear_inserts = tsc.cycles();

    tsc.start();  // create LinearNeighborSearchV {{{2
    LinearNeighborSearchV<Point> linearSearchV(SetSize);
    for (const Point &p : randomPoints) {
        linearSearchV.insert(p);
    }
    tsc.stop();
    const auto linearv_inserts = tsc.cycles();

    // print insert timings {{{2
    std::cout << "             KdTree    KdTreeV  LinearNeighborSearch  LinearNeighborSearchV  KdTree/KdTreeV  "
                 "Linear/KdTree  Linear/LinearV\n";
    std::cout << "inserts "
              << std::setw(11) << kdtree_inserts
              << std::setw(11) << kdtreev_inserts
              << std::setw(22) << linear_inserts
              << std::setw(23) << linearv_inserts
              << std::setw(16) << double(kdtree_inserts) / double(kdtreev_inserts) << std::endl;

    // random search points {{{2
    std::vector<Point> searchPoints;
    searchPoints.reserve(NumberOfSearches);
    for (int i = 0; i < NumberOfSearches; ++i) {
        searchPoints.emplace_back(uniform(randomEngine), uniform(randomEngine), uniform(randomEngine));
    }

    tsc.start();  // KdTree searches {{{2
    for (int i = 0; i < NumberOfSearches; ++i) {
        const auto &p = searchPoints[i];
        const auto &p2 = pointsTree.findNearest(p);
        asm(""::"m"(p2));
        //std::cout << "looking near " << p << ", found " << p2 << ", distance " << get_kdtree_distance(p, p2) << '\n';
    }
    tsc.stop();
    const auto time_kdtree = tsc.cycles();

    tsc.start();  // KdTreeV searches {{{2
    for (int i = 0; i < NumberOfSearches; ++i) {
        const auto &p = searchPoints[i];
        const auto &p2 = pointsTreeV.findNearest(p);
        asm(""::"m"(p2));
        //std::cout << "looking near " << p << ", found " << p2 << ", distance " << get_kdtree_distance(p, p2) << '\n';
    }
    tsc.stop();
    const auto time_kdtreev = tsc.cycles();

    tsc.start();  // LinearNeighborSearch searches {{{2
    for (int i = 0; i < NumberOfSearches; ++i) {
        const auto &p = searchPoints[i];
        const auto &p2 = pointsVector.findNearest(p);
        asm(""::"m"(p2));
#if 0
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

    tsc.start();  // LinearNeighborSearchV searches {{{2
    for (int i = 0; i < NumberOfSearches; ++i) {
        const auto &p = searchPoints[i];
        const auto &p2 = linearSearchV.findNearest(p);
        asm(""::"m"(p2));
#if 0
        const auto &p4 = pointsTreeV.findNearest(p);
        if (get_kdtree_distance(p, p2) != get_kdtree_distance(p, p4)) {
            std::cerr << p << " failed " << p2 << " (" << get_kdtree_distance(p, p2)
                      << ") vs. KdTreeV " << p4 << " (" << get_kdtree_distance(p, p4) << ")\n";
        }
#endif
    }
    tsc.stop();
    const auto time_linearv = tsc.cycles();

    // print search timings {{{2
    std::cout << "searches"
              << std::setw(11) << time_kdtree
              << std::setw(11) << time_kdtreev
              << std::setw(22) << time_linear
              << std::setw(23) << time_linearv
              << std::setw(16) << double(time_kdtree) / double(time_kdtreev)
              << std::setw(15) << double(time_linear) / double(time_kdtree)
              << std::setw(16) << double(time_linear) / double(time_linearv) << '\n';

    return 0;
} // }}}1

// vim: foldmethod=marker
