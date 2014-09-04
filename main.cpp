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
#include "point.h"
#include "linearsearch.h"
#include "kdtree.h"

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
