#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division

from collections import defaultdict

from six.moves import xrange

def fastdtw(x, y, int radius=1, dist=lambda a, b: abs(a - b)):
    cdef int min_time_size
    min_time_size = radius + 2

    if len(x) < min_time_size or len(y) < min_time_size:
        return dtw(x, y, dist=dist)

    x_shrinked = __reduce_by_half(x)
    y_shrinked = __reduce_by_half(y)
    distance, path = fastdtw(x_shrinked, y_shrinked, radius=radius, dist=dist)
    window = __expand_window(path, len(x), len(y), radius)
    return dtw(x, y, window, dist=dist)


def dtw(x, y, window=None, dist=lambda a, b: abs(a - b)):
    cdef int len_x, len_y, len_y_p1, ip, jp, i, j
    cdef double dt, d_00, d_10, d_01, inf

    len_x, len_y = len(x), len(y)
    if window is None:
        window = [(i, j) for i in xrange(len_x) for j in xrange(len_y)]

    inf = float('inf')
    D = defaultdict(lambda: (inf,))
    D[0, 0] = (0, 0, 0)
    len_y_p1 = len_y + 1
    for i, j in window:
        dt = dist(x[i], y[j])

        ip = i + 1
        jp = j + 1
        d_00 = D[i, j][0]
        d_10 = D[ip, j][0]

        if j == len_y_p1:
            if d_00 < d_10:
                D[ip, jp] = (d_00 + dt, i, j)
            else: 
                D[ip, jp] = (d_10 + dt, i, jp)
        else: 
            d_01 = D[i, jp][0]

            if d_00 < d_10 and d_00 < d_01:
                D[ip, jp] = (d_00 + dt, i, j)
            elif d_10 < d_01:
                D[ip, jp] = (d_10 + dt, i + 1, j)
            else:
                D[ip, jp] = (d_01 + dt, i, j + 1)
 

    path = []
    i, j = len_x, len_y
    while not (i == j == 0):
        path.append((i-1, j-1))
        i, j = D[i, j][1], D[i, j][2]
    path.reverse()
    return (D[len_x, len_y][0], path)


cdef __reduce_by_half(x):
    cdef int i
    return [(x[i] + x[1+i]) / 2 for i in xrange(0, len(x) - len(x) % 2, 2)]


cdef __expand_window(path, int len_x, int len_y, int radius):
    cdef int i, j, a, b, i2, j2, i2p1, j2p1, start_j, new_start_j
    path_ = set(path)

    lst = xrange(-radius, radius+1)
    lst = [(a, b) for a in lst for b in lst]

    for i, j in path:
        for a, b in lst:
            path_.add((i+a, j+b))

    window_ = set()
    for i, j in path_:
        i2 = i * 2
        j2 = j * 2
        i2p1 = i2 + 1
        j2p1 = j2 + 1
        window_.add((i2, j2))
        window_.add((i2, j2p1))
        window_.add((i2p1, j2))
        window_.add((i2p1, j2p1))


    window = []
    start_j = 0
    for i in xrange(0, len_x):
        new_start_j = -1
        for j in xrange(start_j, len_y):
            if (i, j) in window_:
                window.append((i, j))
                if new_start_j == -1:
                    new_start_j = j
            elif new_start_j != -1:
                break
        start_j = new_start_j

    return window

