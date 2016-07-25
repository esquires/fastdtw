#!python
#cython: boundscheck=False, cdivision=True, wraparound=False

from __future__ import absolute_import, division

from pprint import pprint

from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.utility cimport pair
from cython.operator cimport dereference

cdef struct CostEntry:
    double cost
    int prev_x_idx
    int prev_y_idx


def fastdtw(x, y, radius=1, dist=lambda a, b: abs(a - b)):
    min_time_size = radius + 2

    if len(x) < min_time_size or len(y) < min_time_size:
        return dtw(x, y, dist=dist)

    x_shrinked = __reduce_by_half(x)
    y_shrinked = __reduce_by_half(y)
    distance, path = fastdtw(x_shrinked, y_shrinked, radius=radius, dist=dist)
    window = __expand_window(path, len(x), len(y), radius)
    return dtw(x, y, window, dist=dist)


def dtw(x, y, window=None, dist=lambda a, b: abs(a - b)):

    cdef int len_x, len_y, i, j
    cdef double inf, cost_up, cost_left, cost_cnr, dt
    len_x, len_y = len(x), len(y)
    if window is None:
        window = [(i, j) for i in range(len_x) for j in range(len_y)]

    # window could be a vector but profiling indicated little improvement
    window = [(i + 1, j + 1) for i, j in window]

    cdef map[pair[int, int], CostEntry] D
    cdef CostEntry cost_entry

    cost_entry.cost = 0
    cost_entry.prev_x_idx = 0
    cost_entry.prev_y_idx = 0
    inf = float('inf')

    D[0, 0] = cost_entry

    for i, j in window:

        dt = dist(x[i-1], y[j-1])

        # note: cython will converting tuples to a c++ pair. Profiling
        # indicates this makes about a 10% performance so it was not performed
        cost_up = D[i-1, j].cost if D.find((i-1, j)) != D.end() else inf
        cost_left = D[i, j-1].cost if D.find((i, j-1)) != D.end() else inf
        cost_cnr = D[i-1, j-1].cost if D.find((i-1, j-1)) != D.end() else inf

        if cost_up < cost_left and cost_up < cost_cnr:
            cost_entry.cost = cost_up + dt
            cost_entry.prev_x_idx = i - 1
            cost_entry.prev_y_idx = j

        elif cost_left < cost_cnr:
            cost_entry.cost = cost_left + dt
            cost_entry.prev_x_idx = i
            cost_entry.prev_y_idx = j - 1

        else:
            cost_entry.cost = cost_cnr + dt
            cost_entry.prev_x_idx = i - 1
            cost_entry.prev_y_idx = j - 1

        D[i, j] = cost_entry

    path = []
    i, j = len_x, len_y
    while not (i == j == 0):
        path.append((i-1, j-1))
        cost_entry = D[i, j]
        i = cost_entry.prev_x_idx
        j = cost_entry.prev_y_idx
    path.reverse()
    return (D[len_x, len_y].cost, path)


cdef __reduce_by_half(x):
    cdef int i
    return [(x[i] + x[1+i]) / 2.0 for i in range(0, len(x) - len(x) % 2, 2)]


cdef __expand_window(path, int len_x, int len_y, int radius):
    cdef set[pair[int, int]] path_, window_
    cdef pair[int, int] key
    cdef int i, j, a, b, start_j, new_start_j

    # this could be made into a vector. Profiling indicated little performance
    # improvement
    lst = [(a, b)
           for a in range(-radius, radius+1)
           for b in range(-radius, radius+1)]

    for i, j in path:
        for a, b in lst:
            key.first = i + a
            key.second = j + b
            path_.insert(key)

    for i, j in path_:
        for a, b in ((i * 2, j * 2), (i * 2, j * 2 + 1),
                     (i * 2 + 1, j * 2), (i * 2 + 1, j * 2 + 1)):
            key.first = a
            key.second = b
            window_.insert(key)

    window = []
    start_j = 0
    for i in range(0, len_x):
        new_start_j = -1
        for j in range(start_j, len_y):
            # we could use a key rather than a tuple to do the find lookup.
            # Profiling indicated little performance improvement
            if window_.find((i, j)) != window_.end():
                window.append((i, j))
                if new_start_j == -1:
                    new_start_j = j
            elif new_start_j != -1:
                break
        start_j = new_start_j

    return window
