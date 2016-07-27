#!python
#cython: boundscheck=False, cdivision=True, wraparound=False

from __future__ import absolute_import, division

from libcpp.vector cimport vector  # noqa
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.utility cimport pair
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement

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
    cdef unsigned int len_x, len_y, i, j, k, len_window

    cdef double inf, cost_up, cost_left, cost_cnr, dt
    cdef pair[int, int] key_up, key_left, key_cnr, key_next
    cdef vector[pair[int, int]] window_vec
    len_x, len_y = len(x), len(y)
    if window is None:
        window = [(i, j) for i in range(len_x) for j in range(len_y)]

    len_window = len(window)
    window_vec.resize(len_window)
    for i in range(len_window):
        window_vec[i].first = window[i][0] + 1
        window_vec[i].second += window[i][1] + 1

    cdef map[pair[int, int], CostEntry] D
    cdef CostEntry cost_entry

    cost_entry.cost = 0
    cost_entry.prev_x_idx = 0
    cost_entry.prev_y_idx = 0
    inf = float('inf')

    D[__get_key(0, 0)] = cost_entry

    for k in range(len_window):
        i = window_vec[k].first
        j = window_vec[k].second

        dt = dist(x[i-1], y[j-1])

        cost_up = __get_cost(D, i - 1, j, inf)
        cost_left = __get_cost(D, i, j - 1, inf)
        cost_cnr = __get_cost(D, i - 1, j - 1, inf)

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

        D[__get_key(i, j)] = cost_entry

    path = []
    i, j = len_x, len_y
    while not (i == j == 0):
        path.append((i-1, j-1))
        cost_entry = D[__get_key(i, j)]
        i = cost_entry.prev_x_idx
        j = cost_entry.prev_y_idx
    path.reverse()
    return (D[len_x, len_y].cost, path)


cdef __reduce_by_half(x):
    cdef int i, len_x = len(x)
    out = []
    for i in range(0, len_x - len_x % 2, 2):
        out.append((x[i] + x[1+i]) / 2.0)
    return out


cdef __expand_window(path, int len_x, int len_y, int radius):
    cdef set[pair[int, int]] path_, window_
    cdef vector[pair[int, int]] vec
    cdef pair[int, int] key, val
    cdef int i, j, k, a, b, start_j, new_start_j, len_vec

    len_vec = (2 * radius + 1)**2
    vec.reserve(len_vec)
    i = 0
    for a in range(-radius, radius + 1):
        for b in range(-radius, radius + 1):
            vec[i] = __get_key(a, b)
            i += 1

    for i, j in path:
        for k in range(len_vec):
            key.first = i + vec[k].first
            key.second = j + vec[k].second
            path_.insert(key)

    it = path_.begin()
    while it != path_.end():
        i = deref(it).first
        j = deref(it).second
        preincrement(it)

        window_.insert(__get_key(i * 2, j * 2))
        window_.insert(__get_key(i * 2, j * 2 + 1))
        window_.insert(__get_key(i * 2 + 1, j * 2))
        window_.insert(__get_key(i * 2 + 1, j * 2 + 1))

    window = []
    start_j = 0
    for i in range(0, len_x):
        new_start_j = -1
        for j in range(start_j, len_y):
            if window_.find(__get_key(i, j)) != window_.end():
                window.append((i, j))
                if new_start_j == -1:
                    new_start_j = j
            elif new_start_j != -1:
                break
        start_j = new_start_j

    return window

cdef pair[int, int] __get_key(int i, int j):
    cdef pair[int, int] out
    out.first = i
    out.second = j
    return out

cdef double __get_cost(map[pair[int, int], CostEntry] &D,  # noqa
                       int i, int j, double inf):
    cdef pair[int, int] key = __get_key(i, j)
    cdef double cost

    it = D.find(key)
    if it != D.end():
        cost = deref(it).second.cost
    else:
        cost = inf
    return cost
