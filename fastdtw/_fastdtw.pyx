#!python
#cython: boundscheck=False, cdivision=True, wraparound=False

from __future__ import absolute_import, division
import numbers
import numpy as np

from libc.math cimport INFINITY, pow, fabs
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


def fastdtw(x, y, int radius=1, dist=None):
    ''' return the approximate distance between 2 time series with O(N)
        time and memory complexity

        Parameters
        ----------
        x : array_like
            input array 1
        y : array_like
            input array 2
        radius : int
            size of neighborhood when expanding the path. A higher value will
            increase the accuracy of the calculation but also increase time
            and memory complexity. A radius equal to the size of x and y will
            yield an exact dynamic time warping calculation.
        dist : function or int
            The method for calculating the distance between x[i] and y[j]. If
            dist is an int of value p > 0, then the p-norm will be used. If
            dist is a function then dist(x[i], y[j]) will be used. If dist is
            None then abs(x[i] - y[j]) will be used.

        Returns
        -------
        float
            the approximate distance between the 2 time series

        >>> import numpy as np
        >>> import fastdtw
        >>> x = np.array([1, 2, 3, 4, 5.0])
        >>> y = np.array([2, 3, 4.0])
        >>> fastdtw.fastdtw(x, y)
        (2.0, [(0, 0), (1, 0), (2, 1), (3, 2), (4, 2)])

        Optimization Considerations
        ---------------------------
        This function runs fastest if the following conditions are satisfied:
            1) x and y are either 1 or 2d numpy arrays whose dtype is a
               subtype of np.float
            2) The dist input is a positive integer or None
    '''

    if isinstance(dist, numbers.Number) and dist <= 0:
        raise ValueError('dist cannot be a negative integer')
    elif (isinstance(x, np.ndarray) and
          isinstance(y, np.ndarray) and
          x.ndim > 1 and y.ndim > 1 and
          x.shape[1] != y.shape[1]):

        raise ValueError('x and y must have the same width')

    return __fastdtw(x, y, radius, dist)


cdef __fastdtw(x, y, int radius, dist):
    cdef int min_time_size

    min_time_size = radius + 2

    if len(x) < min_time_size or len(y) < min_time_size:
        return dtw(x, y, dist=dist)

    x_shrinked = __reduce_by_half(x)
    y_shrinked = __reduce_by_half(y)
    distance, path = fastdtw(x_shrinked, y_shrinked, radius, dist)

    window = __expand_window(path, len(x), len(y), radius)
    return dtw(x, y, window, dist)


cdef inline double manhattan(double a, double b):
    return fabs(a - b)


def __create_linalg_norm(p):
    def f(a, b):
        return np.linalg.norm(a - b, p)

    return f


def dtw(x, y, window=None, dist=lambda a, b: abs(a - b)):
    cdef unsigned int len_x, len_y, i, j, k, len_window

    cdef double cost_up, cost_left, cost_cnr, dt, diff, sm
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

    # initializing to avoid compiler warnings although if these variables are
    # used they will always be set below
    cdef unsigned int use_1d = 0, use_2d = 0, width = 0

    cdef double[:] x_arr1d, y_arr1d
    cdef double[:, :] x_arr2d, y_arr2d

    # define the dist function
    # if it is numpy array we can get an order of magnitude improvement
    # by calculating the p-norm ourselves.
    cdef double pnorm = -1.0 if not isinstance(dist, numbers.Number) else dist
    if ((dist is None or pnorm > 0) and
        (isinstance(x, np.ndarray) and isinstance(y, np.ndarray) and
         np.issubdtype(x.dtype, np.float) and
         np.issubdtype(y.dtype, np.float))):

        if x.ndim == 1:
            use_1d = 1
            x_arr1d = x
            y_arr1d = y
        elif x.ndim == 2:
            use_2d = 1
            x_arr2d = x
            y_arr2d = y
            width = x.shape[1]

    if not use_1d and not use_2d:
        if dist is None:
            dist = manhattan
        elif pnorm > 0:
            dist = __create_linalg_norm(pnorm)

    cdef map[pair[int, int], CostEntry] D
    cdef CostEntry cost_entry

    cost_entry.cost = 0
    cost_entry.prev_x_idx = 0
    cost_entry.prev_y_idx = 0

    D[__get_key(0, 0)] = cost_entry

    for k in range(len_window):
        i = window_vec[k].first
        j = window_vec[k].second

        if use_1d:
            dt = abs(x_arr1d[i - 1] - y_arr1d[j - 1])
        elif use_2d:
            sm = 0
            for k in range(width):
                diff = abs(x_arr2d[i - 1, k] - y_arr2d[j - 1, k])
                sm += pow(diff, pnorm)
            dt = pow(sm, 1 / pnorm)
        else:
            dt = dist(x[i - 1], y[j - 1])

        cost_up = __get_cost(D, i - 1, j, INFINITY)
        cost_left = __get_cost(D, i, j - 1, INFINITY)
        cost_cnr = __get_cost(D, i - 1, j - 1, INFINITY)

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
    ''' return x whose size is floor(len(x) / 2) with linear interpolation
        between values
    '''
    cdef int i, mx
    try:
        mx = x.shape[0] - x.shape[0] % 2
        return (x[:mx:2] + x[1:mx:2]) / 2
    except AttributeError:
        return [(x[i] + x[1+i]) / 2 for i in range(0, len(x) - len(x) % 2, 2)]


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
