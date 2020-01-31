import time
import random
from reference_impl import *
def naive_miniception(s):
    '''
    This is a naive implementation of the Miniception, intended to check
    the correctness of the actual miniception algorithm.
    '''
    def to_num(s):
        ret = 0
        for c in s:
            ret = ret * 4 + chmap[c]
        return ret
    def get_minimizer_loc(x, w, k, valid_func, order_func):
        assert len(x) == w + k - 1
        ret = -1
        for i in range(w):
            if valid_func(x[i:i+k]):
                if ret == -1:
                    ret = i
                elif order_func(to_num(x[i:i+k])) < order_func(to_num(x[ret:ret+k])):
                    ret = i
        return ret
    def uhs_checker(s):
        assert len(s) == k
        p = get_minimizer_loc(s, w0 + 1, k0, lambda _ : True, seed_order_function)
        return (p == 0) or (p == w0)

    if len(s) < w + k - 1:
        return
    last_pick = -1
    for i in range(len(s) - (w + k - 1) + 1):  # last i is len(s) - (w+k-1)
        current_window = s[i:i+(w+k-1)]
        current_pick = i + get_minimizer_loc(current_window, w, k, uhs_checker, order_function)
        if current_pick != last_pick:
            yield current_pick
            last_pick = current_pick


if __name__ == "__main__":
    s = ""
    for i in range(10000):
        s += random.choice('ACTG')
    print(s)
    l1 = []
    l2 = []
    c = time.perf_counter()
    for val in miniception(s):
        l1.append(val)
    print("Original Miniception finishes in {:.2f} s".format(time.perf_counter() - c))
    c = time.perf_counter()
    for val in naive_miniception(s):
        l2.append(val)
    print("Naive Miniception finishes in {:.2f} s".format(time.perf_counter() - c))
    print(l1)
    print(len(l1) / len(s) * (w+1))
    print(l2)
    print(len(l2) / len(s) * (w+1))
    assert len(l1) == len(l2)
    for i in range(len(l1)):
        assert l1[i] == l2[i]
