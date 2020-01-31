# This is a reference implementation of the Miniception, achieving linear
# overhead per character in input sequence assuming kmers fit into a word.
from collections import deque
import random
chmap = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
# Specify the parameters of the Miniception here.
w = 30
k = 30
k0 = 10
w0 = k - k0


class MonotoneQueue:
    def __init__(self, expire_time):
        self.q = deque()
        self.expire_time = expire_time

    def expire(self, time):
        # from left (front): remove expired items
        while len(self.q) > 0:
            ctime, citem = self.q[0]
            if ctime <= time - self.expire_time:
                self.q.popleft()
            else:
                return

    def insert(self, time, item):
        self.expire(time)
        # from right (rear): remove dominated items
        while len(self.q) > 0:
            ctime, citem = self.q[-1]
            # We use strict inequality here, as we allow items with same
            # value in the queue as the minimizer break ties by favoring
            # leftmost kmer.
            if item < citem:
                self.q.pop()
            else:
                break
        self.q.append((time, item))

    def peek(self):
        return self.q[0]

def seed_order_function(mer):
    '''
    This is the ordering function for the seed minimizer.
    It is supposed to be a function that takes an integer in [4^k_0] as
    input, and output something comparable.
    By returning the input as untouched, we implement the lexicographical order.
    Note:
    If you are implementing a hashing function, the behavior of the miniception
    process is undefined if collisions happen. In that case, one band-aid fix
    is to return (hash_value, input_value) so it defaults to lexi-order in
    case of a collision.
    '''
    return mer

def order_function(mer):
    '''
    This is the ordering function for the miniception.
    See comments above for description of this function.
    '''
    return mer

def miniception(s):
    '''
    This is the main Miniception function.
    The function will yield i for every selected k-mer starting at i
    (all indices are 0-based).
    '''
    # We requre the sequence to contain at least one window.
    if len(s) <= w + k - 1:
        return
    # Warm-up phase: parse the first window, and output first selected
    # location.
    last_yield = -1
    seed_modulus = 4 ** k0
    main_modulus = 4 ** k
    seed_queue = MonotoneQueue(w0 + 1)
    main_queue = MonotoneQueue(w)
    small_mer = main_mer = 0
    for i in range(len(s)):
        small_mer = (small_mer * 4 + chmap[s[i]]) % seed_modulus
        main_mer = (main_mer * 4 + chmap[s[i]]) % main_modulus
        # Runs the seed queue after first complete k0-mer observed
        if i >= k0 - 1:
            seed_queue.insert(i, seed_order_function(small_mer))
        # Runs the main queue after first complete k-mer observed
        if i >= k - 1:
            seed_pick = seed_queue.peek()[0]
            # This checks if the k-mer is a charged context of the seed
            # minimizer by retriving the minimizer pick.
            if (seed_pick == i) or (seed_pick == i - w0):
                main_queue.insert(i, order_function(main_mer))
            else:
                # We need to clean expired elements in the main monotone queue
                main_queue.expire(i)
        # Runs the acutal miniception after first complete window observed
        if i >= w + k - 2:
            # Get the current minimizer picks and compare against last pick
            pick_loc = main_queue.peek()[0]
            if pick_loc != last_yield:
                # The pick_loc is the index of the last character in the picked
                # k-mer. The starting location shift by (k-1) accordingly
                yield pick_loc - (k - 1)
                last_yield = pick_loc

