# This code produces all data for visualization later.

import random
import numpy as np
import sys
import multiprocessing
from functools import wraps
data_dir = "../data/"
result_dir = "../results/"
data_k13 = data_dir + "PASHA_k/"
data_w100 = data_dir + "PASHA_w/"
seq_file = "hg38_all.seq"
chmap = {'A': 2, 'C': 0, 'G': 1, 'T': 3}
samples_override = 1000000

def random_mer_generator(w, k):
    modulus = 4 ** k
    cur_kmer = random.randrange(modulus)
    for i in range(w):
        cur_kmer = (cur_kmer * 4 + random.randrange(4)) % modulus
        yield cur_kmer

def prepare_sequence_context_sampler(seq_file):
    with open(seq_file) as f:
        seq = f.readline().strip()
    slen = len(seq)
    def sequence_mer_sampler(w, k):
        start_idx = random.randrange(slen - w - k + 1)
        cseq = seq[start_idx:start_idx + w + k - 1]
        current = 0
        modulus = 4 ** k
        for i in range(k-1):
            current = current * 4 + chmap[cseq[i]]
        for i in range(k-1, w+k-1):
            current = (current * 4 + chmap[cseq[i]]) % modulus
            yield current
    global sampler_func
    sampler_func = sequence_mer_sampler

def const_func(mer):
    return True

def lexi_charge_checker(mers, uhs_func = const_func):
    # Given list of k-mers, check if the context is charged
    v = max(mers) + 1
    for m in mers:
        if uhs_func(m):
            v = min(m, v)

    #  if (v == max(mers) + 1):
        #  with open("results/misses.log", "a") as f:
            #  print(uhs_func, file=f)
            #  print(mers, file=f)
    if v == mers[0]:
        return True
    elif (v == mers[-1]) and (mers.count(v) == 1):
        return True
    else:
        return False

def lexi_wrapper(uhs_func):
    return lambda mer: lexi_charge_checker(mer, uhs_func)

def miniception_uhs(w, k, k0):
    submod = 4 ** k0
    threshold = 14
    if k0 <= threshold:
        min_mapper = list(range(submod))
        random.shuffle(min_mapper)
    else:
        m = 2 ** 67 - 1  # A Mersenne Prime!
        p = random.randrange(m)
        def simple_hash(mer):
            return (mer * p) % m
    assert (k0 < k)
    @wraps(miniception_uhs)
    def func(mer):
        mers = []
        for i in range(k - k0 + 1):
            if k0 <= threshold:
                mers.append(min_mapper[mer % submod])
            else:
                mers.append(simple_hash(mer % submod))
            mer = mer // 4
        return lexi_charge_checker(list(reversed(mers)))
    return func

def external_UHS(k, filename, filename2 = None):
    bpd = 64
    modulus = 4 ** k
    data = [0] * (modulus // bpd)
    def process_line(s):
        assert len(s) == k
        cm = 0
        for c in s:
            cm = cm * 4 + chmap[c]
        data[cm // bpd] |= 2 ** (cm % bpd)
    with open(filename) as f:
        for line in f:
            s = line.strip()
            process_line(s)
    if filename2 is not None:
        with open(filename2) as f:
            for line in f:
                s = line.strip()
                process_line(s)
    @wraps(external_UHS)
    def func(mer):
        return (data[mer // bpd] & (2 ** (mer % bpd))) > 0
    return func

def calculate_energy(mers, uhs_func):
    dv = 0
    seen_mers = set()
    if uhs_func(mers[0]):
        dv += 1
    for mer in mers[:-1]:
        if uhs_func(mer):
            seen_mers.add(mer)
    mer = mers[-1]
    if uhs_func(mer):
        if mer not in seen_mers:
            dv += 1
        seen_mers.add(mer)
    if len(seen_mers) == 0:
        return 0
    else:
        return dv / len(seen_mers)

def energy_wrapper(func):
    return lambda mers: calculate_energy(mers, func)

def calculate_density(w, k, cf, samples):
    count = 0
    for i in range(samples):
        mers = []
        if specific_densities:
            for mer in sampler_func(w+1, k):
                mers.append(mer)
        else:
            for mer in random_mer_generator(w+1, k):
                mers.append(mer)
        count += cf(mers)
    return count / samples

def mc_param_search(w, k, k0_low, k0_high, samples):
    best_k0 = 0
    best_density = 1
    for k0 in range(k0_low, k0_high + 1):
        print("Miniception: w={}, k={}, k0={}".format(w, k, k0))
        cf = energy_wrapper(miniception_uhs(w, k, k0))
        density = calculate_density(w, k, cf, samples)
        if (density < best_density):
            best_k0 = k0
            best_density = density
    return (best_k0, best_density)


def k13_subprocess(v):
    samples = samples_override
    spec_file = data_k13 + "PASHA13_{}.txt"
    decycle_file = data_k13 + "decyc13.txt"
    print("w =", v, " - Lexi")
    lexi_density = calculate_density(v, 13, lexi_wrapper(const_func), samples)
    print("w =", v, " - Random")
    random_density = calculate_density(v, 13, energy_wrapper(const_func), samples)
    mc_choice, mc_density = mc_param_search(v, 13, 3, 6, samples)
    print("w =", v, " - External Parse")
    pasha_uhs = external_UHS(13, spec_file.format(v), decycle_file)
    print("w =", v, " - External Lexi")
    pasha_lex_density = calculate_density(v, 13, lexi_wrapper(pasha_uhs), samples)
    print("w =", v, " - External Random")
    pasha_density = calculate_density(v, 13, energy_wrapper(pasha_uhs), samples)
    return (v, lexi_density, random_density, mc_density, pasha_lex_density, pasha_density, mc_choice)

def calculate_fixed_k13():
    values = list(range(20, 210, 10))
    p = multiprocessing.Pool(processes=20)
    returns = p.map(k13_subprocess, values)
    with open(result_dir + "/result_k13" + ("_seq" if specific_densities else "") + ".txt", "w") as output_file:
        for data in returns:
            print(*data, sep=',', file=output_file)

def w100_subprocess(v):
    samples = samples_override
    spec_file = data_w100 + "PASHA{}_100.txt"
    decycle_file = data_w100 + "decyc{}.txt"
    print("k =", v, " - Lexi")
    lexi_density = calculate_density(100, v, lexi_wrapper(const_func), samples)
    print("k =", v, " - Random")
    random_density = calculate_density(100, v, energy_wrapper(const_func), samples)
    mc_choice, mc_density = mc_param_search(100, v, 2, min(v // 2, 6), samples)
    print("k =", v, " - External Parse")
    pasha_uhs = external_UHS(v, spec_file.format(v), decycle_file.format(v))
    print("k =", v, " - External Lexi")
    pasha_lex_density = calculate_density(100, v, lexi_wrapper(pasha_uhs), samples)
    print("k =", v, " - External Random")
    pasha_density = calculate_density(100, v, energy_wrapper(pasha_uhs), samples)
    return (v, lexi_density, random_density, mc_density, pasha_lex_density, pasha_density, mc_choice)

def calculate_fixed_w100():
    values = [7,8,9,10,12,13,14,15]
    p = multiprocessing.Pool(processes=20)
    returns = p.map(w100_subprocess, values)
    with open(result_dir + "/result_w100" + ("_seq" if specific_densities else "") + ".txt", "w") as output_file:
        for data in returns:
            print(*data, sep=',', file=output_file)

def free_subprocess(param):
    samples = samples_override
    w, k, k0_low, k0_high = param
    print("w = {}, k = {} - Lexi".format(w, k))
    lexi_density = calculate_density(w, k, lexi_wrapper(const_func), samples)
    print("w = {}, k = {} - Random".format(w, k))
    random_density = calculate_density(w, k, energy_wrapper(const_func), samples)
    mc_choice, mc_density = mc_param_search(w, k, k0_low, k0_high, samples)
    return (w, k, lexi_density, random_density, mc_density, mc_choice)

def calculate_sample_free():
    # small fixed w = 10, large fixed w = 100. k = 10 to 30.
    # k' = 3 to 7
    params = []
    for k in range(30, 9, -1):
        if k < 15:
            params.append((10, k, max(1,k-10), 8))
        else:
            params.append((10, k, k-10, k-10))
        params.append((100, k, 3, 7))
    p = multiprocessing.Pool(processes = 20)
    returns = p.map(free_subprocess, params)
    with open(result_dir + "/free_result" + ("_seq" if specific_densities else "") + ".txt", "w") as output_file:
        for data in returns:
            print(*data, sep=',', file=output_file)

if __name__ == "__main__":
    global specific_densities
    if len(sys.argv) == 2:
        specific_densities = True
    else:
        specific_densities = False
    if specific_densities:
        print("Parsing sequence now")
        prepare_sequence_context_sampler(seq_file)
    calculate_fixed_k13()
    calculate_fixed_w100()
    calculate_sample_free()
