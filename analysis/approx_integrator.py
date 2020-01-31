'''
This is the code to implement approximate integration, to plot D(x) for x>2
via dynamic programming. Use integrator_plot.py to read the data and plot.
'''

import numpy as np
import pickle
import matplotlib.pyplot as plt
from tqdm import trange
results_dir = "../results/"
w = 2500
maxc = 8
maxn = 20

s = w * maxc
d = np.zeros(shape=(maxn+1, s), dtype=np.float)
dp = np.zeros(shape=(maxn+1, s), dtype=np.float)
dm = np.zeros(shape=(maxn+1, s), dtype=np.float)
# get d[0]
for i in range(s):
    if i < w:
        d[0][i] = 1
    else:
        f = i - 2 * (i-(w-1))
        d[0][i] = max(0, f/i)
for m in range(1, maxn+1):
    print("m=", m)
    for i in trange(w, s):  # i slots
        total = 0
        for j in range(i):  # min element slot
            seg_left = j
            seg_right = i - j - 1
            conv_idx = m
            if seg_left >= w-1:
                conv_idx -= 1
            if seg_right >= w-1:
                conv_idx -= 1
            # do the convolution now
            for k in range(conv_idx + 1):
                total += d[k][seg_left] * d[conv_idx - k][seg_right]
        d[m][i] = total / i

for m in range(1, maxn+1):
    print("[Q+] m=", m)
    for i in trange(w, s):
        total = w * d[m-1][i-1]
        for j in range(w, i):
            seg_left = j
            seg_right = i - j - 1
            conv_idx = m - 1
            if seg_right >= w-1:
                conv_idx -= 1
            for k in range(conv_idx + 1):
                total += dp[k][seg_left] * d[conv_idx - k][seg_right]
        dp[m][i] = total / i

for m in range(1, maxn+1):
    print("[Q-] m=", m)
    for i in trange(w, s):
        if i >= 2*w-1:
            total = w * d[m-2][i-w]
        else:
            total = w * d[m-1][i-w]
        for j in range(w, i):
            seg_left = j
            seg_right = i - j - 1
            conv_idx = m - 1
            if seg_right >= w-1:
                conv_idx -= 1
            for k in range(conv_idx + 1):
                total += dm[k][seg_left] * d[conv_idx - k][seg_right]
        dm[m][i] = total / i

#  sanity check
for i in range(maxn):
    print(i, d[i][2*w])

densities = np.zeros(shape=s)
xs = []
ys = []
zs = []
# density bound calculation
for i in range(2 * w, s):
    rp = 1
    total = 0
    for j in range(1, maxn + 1):
        v = dp[j][i]
        rp -= v
        total += v/j
    assert (rp > -0.01)
    total += rp / (maxn + 1)
    rp = 1
    for j in range(1, maxn + 1):
        v = dm[j][i]
        rp -= v
        total += v / j
    assert (rp > -0.01)
    total += rp / (maxn + 1)
    cc = i / w - 1
    densities[i] = total * 2 * cc
    xs.append(cc)
    ys.append(densities[i])
    #  zs.append(1.5 * cc /(cc + 1))

plt.plot(xs, ys)
#  plt.plot(xs, zs)
plt.hlines(2, min(xs), max(xs))
#  plt.savefig("temp.jpg")

with open(results_dir + "density_dump", "bw") as f:
    pickle.dump(densities, f)

#  with open("dp_dump", "bw") as f:
    #  pickle.dump(d, f)
