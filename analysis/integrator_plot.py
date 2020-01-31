import pickle
import matplotlib
matplotlib.rcParams.update({'font.size':8})
import matplotlib.pyplot as plt
with open("../results/density_dump", "rb") as f:
    ds = pickle.load(f)

w = 2500
s = len(ds)
xs = []
ys = []
for i in range(2 * w, s):
    xs.append(i / w - 1)
    ys.append(ds[i])

plt.figure(figsize=(3.5, 2.8))
plt.plot(xs, ys, label="Miniception")
plt.hlines(2, min(xs), max(xs), label="Random Bound")
plt.xlabel("Window length / k-mer length")
plt.ylabel("Asymptotic Density Factor")
plt.legend(loc="lower right")
plt.tight_layout(pad=0.5)
plt.savefig("../results/asympto.pdf")
