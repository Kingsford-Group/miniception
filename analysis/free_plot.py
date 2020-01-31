import matplotlib
matplotlib.rcParams.update({'font.size':6})
import matplotlib.pyplot as plt
plt.figure(figsize=(3.5, 2))
import csv
import sys
specific_densities = len(sys.argv) == 2
sid = "_seq" if specific_densities else ""
f = open("../results/free_result" + sid + ".txt")
dat = csv.reader(f)
columns = ["", "", "Lexicographic", "Random", "Miniception"]
datas_10 = [[], [], [], [], []]
datas_100 = [[], [], [], [], []]
for row in dat:
    w = int(row[0])
    k = int(row[1])
    if w == 100:
        datas_100[0].append(k)
        for i in range(2, 5):
            datas_100[i].append(float(row[i]) * 101)
    else:
        datas_10[0].append(k)
        for i in range(2, 5):
            datas_10[i].append(float(row[i]) * 11)
plt.subplot(121)
for i in range(2, 5):
    plt.plot(datas_100[0], datas_100[i], label=columns[i])
plt.xlabel("Value of k (w = 100)")
plt.ylabel("Density Factor")
#  plt.savefig("results/free_100" + sid + ".pdf")
#  plt.clf()
plt.subplot(122)
for i in range(2, 5):
    plt.plot(datas_10[0], datas_10[i], label=columns[i])
plt.xlabel("Value of k (w = 10)")
plt.ylabel("Density Factor")
plt.legend(loc='lower left', bbox_to_anchor=(-1.32, -0.4), ncol=3)
plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig("../results/free" + sid + ".pdf")
