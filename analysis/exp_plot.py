import matplotlib
matplotlib.rcParams.update({'font.size':6})
import matplotlib.pyplot as plt
import csv
import sys

specific_densities = len(sys.argv) == 2
sid = "_seq" if specific_densities else ""

plt.figure(figsize=(3.5, 2))
plt.subplot(121)

f = open("../results/result_k13" + sid + ".txt")
dat = csv.reader(f)
columns = ["", "Lexicographic", "Random", "Miniception", "PASHA", "PASHA-Random"]
datas = [[], [], [], [], [], []]
for row in dat:
    w = int(row[0])
    datas[0].append(w)
    for i in range(1, 6):
        datas[i].append(float(row[i]) * (w+1))

for i in range(1, 6):
    if (i == 5) and not specific_densities:
        continue
    plt.plot(datas[0], datas[i], label=columns[i])

#  plt.legend()
plt.xlabel("Value of w (k = 13)")
plt.ylabel("Density Factor")
#  plt.tight_layout(pad=0.5)
#  plt.savefig("results/pasha_k" + sid + ".pdf")
#  plt.clf()

f.close()
plt.subplot(122)
f = open("../results/result_w100" + sid + ".txt")
dat = csv.reader(f)
datas = [[], [], [], [], [], []]
for row in dat:
    w = int(row[0])
    datas[0].append(w)
    for i in range(1, 6):
        datas[i].append(float(row[i]) * 101)

for i in range(1, 6):
    if (i == 5) and not specific_densities:
        continue
    plt.plot(datas[0], datas[i], label=columns[i])
if specific_densities:
    plt.legend(loc='lower left', bbox_to_anchor=(-1.45, -0.52), ncol=3)
else:
    plt.legend(loc='lower left', bbox_to_anchor=(-1.70, -0.44), ncol=4)
plt.xlabel("Value of k (w = 100)")
plt.ylabel("Density Factor")
plt.tight_layout(rect=[0, 0.1, 1, 1])
plt.savefig("../results/pasha" + sid + ".pdf")
f.close()
