#import subprocess
#subprocess.call(["./clean_data.sh"])

import collections

mut_list = []
with open("./data/ssm.open.SKCM-US_intersect.bed") as f:
    for line in f:
        _,_,_,_,pos,tf = line.strip().split()
        mut_list.append((int(pos), tf))

pos_dict = collections.defaultdict(int)
for pos, _ in mut_list:
    pos_dict[pos] += 1

pos_dict_tf = {}
for pos, tf in mut_list:
    if tf not in pos_dict_tf:
        pos_dict_tf[tf] = collections.defaultdict(int)
    pos_dict_tf[tf][pos] += 1

from matplotlib import pyplot as plt

X = sorted([int(key) for key in pos_dict.keys()])
y = [pos_dict[key] for key in X]

plt.plot(X, y)
plt.xlim(min(X), max(X))
plt.xlabel("Distance from TFBS center (bp)")
plt.ylabel("Number of mutations")
plt.title("Mutation profile for all TFs")
plt.savefig("mut-profile.png")
