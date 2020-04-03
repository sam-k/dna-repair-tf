#!/usr/bin/env python

import matplotlib, os

os.environ["QT_QPA_PLATFORM"] = "offscreen"
matplotlib.use("agg")

import collections, sys
from matplotlib import pyplot as plt
from matplotlib import gridspec as gs


### Set global variables ###


WORKSPACE = "/data/gordanlab/samkim/dna-repair-tf"
RUN_ID = sys.argv[1]  # run ID
WHICH_DATA = sys.argv[2]  # data group name
SUFFIX = "" if len(sys.argv) <= 3 else "_" + sys.argv[3]
# suffix to filename (currently only for mut-profile_merged.cl.sh)


### Get mutation counts per position, and per position per TF ###


def get_dists(mut_dataset_name):
    mut_list = []
    with open(
        "{}/data/ssm.open.{}_{}{}_centered.bed".format(
            WORKSPACE, RUN_ID, mut_dataset_name, SUFFIX
        )
    ) as f:
        for line in f:
            _, dist, _, mut, tf = line.strip().split()
            mut_list.append((int(dist), mut, tf))

    counts = collections.defaultdict(int)
    for dist, _, _ in mut_list:
        counts[dist] += 1

    counts_by_tf = {}
    for dist, mut, tf in mut_list:
        if tf not in counts_by_tf:
            counts_by_tf[tf] = collections.defaultdict(int)  # initialize
        counts_by_tf[tf][dist] += 1

    return counts, counts_by_tf


### Create mutation profile plots ###


def plot_dists(
    counts,
    counts_by_tf,
    mut_dataset_name=None,
    figsize=(10, 14),
    h1=0.2,
    h2=0.6,
    w2=0.5,
):
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(hspace=h1)
    outer = gs.GridSpec(2, 1, height_ratios=[2, 5])
    gs1 = gs.GridSpecFromSubplotSpec(1, 6, subplot_spec=outer[0])
    gs2 = gs.GridSpecFromSubplotSpec(5, 6, subplot_spec=outer[1], hspace=h2, wspace=w2)

    X1 = sorted([int(dist) for dist in counts.keys()])
    y1 = [counts[dist] for dist in X1]

    ax1 = plt.subplot(gs1[:, 1:-1])
    ax1.plot(X1, y1)
    ax1.set_xlim(-1000, 1000)
    ax1.set_ylim(0, None)
    ax1.set_xlabel("Distance from TFBS center (bp)")
    ax1.set_ylabel("Number of mutations")
    if mut_dataset_name is None:
        ax1.set_title("Mutation profile for all TFs")
    else:
        ax1.set_title("{} mutation profile for all TFs".format(mut_dataset_name))

    ordered_counts = sorted(
        counts_by_tf.items(), key=lambda t: sum(t[1].values()), reverse=True
    )
    row, col = 0, 0
    for tf, counts in ordered_counts:
        X2 = sorted([int(dist) for dist in counts.keys()])
        y2 = [counts[dist] for dist in X2]

        ax2 = plt.subplot(gs2[row, col])
        ax2.plot(X2, y2)
        ax2.set_xlim(-1000, 1000)
        ax2.set_ylim(0)
        ax2.set_title(tf)

        if not ax2.is_last_row():
            plt.setp(ax2.get_xticklabels(), visible=False)
        if not ax2.is_first_col():
            plt.setp(ax2.get_yticklabels(), visible=False)

        col += 1
        if col >= 6:
            row += 1
            col = 0
        if row >= 5:
            break

    fig.savefig(
        "{}/figures/temp/{}_{}{}.png".format(
            WORKSPACE, RUN_ID, mut_dataset_name, SUFFIX
        ),
        dpi="figure",
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )


### Actually run the datasets ###


if WHICH_DATA == "small":
    # small datasets, for testing purposes
    all_names = ["BLCA", "COAD", "HNSC", "LUAD", "READ"]
elif WHICH_DATA == "skcm":
    # melanomas, for use with distalTFBS-DHS_skcm.bed
    all_names = ["MELA", "SKCA", "SKCM"]
elif WHICH_DATA == "dhs":
    # cancers w/ DHS data, for use with merged_ENCODE.tf.bound.union.bed
    all_names = [
        "BRCA",
        "COAD",
        "COCA",
        "LUAD",
        "LUSC",
        "MELA",
        "READ",
        "SKCA",
        "SKCM",
    ]
elif WHICH_DATA == "all":
    # everything
    all_names = [
        "BLCA",
        "BRCA",
        "COAD",
        "COCA",
        "HNSC",
        "LUAD",
        "LUSC",
        "MELA",
        "READ",
        "SKCA",
        "SKCM",
    ]
else:
    # single cancer type
    all_names = [WHICH_DATA]

all_counts = {}
all_counts_by_tf = {}
for name in all_names:
    all_counts[name], all_counts_by_tf[name] = get_dists(name)
    plot_dists(all_counts[name], all_counts_by_tf[name], name)
