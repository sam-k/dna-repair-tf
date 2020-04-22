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
SUFFIX = "" if len(sys.argv) <= 3 else sys.argv[3]
# suffix to filename (currently only for merged and merged-bg)


### Get mutation counts per position, and per position per TF ###


def get_dists(mut_dataset_name=None, suffix=None):
    if mut_dataset_name is None:
        return {}, {}

    mut_list = []
    with open(
        "{}/data/ssm.open.{}_{}{}_centered.bed".format(
            WORKSPACE, RUN_ID, mut_dataset_name, "" if suffix is None else "_" + suffix,
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
    counts_2=None,
    counts_2_by_tf=None,
    name=None,
    figsize=(10, 14),
    h1=0.2,
    h2=0.6,
    w2=0.5,
    center_lims=True,
    save_fig=True,
):
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(hspace=h1)
    outer = gs.GridSpec(2, 1, height_ratios=[2, 5])
    gs1 = gs.GridSpecFromSubplotSpec(1, 6, subplot_spec=outer[0])
    gs2 = gs.GridSpecFromSubplotSpec(5, 6, subplot_spec=outer[1], hspace=h2, wspace=w2)

    # Create large graph of all TFs
    ax1 = plt.subplot(gs1[:, 1:-1])

    # Default, or bound DHS
    X1 = sorted([int(dist) for dist in counts])
    y1 = [counts[dist] for dist in X1]
    ax1.plot(X1, y1)

    # Unbound DHS
    X1_2 = []
    y1_2 = []
    if counts_2 is not None and counts_2_by_tf is not None:
        X1_2 = sorted([int(dist) for dist in counts_2])
        y1_2 = [counts_2[dist] for dist in X1_2]
        ax1.plot(X1_2, y1_2)

    # Style large graph
    if center_lims:
        ax1.set_xlim(-1000, 1000)
    # else:
    #     min_lim = min(min(X1), min(X1_2)) if len(X1_2) > 0 else min(X1)
    #     max_lim = max(max(X1), max(X1_2)) if len(X1_2) > 0 else max(X1)
    #     ax1.set_xlim(min_lim, max_lim)
    ax1.set_ylim(0)
    ax1.set_xlabel("Distance from TFBS center (bp)")
    ax1.set_ylabel("Number of mutations")
    if name is None:
        ax1.set_title("Mutation profile for all TFs")
    else:
        ax1.set_title("{} mutation profile for all TFs".format(name))

    # Create small graphs of every TF
    ordered_counts = sorted(
        counts_by_tf.items(), key=lambda t: sum(t[1].values()), reverse=True
    )
    row, col = 0, 0
    for tf, tf_counts in ordered_counts:
        ax2 = plt.subplot(gs2[row, col])

        # Default, or bound DHS
        X2 = sorted([int(dist) for dist in tf_counts])
        y2 = [tf_counts[dist] for dist in X2]
        ax2.plot(X2, y2)

        # Unbound DHS
        X2_2 = []
        y2_2 = []
        if counts_2 is not None and counts_2_by_tf is not None and tf in counts_2_by_tf:
            X2_2 = sorted([int(dist) for dist in counts_2_by_tf[tf]])
            y2_2 = [counts_2_by_tf[tf][dist] for dist in X2_2]
            ax2.plot(X2_2, y2_2)

        # Style small graph
        if center_lims:
            ax2.set_xlim(-1000, 1000)
        # else:
        #     min_lim = min(min(X2), min(X2_2)) if len(X2_2) > 0 else min(X2)
        #     max_lim = max(max(X2), max(X2_2)) if len(X2_2) > 0 else max(X2)
        #     ax2.set_xlim(min_lim, max_lim)
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

    # Save figure
    if save_fig:
        fig.savefig(
            "{}/figures/temp/{}_{}{}.png".format(
                WORKSPACE, RUN_ID, name, ("-" if len(SUFFIX) > 0 else "") + SUFFIX
            )
            .replace("proximal", "prox")
            .replace("distal", "dist")
            .replace("DHS_", "DHS-"),
            dpi="figure",
            transparent=True,
            bbox_inches="tight",
            pad_inches=0,
        )
    else:
        plt.show()


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
all_counts_2 = {}
all_counts_2_by_tf = {}
for name in all_names:
    if RUN_ID.endswith("mergedbg"):
        all_counts[name], all_counts_by_tf[name] = get_dists(
            name, SUFFIX + ("_" if len(SUFFIX) > 0 else "") + "bound"
        )
        all_counts_2[name], all_counts_2_by_tf[name] = get_dists(
            name, SUFFIX + ("_" if len(SUFFIX) > 0 else "") + "unbound"
        )
        plot_dists(
            all_counts[name],
            all_counts_by_tf[name],
            all_counts_2[name],
            all_counts_2_by_tf[name],
            name=name,
            center_lims=False,
        )
    else:
        all_counts[name], all_counts_by_tf[name] = get_dists(name)
        plot_dists(all_counts[name], all_counts_by_tf[name], name=name)
