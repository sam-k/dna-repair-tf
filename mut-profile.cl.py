import collections
from matplotlib import pyplot as plt
from matplotlib import gridspec as gs


WORKSPACE = "/data/gordanlab/samkim/dna-repair-tf"
FILTER_BY = "WGS"
TFBS_DHS = True
TFBS_TYPE = "proximal"
WHICH_DATA = "skcm"


def get_dists(mut_dataset_name, filter_by="", tfbs_dhs=True, tfbs_type="proximal"):
    mut_list = []
    with open(
        "{}/data/ssm.open.{}{}{}_{}_centered.bed".format(
            WORKSPACE,
            ("" if tfbs_type == "proximal" else tfbs_type + "TFBS_"),
            (filter_by + "_" if filter_by else ""),
            ("DHS" if tfbs_dhs else "noDHS"),
            mut_dataset_name,
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


def plot_dists(
    counts,
    counts_by_tf,
    mut_dataset_name=None,
    figsize=(10, 14),
    h1=0.2,
    h2=0.6,
    w2=0.5,
):
    plt.figure(figsize=figsize)
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

    plt.show()


def make_plots(all_counts, all_counts_by_tf, all_names, name):
    if isinstance(name, int) or name.isnumeric():
        if int(name) < len(all_names):
            name = all_names[int(name)]
        else:
            print("There are no more datasets.")
            return
    plot_dists(all_counts[name], all_counts_by_tf[name], name)


if WHICH_DATA == "small":
    all_names = ["BLCA", "COAD", "HNSC", "LUAD", "READ"]
elif WHICH_DATA == "skcm":
    all_names = ["MELA", "SKCA", "SKCM"]
else:
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

all_counts = {}
all_counts_by_tf = {}
for name in all_names:
    all_counts[name], all_counts_by_tf[name] = get_dists(
        name, FILTER_BY, TFBS_DHS, TFBS_TYPE
    )

name_index = 0
