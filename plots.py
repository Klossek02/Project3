import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

try:
    from upsetplot import UpSet, from_contents
except ImportError:
    UpSet = None
    from_contents = None

def plot_len_by_type(all_sv, out_dir):
    for tool, svs in all_sv.items():
        lengths_by_type = defaultdict(list)
        for _, svtype, start, end in svs:
            lengths_by_type[svtype].append(abs(end - start))

        plt.figure(figsize=(8,5))
        for svtype, lens in lengths_by_type.items():
            plt.hist(lens, bins=30, alpha=0.5, label=svtype)
        plt.title(f"SV length distribution – {tool}")
        plt.xlabel("length [bp]"); plt.ylabel("#SV")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"hist_by_type_{tool}.png"))
        plt.close()

def jaccard(a, b, tol):
    def key(t):
        c, t, s, e = t ; return f"{c}:{t}:{round(s/tol)}:{round(e/tol)}"
    sa = set(map(key, a)); sb = set(map(key, b))
    return len(sa & sb) / len(sa | sb) if sa or sb else 0.0

def plot_jaccard_heat(all_sv, tol, out_dir):
    tools = list(all_sv.keys())
    m = np.zeros((len(tools), len(tools)))
    for i,t1 in enumerate(tools):
        for j,t2 in enumerate(tools):
            m[i,j] = jaccard(all_sv[t1], all_sv[t2], tol)

    df = pd.DataFrame(m, index=tools, columns=tools)
    plt.figure(figsize=(6,5))
    plt.imshow(df, interpolation="nearest")
    plt.xticks(range(len(tools)), tools, rotation=45, ha="right")
    plt.yticks(range(len(tools)), tools)
    for i in range(len(tools)):
        for j in range(len(tools)):
            val = m[i, j]
            color = "black" if val > 0.5 else "white"
            plt.text(j, i, f"{val:.2f}",
                     ha="center", va="center", fontsize=8, color=color)
    plt.colorbar(label="Jaccard index")
    plt.title("Overlap similarity (Jaccard)")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "jaccard_heatmap.png"))
    plt.close()

def plot_upset(all_sv, tol, out_dir):
    contents = {}
    for tool, svs in all_sv.items():
        contents[tool] = {f"{sv[0]}:{sv[1]}:{round(sv[2]/tol)}:{round(sv[3]/tol)}"
                          for sv in svs}
    upset_data = from_contents(contents)
    plt.figure(figsize=(8,5))
    UpSet(upset_data, subset_size='count').plot()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "upset.png"))
    plt.close()

def plot_svtype_bar(all_sv, out_dir):
    svtypes = sorted({sv[1] for lst in all_sv.values() for sv in lst})
    tool_names = list(all_sv.keys())
    counts = np.zeros((len(tool_names), len(svtypes)), dtype=int)

    for i, tool in enumerate(tool_names):
        for sv in all_sv[tool]:
            counts[i, svtypes.index(sv[1])] += 1

    x = np.arange(len(tool_names))
    width = 0.15
    plt.figure(figsize=(9,5))
    for j, svt in enumerate(svtypes):
        plt.bar(x + j*width, counts[:,j], width, label=svt)
    plt.xticks(x + width*len(svtypes)/2, tool_names, rotation=45, ha="right")
    plt.ylabel("#SV")
    plt.legend(title="SVTYPE")
    plt.title("SVTYPE composition per tool")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "type_bar.png"))
    plt.close()


