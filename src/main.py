

import os
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image

from typing import Generator, Any


def get_fontpath() -> Path:
    fpath = Path(mpl.get_data_path(), "resources/Cambria.ttf")
    return fpath


def read_inp_names(filename: str) -> list[list[str]]:
    name_sets: list[list[str]] = []
    with open(filename, "r") as fh:
        lines: list[str] = fh.readlines()
    for line in lines:
        name_sets.append(line.split()[::-1])
    return name_sets


def make_avg(txt: str) -> str:
    avg: str = txt
    try:
        float(txt)
        avg = txt
    except ValueError:
        avg = "0"
        a: list[str] = []
        b: list[str] = []
        for char in txt:
            if char.isnumeric():
                a.append(char)
            else:
                break
        txt = txt[::-1]
        for char in txt:
            if char.isnumeric():
                b.append(char)
            else:
                break
        b.reverse()
        ax = float("".join(a))
        bx = float("".join(b))
        avg = str((ax + bx) / 2)
    return avg


def mock_input() -> str:
    data_name = r"resources\Data7.csv"
    return data_name


def mock_input_2() -> dict[str, list[float]]:
    key_names: list[str] = ['genom', 'gc', 'n2', 'velChrom', 'vaha', 'delka']
    data_per_category: dict[str, list[float]] = {}
    for key in key_names:
        data_per_category[key] = []
    return data_per_category
    

def populate_data_dct(lines: list[str], col_count: int) -> dict[str, list[list[float]]]:
    data_dct: dict[str, list[list[float]]] = {}
    for _, line in enumerate(lines):
        lst: list[str] = line.split(";")
        tax: str = lst[1]
        tmp: list[float] = []
        for j in range(2, 2+col_count):
            val: str = lst[j]
            if len(val) == 0 or val == "\n":
                tmp.append(np.nan)
            else:
                if not val.isnumeric() and val[-1:] != "\n":
                    val = make_avg(val)
                tmp.append(float(val))
        if tax not in data_dct.keys():
            data_dct[tax] = [[] for _ in range(col_count)]
        for j in range(col_count):
            data_dct[tax][j].append(tmp[j])
    return data_dct


def check_dict_validity(data_dct: dict[str, list[list[float]]]) -> None:
    keys: list[str] = list(data_dct.keys())
    for key in keys:
        if len(key) < 1:
            del data_dct[key]


def add_subfamily(data_dct: dict[str, list[list[float]]], col_count: int, name: str) -> None:
    data_dct[name] = [[] for _ in range(col_count)]
    for key in data_dct.keys():
        if key.split("-")[0] == name and len(key.split("-")) > 1:
            for i in range(col_count):
                data_dct[name][i] += data_dct[key][i]


def add_tribe(data_dct: dict[str, list[list[float]]], col_count: int, name: str) -> None:
    data_dct[name] = [[] for _ in range(col_count)]
    for key in data_dct.keys():
        if len(key.split("-")) == 1:
            continue
        if key.split("-")[0] == name.split("-")[0] and key.split("-")[1] == name.split("-")[1] and len(key.split("-")) > 2:
            for i in range(col_count):
                data_dct[name][i] += data_dct[key][i]


def check_inp_valid(name_sets: list[list[str]], graph_labels: list[list[str]]) -> bool:
    return True if [len(x) for x in name_sets] == [len(x) for x in graph_labels] else False


def get_data_per_category(data_per_category: dict[str, list[float]], name_sets: list[list[str]], data_dct: dict[str, list[float]]) -> Generator[dict[str, list[float]], Any, None]:
    val_sets: list[list[list[float]]] = [[] for _ in name_sets]
    lbl_sets: list[list[str]] = [[] for _ in name_sets]
    for i, name_set in enumerate(name_sets):
        for key in data_per_category.keys():
            data_per_category[key] = []
        for name in name_set:
            if name in data_dct.keys():
                val_sets[i].append(data_dct[name])
                lbl_sets[i].append(name)
                for j, kw in enumerate(data_per_category):
                    data_per_category[kw].append(data_dct[name][j])
            else:
                for key in data_dct.keys():
                    if name in key:
                        val_sets[i].append(data_dct[key])
                        lbl_sets[i].append(key)
                        for j, kw in enumerate(data_per_category):
                            data_per_category[kw].append(data_dct[key][j])
                        break
                else:
                    raise KeyError(f"{name}")
        yield data_per_category


def custom_plot(data_per_category, graph_num, legend_lst) -> None:

    nm: int = len(data_per_category["genom"])
    pos: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = np.linspace(1, nm, nm)

    for i, key in enumerate(data_per_category.keys()):
        for j in range(nm):
            off = 0
            for k in range(len(data_per_category[key][j])):
                if np.isnan(data_per_category[key][j][k - off]):
                    data_per_category[key][j].pop(k - off)
                    off += 1

    fig_sz: tuple[tuple[int, int], tuple[int, int], tuple[int, float]] = ((13, 6), (13, 3), (13, 10.5))

    categories: int = len(legend_lst)
    fig, ax = plt.subplots(1, categories, figsize=fig_sz[graph_num]) # 10/6 (8), 10/3 (4) , 10/10.5 (9)
    
    fig.tight_layout(pad=1.5)
    # fontsize: int = 12

    fig.subplots_adjust(wspace=0.25)

    xlims: tuple[tuple[int | float, int | float], ...] = ((1, 100), (34, 44), (1, 100), (0, 500), (0.1, 100), (0, 18))
    ticknum: tuple[int, ...] = (2, 11, 2, 11, 3, 10)
    heights: list[float] = [0.5, 0.46, 0.5]

    graph_ylabels: list[list[int]] = [[],[],[],[],[],[]]
    for i, kw in enumerate(data_per_category.keys()):
        for j in data_per_category[kw]:
            graph_ylabels[i].append(len(j))

    color_palette: list[str] = ["lightcoral", "skyblue", "wheat", "mediumspringgreen", "plum", "sandybrown"]
    color_palette_edge: list[str] = ["indianred", "steelblue", "tan", "mediumseagreen", "mediumorchid", "peru"]

    font_path = ".\\resources\\Cambria.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()

    ax[0].set_xscale("log")
    ax[2].set_xscale("log")
    ax[4].set_xscale("log")

    for i, kw in enumerate(data_per_category):

        flierprops: dict[str, str | int | float] = dict(marker='.', markerfacecolor=color_palette[i], markersize=4,
                    linestyle='none', markeredgecolor=color_palette_edge[i], markeredgewidth=0.2)
        boxprops: dict[str, str | float] = dict(linestyle='-', linewidth=0.6, color=color_palette_edge[i])
        capprops: dict[str, str | float] = dict(linestyle='-', linewidth=0.6, color=color_palette_edge[i])
        whiskerprops: dict[str, str | float] = dict(linestyle='-', linewidth=0.6, color=color_palette_edge[i])
        medianprops: dict[str, str | int] = dict(linestyle='-', linewidth=1, color=color_palette_edge[i])
        ax[i].boxplot(data_per_category[kw], positions=pos, vert=False, flierprops=flierprops, boxprops=boxprops, capprops=capprops, whiskerprops=whiskerprops, medianprops=medianprops, tick_labels=graph_ylabels[i], capwidths=0.1)
        if i in (0, 2, 4):
            if i == 4: 
                ticklst: list[int | float] = [0.1]
            else:
                ticklst: list[int | float] = [1]
            for ii in range(ticknum[i]):
                if i == 4:
                    y: list[int] = [(x) for x in np.linspace(10 ** (ii - 1), 10 ** (ii), 10)]
                    ticklst += y[1:]
                else:
                    y: list[int] = [round(x) for x in np.linspace(10 ** ii, 10 ** (ii + 1), 10)]
                    ticklst += y[1:]
            ax[i].set_xticks(ticklst)
        else:
            ax[i].set_xticks(np.linspace(xlims[i][0], xlims[i][1], ticknum[i]))
        ax[i].tick_params(axis='x', color="grey", pad=1.5, length=0.5, width=0.1, labelsize=9, grid_color="grey", grid_linewidth=0.1, labelcolor="grey", labelrotation=90)
        ax[i].tick_params(axis='y', length=0, pad=1.5, labelsize=12, grid_color="grey", grid_linewidth=0.1, labelcolor=color_palette_edge[i], labelrotation=0, labelfontfamily="sans-serif")
        ax[i].set_title(legend_lst[i], fontsize=12, color=color_palette_edge[i])
        ax[i].grid(axis="x")
        ax[i].minorticks_off()

        ax[i].spines["top"].set(linewidth=0.1, color="grey")
        ax[i].spines["bottom"].set(linewidth=0.1, color="grey")
        ax[i].spines["right"].set(linewidth=0.1, color="grey")
        ax[i].spines["left"].set(linewidth=0.1, color="grey")

        ax[i].set_xlim(xlims[i])
        ax[i].set_ylim(0.5, nm + 0.5)
        
        meds: list[int] = [0]
        for lst in data_per_category[kw]:
            meds.append(np.median(np.asarray(lst))) # type: ignore

        ax[i].barh(list(range(len(data_per_category[kw]) + 1)), meds, height=heights[graph_num], color=color_palette[i])

    savepath: str = os.path.join(os.getcwd(), f"V7b_{graph_num}_log.png")
    fig.savefig(savepath, dpi=600, format="png")

    img: Image.Image = Image.open(savepath)
    img.show()


def main() -> None:

    data_name: str = mock_input()
    data_per_category: dict[str, list[float]] = mock_input_2()
    name_sets: list[list[str]] = read_inp_names(r"src\inp.inp")
    col_count: int = 6
    legend_lst: list[str] = ["Velikost genomu [pg]", "Obsah GC [%]", "Chromozomový počet 2n", "Velikost chromozomu [Mbp]", "Váha nažky [mg]", "Délka nažky [mm]"]

    with open(data_name, "r") as fh:
        lines: list[str] = fh.readlines()

    data_dct: dict[str, list[list[float]]] = populate_data_dct(lines, col_count)
    
    check_dict_validity(data_dct) # family - subfamily - tribe - subtribe
    add_subfamily(data_dct, col_count, name="Carduoideae")
    add_tribe(data_dct, col_count, name="Carduoideae-Cardueae")

    data_per_category_gen: Generator[tuple[list[list[list[float]]], list[str]], Any, None] = get_data_per_category(data_per_category, name_sets, data_dct)
    for k in range(len(name_sets)):
        custom_plot(next(data_per_category_gen), k, legend_lst)


if __name__ == "__main__":
    main()


