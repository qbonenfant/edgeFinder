# coding=utf-8

import sys
from collections import defaultdict as dd
from math import sqrt
from functools import reduce

import matplotlib.pyplot as plt


def get_med(elem_list):
    elem_list.sort()
    ln = len(elem_list)
    b = round(ln / 2)
    s = ln // 2
    return(float((elem_list[s] + elem_list[b]) / 2))


def get_dist(elem_list):
    res = 0
    last = elem_list[0]
    for el in elem_list:
        res += abs(last - el)
        last = el
    return(res)


def get_stdev(elem_list):
    res = 0
    mean = float(sum(elem_list) / len(elem_list))
    for el in elem_list:
        res += abs(el - mean)**2
    return(res / len(elem_list))


def cover(l, coverage):
    return(len(coverage) / float(l))


def check_sum(it):

    it = str(it)
    res = 0
    for c in it:
        res ^= (ord(c)**2 * (res + 1) // 3) % 256
    return(res)


def get_coverage(pos_list, length, k=15):
    """get coverage for a position list with known read and kmer length"""
    coverage = [0] * length
    for p in pos_list:
        # computing coverage
        for j in range(k):
            if(p + j < length):
                coverage[p + j] = 1
    return(coverage)


def get_coherence(pos_pair_list):
    """Testing seed coherence"""

    # keeping relative differences, this is out observation set.
    diff_vec = []
    last_ref = pos_pair_list[0][0]
    last_tgt = pos_pair_list[0][1]
    for p in pos_pair_list:
        # how much did we advanced on reference and target ?
        dif_ref = abs(p[0] - last_ref)
        dif_tgt = abs(p[1] - last_tgt)

        # offset difference  between the two read
        relative_diff = abs(dif_ref - dif_tgt)
        diff_vec.append(relative_diff)

        last_ref = p[0]
        last_tgt = p[1]
    return(diff_vec)


def get_medians_dist(seed_list, length):
    dist = {}
    for i in range(length):
        dist[i] = 0
    for r in seed_list.keys():
        md = round(get_med(list(seed_list[r])))
        dist[md] += 1
    return(dist)


def AB_exp_diff(seed_dict, r_len):
    reads = list(seed_dict.keys())

    # getting all unique seeds occurences
    seeds = set()
    for r in reads:
        seeds |= set(seed_dict[r])

    score = {}
    # Using seeds as pivot.
    for s in sorted(list(seeds))[1:]:
        exp = abs(s - (r_len - s)) / float(r_len)
        c_score = 0
        for r in reads:
            fl = float(len(seed_dict[r]))
            A = len(list(el for el in seed_dict[r] if el >= s)) / fl
            B = len(list(el for el in seed_dict[r] if el < s)) / fl
            c_score += abs((A - B)**2 - exp**2)
        score[s] = c_score / float(len(reads))
        # score[s] = 1 if len((L & R)) <= 1 else 0
    return(score)


def simple_AB_exp_diff(seed_pairs, l1, l2, limit=0.8):

    # Assuming both are not chimeras
    source_chimera = False
    target_chimera = False
    s_seeds, t_seeds = zip(*seed_pairs)

    # checking source first
    for s in s_seeds[1:]:
        exp = abs(s - (l1 - s)) / float(l1)
        fl = float(len(s_seeds))
        A = len(list(el for el in s_seeds if el >= s)) / fl
        B = len(list(el for el in s_seeds if el < s)) / fl
        c_score = abs((A - B)**2 - exp**2)
        if(c_score >= limit):
            source_chimera = True
            break

    # checking target
    for s in t_seeds[1:]:
        exp = abs(s - (l2 - s)) / float(l2)
        fl = float(len(t_seeds))
        A = len(list(el for el in t_seeds if el >= s)) / fl
        B = len(list(el for el in t_seeds if el < s)) / fl
        c_score = abs((A - B)**2 - exp**2)
        if(c_score >= limit):
            target_chimera = True
            break

    return(source_chimera, target_chimera)


def lr_local_prop(seed_dict, read, r_len):
    reads = list(seed_dict.keys())

    # getting all unique seeds occurences
    seeds = set()
    for r in reads:
        seeds |= set(seed_dict[r])

    score = {}

    # Using seeds as pivot.
    for s in sorted(list(seeds))[1:-1]:

        exp = abs(s - (r_len - s)) / float(r_len)

        # Seed # on the left
        A_loc = len(list(el for el in seed_dict[read] if el <= s))
        A_tot = len(list(el for el in seeds if el <= s))
        # Seed # on the right
        B_loc = len(list(el for el in seed_dict[read] if el > s))
        B_tot = len(list(el for el in seeds if el > s))

        score[s] = abs((float(A_loc) / A_tot - float(B_loc) / B_tot) )
        # score[s] = A_loc / float(A_tot)

    return(score)


def simple_R(seed_dict):
    """ Test for the proportion of seeds on the right of a positon"""

    reads = list(seed_dict.keys())

    # getting all unique seeds occurences
    seeds = set()
    for r in reads:
        seeds |= set(seed_dict[r])

    score = {}
    # Using seeds as pivot.
    for s in sorted(list(seeds))[1:]:
        c_score = 0.0
        for r in reads:
            fl = float(len(seed_dict[r]))
            R = len(list(el for el in seed_dict[r] if el > s)) / fl
            c_score += R
        score[s] = c_score / float(len(reads))
        # score[s] = 1 if len((L & R)) <= 1 else 0
    return(score)


def simple_diff(seed_dict):

    reads = list(seed_dict.keys())

    # getting all unique seeds occurences
    seeds = set()
    for r in reads:
        seeds |= set(seed_dict[r])

    score = {}
    # Using seeds as pivot.
    for s in sorted(list(seeds))[1:]:
        c_score = 0.0
        score_dict = {}
        for i, r in enumerate(reads):
            fl = float(len(seed_dict[r]))
            R = len(list(el for el in seed_dict[r] if el > s))
            L = len(list(el for el in seed_dict[r] if el <= s))

            R = 1 if R / fl > 0.05 else 0
            L = 1 if L / fl > 0.05 else 0

            score_dict[i] = (R - L)**2

        c_score = 0
        for rid in score_dict:
            if(score_dict[rid]):
                c_score ^= rid

        score[s] = c_score
        # score[s] = 1 if len((L & R)) <= 1 else 0
    return(score)


def plot_chimera(read_list, edge_data):

    fig, ax = plt.subplots(2, 5, figsize=(40, 5))
    plt.subplots_adjust(wspace=.3, hspace=.5, left=.05, right=.95)
    nb_plot = 0
    nb_chimera = 0
    nb_normal = 0
    for read in read_list:
        if(("chimera" in read and nb_chimera < 5) or
          ("chimera" not in read and nb_normal < 5)) and\
            len(edge_data[read]) == 2:

            # short_read_name = read.split("_")[0]
            short_read_name="chimeric" if "chimera" in read else "normal"
            col="r" if "chimera" in read else "b"

            # pl_x, pl_y = nb_plot % 2, nb_plot // 2
            pl_x=1 if col == "r" else 0
            pl_y=nb_chimera if pl_x == 1 else nb_normal

            # Computing data dict
            #######################################################
            # # AB-EXP delta squared
            # res = AB_exp_diff(edge_data[read], read_length[read])
            # x, y = zip(*sorted([(k, v) for k, v in res.items()]))
            # ax[pl_x, pl_y].scatter(
            #     x, y, color=col, marker="+", label=short_read_name)
            # ax[pl_x, pl_y].set(xlabel="position",
            #                    ylabel=r"$\sum(\delta(A,B)²,exp²)/N$")

            #######################################################

            #######################################################
            # Local proportion
            # # for r in edge_data[read].keys():
            r=list(edge_data[read].keys())[0]
            res=lr_local_prop(edge_data[read], r, read_length[read])
            x, y=zip(*sorted([(k, v) for k, v in res.items()]))
            # y1, y2=zip(* y)
            ax[pl_x, pl_y].scatter(
                x, y, marker = "+", color = col, label = short_read_name)
            ax[pl_x, pl_y].set(xlabel = "position",
                               ylabel = r"$A(r2)/(A(r2)+A(r3))$")


            # ax[pl_x, pl_y].scatter(
            # x, y2, color = 'b', marker = "+", label = short_read_name + "_B")
            #######################################################

            #######################################################
            # Simple R seed count
            # res = simple_diff(edge_data[read])
            # x, y = zip(*sorted([(k, v) for k, v in res.items()]))
            # # y1, y2=zip(* y)
            # ax[pl_x, pl_y].scatter(
            #     x, y, marker = "+", color= col, label=short_read_name)

            #######################################################

            ax[pl_x, pl_y].legend(loc = 'lower left', bbox_to_anchor = (
                0.0, 1.01), borderaxespad = 0, frameon = False)

            # ax[pl_x, pl_y].set(xlabel = "position",
            #                    ylabel = "L/R prop")

                               # ylim=(0, 1))

            if(col == "r"):
                nb_chimera += 1
            else:
                nb_normal += 1
            nb_plot += 1

        if(nb_plot >= 10):
            break

    plt.show()


def parse_edge(edge_file_name):
    edge_data=dd(dict)
    with open(edge_file_name) as edge_file:
        for line in edge_file:
            line=line.rstrip("\n")
            data=line.split("\t")
            if(len(data) > 2):
                read1=data[0]
                l1=int(data[1])
                read2=data[2]
                l2=int(data[3])
                # orientation = data[4]
                pos=[tuple(int(a) for a in el.split(",")) for el in data[6:]]

                read_set.add(read1)

                # checking if reads are the same isoforms or not.
                # ref_gene = read1.split("_")[-1]
                # tgt_gene = read2.split("_")[-1]

                read_length[read1]=l1
                read_length[read2]=l2

                # exporting respectiv cover of each read.
                r_cover, t_cover=zip(*pos)

                # storing edge data
                # edge_data[read1][read2] = get_coverage(r_cover, l1)
                # edge_data[read2][read1] = get_coverage(t_cover, l2)
                edge_data[read1][read2]=r_cover
                edge_data[read2][read1]=t_cover
    return(edge_data)


# compute stdev score for the two reads of each line.
def single_line_score(edge_file_name):
    chimera_count = dd(int)
    with open(edge_file_name) as edge_file:
        for line in edge_file:
            line = line.rstrip("\n")
            data = line.split("\t")
            if(len(data) > 2):
                read1 = data[0]
                l1 = int(data[1])
                read2 = data[2]
                l2 =int(data[3])
                # orientation = data[4]
                pos = [tuple(int(a) for a in el.split(",")) for el in data[6:]]
                sc, tc = simple_AB_exp_diff(pos, l1, l2)

                chimera_count[read1] += 1 if sc else 0
                chimera_count[read2] += 1 if tc else 0

    # printing results
    for k in sorted(chimera_count.keys(), key= lambda x: chimera_count[x], reverse = True):
        print(k, chimera_count[k] , sep="\t")




##############################################################################
# Main start
edge_file_name=sys.argv[1]



read_length = {}
read_set = set()

# Parsing edge file and keeping seed occurences for each read.
edge_data = parse_edge(edge_file_name)

read_list = []
all_reads = set(edge_data.keys())
# read_list = list(all_reads - read_set)
# read_list = list(read_set)
# read_list = list(all_reads)

# checking if a list of read has been given
if(len(sys.argv) > 2):
    with open(sys.argv[2]) as rl:
        for line in rl:
            read_list.append(line.rstrip("\n"))

# else, using all reads as input
if(len(read_list) == 0):
    read_list=list(all_reads)

##############################################################################
# Data processing

# Testing metric
# TP, TN, FP, FN=0, 0, 0, 0

for read in read_list:

    res=None
    mapped_reads=list(edge_data[read].keys())

    met = 0
    if(len(mapped_reads) > 1):
        res = AB_exp_diff(edge_data[read], read_length[read])
        key, val = zip(*sorted([(k, v) for k, v in res.items()]))
        met = max(val)

    if(met > 0.90):
        if("chimera" in read):
            TP += 1
        else:
            FP += 1
    else:
        if("chimera" not in read):
            TN += 1
        else:
            FN += 1

print("TP", "TN", "FP", "FN")
print(TP, TN, FP, FN)
print("Rand index:")
val=TP + TN
tot=TP + FN + TN + FP
print(round(float(val / tot), 2))

# plot_chimera(read_list, edge_data)
