# coding=utf-8

import sys
import igraph as ig
from collections import defaultdict as dd


def parse_fasta(fasta_file):
    fasta = {}
    seq = ""
    acc = ""
    with open(fasta_file, 'r') as f:
        for l in f:
            line = l.rstrip("\n")
            if(line and l[0] == ">"):
                if(seq != ""):
                    fasta[acc] = seq
                    seq = ""
                acc = line[1:]
            else:
                seq += line
        fasta[acc] = seq
    return(fasta)


def complexity_score(kmer):
    count = dd(int)
    nuc = 2
    for i in range(len(kmer) - nuc + 1):
        count[kmer[i:i + nuc]] += 1
    return(sum(c * (c - 1) for nc, c in count.items()) / float(2 * (len(kmer) - 2)))


def test_cover(l1, l2, coverage_ref, coverage_tgt, min_cover):

    reference_cover = sum([1 for i in coverage_ref.keys()]) / float(l1)
    target_cover = sum([1 for i in coverage_tgt.keys()]) / float(l2)
    #print(reference_cover, target_cover, l1, l2)
    # print(reference_cover,target_cover)
    if(l1 <= l2 and reference_cover >= min_cover):
        return(True)

    elif(l2 <= l1 and target_cover >= min_cover):
        return(True)
    return(False)


# decide if seed ocurences is compatible with both the reads being isoforms.
def is_iso(l1, l2, pos_list, k=16, ks=3, max_diff_rate=10.0, min_cover=0.0):

    coverage_ref = dd(bool)
    coverage_tgt = dd(bool)
    last_ref = pos_list[0][0]
    last_tgt = pos_list[0][1]

    for p in pos_list:
        # computing coverage
        for j in range(k):
            coverage_ref[p[0] + j] = True
            coverage_tgt[p[1] + j] = True

        # testing seed coherence

        # how much did we advanced on reference and target ?
        dif_ref = abs(p[0] - last_ref)
        dif_tgt = abs(p[1] - last_tgt)

        # find the biggest gap
        biggest_dif = max(dif_ref, dif_tgt)
        # and smallest
        smallest_dif = min(dif_ref, dif_tgt)

        # maximum allowed difference
        max_dif = ks
        if(smallest_dif != 0):
            max_dif = max(smallest_dif * max_diff_rate, max_dif)

        # offset difference  between the two read
        relative_dif = abs(dif_ref - dif_tgt)

        if(biggest_dif > max_dif):
            # if not, discard the mapping
            #print(l1,l2, relative_dif,  smallest_dif * max_diff_rate, biggest_dif)
            # print("END")
            # exit()
            return(False)
    # else, check if cover is high enough
    return(test_cover(l1, l2, coverage_ref, coverage_tgt, min_cover))



# Edge graph
g = ig.Graph()

fasta = parse_fasta(sys.argv[2])


lines = []
# parsing edge file
with open(sys.argv[1]) as edge_file:
    for line in edge_file:
        lines.append(line)


def run(lines, fasta, mdr, mc):

    TP = 0
    TN = 0
    FP = 0
    FN = 0

    # current mapping
    for i, line in enumerate(lines):

        if(i > 1):
            line = line.rstrip("\n")
            data = line.split("\t")
            read1 = data[0]
            l1 = int(data[1])
            read2 = data[2]
            l2 = int(data[3])
            # orientation = data[4]
            pos = [tuple(int(a) for a in el.split(",")) for el in data[5:]]

            # checking if reads are the same isoforms or not.
            ref_gene = read1.split("_")[-1]
            tgt_gene = read2.split("_")[-1]
            same_iso = is_iso(l1, l2, pos, max_diff_rate=mdr, min_cover=mc)
            # print(ref_gene,tgt_gene)
            if(same_iso):
                if(ref_gene == tgt_gene):
                    TP += 1
                else:
                    FP += 1
            else:
                if(ref_gene == tgt_gene):
                    FN += 1
                else:
                    TN += 1
    return(TP, FP, TN, FN)
    # print(TP, FP, TN, FN)
    # print("Precision:", round(TP / float(TP + FP), 2))
    # print("Recall   :", round(TP / float(TP + FN), 2))
    # print("Accuracy :", round((TP + TN) / float(TP + TN + FP + FN), 2))


diff_rates = list([float(i / 100) for i in range(100, 141, 2)])
min_covers = list([float(i / 100) for i in range(30,81,1)])
print("diff_rates\\cover\t" + "\t".join(str(el) for el in min_covers), end="\t")
for i in range(3):
    print("\t" + "\t".join(str(el) for el in min_covers), end="\t")
print("")

# going through parameters
for dr in diff_rates:
    print(str(dr) + "\t", end = "")
    TPdict = []
    FPdict = []
    TNdict = []
    FNdict = []
    for mc in min_covers:
        tp, fp, tn, fn = run(lines, fasta, dr, mc)
        TPdict.append(str(tp))
        FPdict.append(str(fp))
        TNdict.append(str(tn))
        FNdict.append(str(fn))

    print(*TPdict, sep="\t", end ="\t\t")
    print(*FPdict, sep="\t", end ="\t\t")
    print(*TNdict, sep="\t", end ="\t\t")
    print(*FNdict, sep="\t", end ="\n")
