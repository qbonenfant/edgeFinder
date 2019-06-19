# coding=utf-8

import sys
import os
from collections import defaultdict as dd
import networkx as nx
from draw_blocks import draw_map, lane
from svgwrite import rgb


error_rate = 0.20
coverage = 0.80


def parse_fasta(fasta_file):
    fasta = {}
    seq = ""
    acc = ""
    with open(fasta_file, 'r') as f:
        for l in f:
            l = l.rstrip("\n")
            if(l and l[0] == ">"):
                if(seq != ""):
                    fasta[acc] = seq
                    seq = ""
                acc = l[1:]
            else:
                seq += l
        fasta[acc] = seq
    return(fasta)


def is_block(last1, last2, p1, p2):
    if(abs(last1 - p1) <= 10 and abs(last2 - p2) <= 10):
        return(True)
    return(False)


def find_blocks(pos_list):
    """
    Process mapped positions in order to produce "coherent contiguous blocks",
    like a simple mapping.
    """
    blocks = []
    current_block = []
    last = None
    for pos in pos_list:
        p1, p2 = pos
        if(last):
            dif1 = abs(p1 - last[0])
            dif2 = abs(p2 - last[1])
            if(dif1 != 0 and dif2 != 0):
                # if the difference between current and last pos for both read
                # is not too different, we increase the bloc size
                if(is_block(*last, p1, p2)):
                    current_block.append((p1, p2))
                elif(current_block):
                    blocks.append(current_block)
                    current_block = [(p1, p2)]
        else:
            current_block.append((p1, p2))
        last = (p1, p2)
    blocks.append(current_block)
    return(blocks)


def check_isoform(pos1, pos2, l1, l2):
    cover1 = float(abs(pos1[0] - pos1[-1])) / l1
    cover2 = float(abs(pos2[0] - pos2[-1])) / l2
    # if( abs(1 - float(l1) / l2 ) <=  error_rate):# and
    if(cover1 >= coverage or cover2 >= coverage):
        return(True)
    return(True)


def compute_pileup(lanes):
    """
    Add all lane local positions into a normaliszed mapping pileup
    """
    mx = 0
    pileup = dd(int)
    # First getting all values
    for l in lanes:
        for k, v in l.pileup.items():
            pileup[k] += v
            if(pileup[k] > mx):
                mx = pileup[k]
    # normalizing, max 50 pixels, relative to max value or small value of 10 px
    mx = max(mx, 10)
    for p in pileup:
        pileup[p] = (pileup[p] * 50) // mx
    return(pileup)


def complexity_score(kmer):
    count = dd(int)
    nuc = 2
    for i in range(len(kmer) - nuc + 1):
        count[kmer[i:i + nuc]] += 1
    return(sum(c * (c - 1) for nc, c in count.items()) / float(2 * (len(kmer) - 2)))


def color_ramp(score, center, mx):

    # B = 100 * ((center - abs(center - score)) / center)
    # B = 100 if score == mx else 0
    R = 75 * (score / mx) # - 1.0 * B
    V = 100 * ((mx - score) / mx) #  - 1.0 * B
    return(rgb(R, V, 0, '%'))


def complexity_legend(dwg, center, mx, height):

    # Linear tricolor scale
    x, y = dwg.origin
    x = int(x // 2)
    y = int(y - height)
    for i in range(101):
        score = mx * i / 100
        col = color_ramp(score, center, mx)
        x_pos = x + i * dwg.h_factor
        st_width = dwg.h_factor
        dwg.get_draw().add(dwg.get_draw().line(start=(x_pos, y),
                                               end=(x_pos, y + height),
                                               stroke_width=st_width,
                                               stroke=col))
    # Drawing axis
    dwg.get_draw().add(dwg.get_draw().line(start=(x, y),
                                           end=(x + 100 * dwg.h_factor, y),
                                           stroke_width=1,
                                           stroke="black"))
    # drawing ticks and text:
    for i, txt in enumerate(["0", str(center), ">" + str(mx)]):
        x_pos = dwg.h_factor * i * 50 + x
        dwg.get_draw().add(dwg.get_draw().line(start=(x_pos, y),
                                               end=(x_pos, y + 5),
                                               stroke_width=1,
                                               stroke="black"))
        dwg.draw_text((x_pos - dwg.char_width, y - 1), txt)

    # Description
    msg = "LC_scale"
    dwg.draw_text((x - dwg.char_width * (len(msg) + 1),
                   y),
                  msg)


def compute_complexity_scale(dwg, pos, height, sequence, k):
    scale = []
    center = 0.6    # center point for scale
    mx = center * 2  # max value
    y = pos[1]
    for i in range(len(sequence) - k + 1):
        score = min(complexity_score(sequence[i:i + k]), mx)
        col = color_ramp(score, center, mx)

        x_pos = dwg.origin[0] + i

        dwg.draw_dash((x_pos, y), height, line_color=col)

    complexity_legend(dwg, center, mx, height)


def print_edges(ref, len_ref, mapping, output_folder="./", display_t="blocks", ref_seq="", k=16):
    """
    Display 'alignements' in SVG format.
    """

    pileup_margin = 50
    # drawing object
    dwg = draw_map(os.path.join(output_folder + ref), font_size=34)

    nb_read = len(mapping)
    name_margin = (len(ref) + 2) * (dwg.char_width) + 20

    width = dwg.LEFT_MARGIN + name_margin + len_ref + dwg.RIGHT_MARGIN

    # Axis origin
    dwg.set_origin((dwg.LEFT_MARGIN + name_margin,
                    dwg.TOP_MARGIN +
                    max(dwg.BOX_HEIGHT, dwg.char_height) +
                    dwg.SPACER +
                    pileup_margin
                    ))

    # Generating lanes
    lane_nb = 1
    # first lane is reference
    lanes = [lane(dwg, ref, "0", {}, 0, is_ref=True)]
    pileup = dd(int)

    sort_map = sorted(mapping.items(), key=lambda x: (x[0][0], len(x[1])))
    for (orientation, read), blocks in sort_map:
        lanes.append(lane(dwg, read, orientation, blocks,
                          lane_nb, display_type=display_t))
        lane_nb += 1

    # Creating background
    dwg.set_bg()

    # Printing lanes
    for l in lanes:
        l.print_lane()

    # Printing complexity scale
    compute_complexity_scale(dwg, dwg.origin, lanes[0].height, ref_seq, k)

    # Drawing the axis
    dwg.draw_x_axis(len_ref)

    # Computing pileup
    pl = compute_pileup(lanes)
    for k, v in pl.items():
        x = dwg.origin[0] + k
        y = dwg.origin[1] - dwg.char_height - 1
        dwg.draw_dash((x, y), -v, line_color=dwg.BORDER_GREEN)

    dwg.save()


# g = nx.Graph()
bdisp = False

if(len(sys.argv) > 2):
    if(sys.argv[2] == "blocks"):
        bdisp = True

fasta = {}
if(len(sys.argv) > 3):
    fasta = parse_fasta(sys.argv[3])

with open(sys.argv[1]) as edge_file:

    # current mapping
    ref = ""
    len_ref = 0
    mapping = {}
    blocks = {}
    folder = "./out/"
    pos = []

    for i, line in enumerate(edge_file):

        if(i > 1):
            line = line.rstrip("\n")
            data = line.split("\t")
            read1 = data[0]
            l1 = int(data[1])
            read2 = data[2]
            l2 = int(data[3])
            orientation = data[4]
            pos = [tuple(int(a) for a in el.split(",")) for el in data[5:]]

            # g.add_edge(read1, read2, nk=len(pos), ori=orientation)

            if(not ref):
                ref = read1
                len_ref = l1
                ref_seq = fasta[ref] if fasta else ""

            if(read1 != ref):

                if(bdisp):
                    print_edges(ref, len_ref, blocks, folder, ref_seq=ref_seq)
                else:
                    print_edges(ref, len_ref, mapping,
                                folder, display_t="dash", ref_seq=ref_seq)

                mapping = {}
                ref = read1
                len_ref = l1
                ref_seq = fasta[ref] if fasta else ""
                # break

            mapping[(orientation, read2)] = pos
            blocks[(orientation, read2)] = find_blocks(pos)
    if(bdisp):
        print_edges(ref, len_ref, blocks, folder, ref_seq=ref_seq)
    else:
        print_edges(ref, len_ref, mapping, folder,
                    display_t="dash", ref_seq=ref_seq)
