# coding=utf-8

import sys
import os
from collections import defaultdict as dd
import networkx as nx
from draw_blocks import draw_map


error_rate = 0.20
coverage = 0.80


# def is_block(last1, last2, p1, p2):
    
#     global error_rate, coverage
   
#     dif1=abs(p1 - last1)
#     dif2=abs(p2 - last2)

#     if( (max(dif1, dif2) <= 10 and abs(dif1 - dif2) <= 3) or (abs(1 - float(dif1) / dif2)) <= error_rate):
#         return(True)
#     return(False)

def is_block(last1, last2, p1, p2):
    return(True)
    



def find_blocks(pos_list):
    """
    Process mapped positions in order to produce "coherent contiguous blocks", like a simple mapping.
    """
    blocks=[]
    current_block=[]
    last=None
    for pos in pos_list:
        p1, p2= pos
        if(last):
            dif1=abs(p1 - last[0])
            dif2=abs(p2 - last[1])
            if(dif1 != 0 and dif2 != 0):
                # if the difference between current and last pos for both read
                # is not too different, we increase the bloc size
                if(is_block(*last, p1, p2)):
                    current_block.append((p1, p2))
                elif(current_block):
                    blocks.append(current_block)
                    current_block=[(p1, p2)]
        else:
            current_block.append( (p1, p2) )
        last=(p1, p2)
    blocks.append(current_block)

    return(blocks)


def check_isoform(pos1, pos2, l1, l2):
    cover1=float(abs(pos1[0] - pos1[-1])) / l1
    cover2=float(abs(pos2[0] - pos2[-1])) / l2
    # if( abs(1 - float(l1) / l2 ) <=  error_rate):# and
    if(cover1 >= coverage or cover2 >= coverage):
        return(True)
    return(True)

def init_dwg():
    pass


def dash_rep(output, len_ref, mapping):

    ref=os.path.basename(output)
    dwg=draw_map(output)
    nb_read=len(mapping)
    name_margin=len(ref) * dwg.font_size

    # Image Size
    dwg.set_size(dwg.LEFT_MARGIN + name_margin + len_ref + dwg.RIGHT_MARGIN,
                   dwg.TOP_MARGIN + \
                       (nb_read + 1) * (dwg.BOX_HEIGHT +
                        dwg.SPACER) + dwg.BOTTOM_MARGIN
                )
    dwg.set_bg()
    # Axis origin
    dwg.set_origin((dwg.LEFT_MARGIN + name_margin,
                    dwg.TOP_MARGIN  + (nb_read + 1) * (dwg.BOX_HEIGHT + dwg.SPACER)
                    ))
    dwg.draw_x_axis(len_ref)
    dwg.draw_y_axis( (nb_read + 1) * (dwg.BOX_HEIGHT + dwg.SPACER) )

    # Starting to draw boxes.
    height = dwg.font_size
    offset = dwg.LEFT_MARGIN + name_margin

    # Ref box:
    dwg.draw_box( (offset, height), len_ref, dwg.BOX_HEIGHT, line_color= dwg.BORDER_GREEN, color = dwg.BOX_GREEN, direction = "")
    dwg.draw_text( (dwg.LEFT_MARGIN, height + dwg.font_size ), ref, text_col = dwg.BORDER_GREEN )

    height += dwg.SPACER + dwg.BOX_HEIGHT
    pileup = dd(int)
    for (read, orientation), positions in mapping.items():

        if(orientation == "1"):
            direction = "left"
            col = dwg.BORDER_RED
        else:
            direction = "right"
            col = dwg.BORDER_BLUE
        
        dwg.draw_text( (dwg.LEFT_MARGIN, height + dwg.font_size ), read, text_col = col )

        for p1,p2 in positions:
            pileup[p1]+=1
            dwg.draw_dash( (  offset + p1, height )  , dwg.BOX_HEIGHT  ,line_color = col )
        height += dwg.SPACER + dwg.BOX_HEIGHT

    height = dwg.origin[1] + dwg.font_size + 5

    for p,value in pileup.items():
        dwg.draw_dash( (  offset + p, height )  , value  ,line_color = dwg.BORDER_GREEN )
    dwg.save()



def block_rep(output, len_ref, mapping):
    
    ref=os.path.basename(output)
    dwg=draw_map(output)
    nb_read=len(mapping)
    name_margin=len(ref) * dwg.font_size

    dwg.set_bg()
    # Axis origin
    dwg.set_origin((dwg.LEFT_MARGIN + name_margin,
                    dwg.TOP_MARGIN  + (nb_read + 1) * (dwg.BOX_HEIGHT + dwg.SPACER)
                    ))
    dwg.draw_x_axis(len_ref)
    dwg.draw_y_axis( (nb_read + 1) * (dwg.BOX_HEIGHT + dwg.SPACER) )

    # Starting to draw boxes.
    height = dwg.font_size
    offset = dwg.LEFT_MARGIN + name_margin

    # Ref box:
    dwg.draw_box( (offset, height), len_ref, dwg.BOX_HEIGHT, line_color= dwg.BORDER_GREEN, color = dwg.BOX_GREEN, direction = "")
    dwg.draw_text( (dwg.LEFT_MARGIN, height + dwg.font_size ), ref, text_col = dwg.BORDER_GREEN )

    height += dwg.SPACER + dwg.BOX_HEIGHT
    pileup = dd(int)
    for (read, orientation), blocks in mapping.items():

        if(orientation == "1"):
            direction = "left"
            col = dwg.BOX_RED
            line_col = dwg.BORDER_RED
        else:
            direction = "right"
            col = dwg.BOX_BLUE
            line_col = dwg.BORDER_BLUE
        
        dwg.draw_text( (dwg.LEFT_MARGIN, height + dwg.font_size ), read, text_col = col )

        for block in blocks:

            start_b = block[0]
            end_b   = block[-1]


            dwg.draw_box( (offset + start_b[0], height),
                abs(end_b[0]-start_b[0]) ,
                dwg.BOX_HEIGHT,
                line_color= line_col,
                color = col,
                direction = direction,
                block_start = start_b[1],
                block_end = end_b[1]
                )

            for pu in range(start_b[0],end_b[0]+1):
                pileup[pu]+=1


        height += dwg.SPACER + dwg.BOX_HEIGHT

    height = dwg.origin[1] + dwg.font_size + 5

    for p,value in pileup.items():
        dwg.draw_dash( (  offset + p, height )  , value  ,line_color = dwg.BORDER_GREEN )
    dwg.save()


g = nx.Graph()

with open(sys.argv[1]) as edge_file:
    # edge_file = ["sequence1\t40\tsequence2\t42\t0\t1,1\t2,2\t4,5\t6,8\t30,39\t33,41"]
    # edge_file = ["ch115_read2745_template_pass_BYK_CB_ONT_1_FAF04998_A\t804 ch273_read4762_template_pass_BYK_CB_ONT_1_FAF04998_A\t3252\t0   9,1007  12,1010 18,1015 21,1018 24,1021 75,1077 78,1080 81,1084 96,1098 99,1100 102,1103    105,1106    108,1109    123,1122    126,1124    132,1130    135,1133    168,1162    174,1169    177,1172    180,1175    183,1178    186,1180    195,1193    198,1197    201,1199    204,1202    207,1205    213,1211    216,1214    219,1217    222,1220    225,1223    228,1227    231,1230    234,1233    237,1237    240,1240    243,1243    246,1246    252,1252    267,1268    270,1270    288,1288    291,1291    294,1294    297,1297    312,1309    315,1313    318,1316    321,1319    324,1322    327,1325    330,1328    333,1331    336,1334    339,1337    342,1340    345,1343    348,1346    411,1406    414,1408    483,1486    486,1488    489,1491    492,1494    495,1497    498,1500    501,1503    504,1506    555,1556    558,1559    561,1562    564,1565    567,1568    606,1604    609,1607    612,1610    615,1613    618,1616    621,1619    624,1622    627,1625    633,1631    636,1634    648,1645    651,1648    654,1650    717,1715    720,1718    723,1720    726,1723    750,1746    753,1748    768,1767    771,1770    774,1773    777,1776    780,1779    783,1782    786,1785"]

    # current mapping
    ref = ""
    len_ref = 0
    mapping = {}
    boxes = {}
    folder = "./out/"
    for i,line in enumerate(edge_file):

        if(i>1):
            line = line.rstrip("\n")
            data = line.split("\t")
            read1 = data[0]
            l1 = int(data[1])
            read2 = data[2]
            l2 = int(data[3])
            orientation = data[4]
            pos = [ tuple(int(a) for a in el.split(",")) for el in  data[5:] ]

            g.add_edge(read1,read2, nk = len(pos), ori= ori)

            if(not ref):
                ref = read1
                len_ref = l1

            if(read1 != ref):
                dash_rep(folder+ref,len_ref,mapping)

                mapping = {}
                ref = read1
                len_ref = l1
                break
                                    

            mapping[(read2,orientation)] = pos
            boxes[(read2,orientation)] = find_blocks(pos)

    #dash_rep(folder+ref,len_ref,mapping)
    block_rep(folder+ref, len_ref, boxes )
    nx.write_graphml(g,"./out.graphml")
    
#            same_iso = check_isoform(pos1,pos2,l1,l2)
            
            #if(same_iso):
            #g.add_edge(read1,read2, same_iso = same_iso, nk = len(pos1), ori= ori)
            
    

