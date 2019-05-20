#coding=utf-8

import sys
import networkx as nx

global error_rate, coverage
error_rate = 0.20
coverage = 0.70

# def check_isoform(pos1,pos2,l1,l2):
	
# 	blocks = []
# 	current_block = []
# 	last = None
# 	for i in  range(len(pos1)):
# 		p1,p2 = pos1[i],pos2[i]
# 		if(last):
# 			dif1 = abs(p1 - last[0])
# 			dif2 = abs(p2 - last[1])
# 			if(dif1!=0 and dif2!=0):

# 				# if the difference between current and last pos for both read
# 				# is not too different, we increase the bloc size
# 				# print(p1,p2,dif1,dif2, abs(1 - float(dif1) / dif2))
# 				if( (max(dif1,dif2) <= 6 and  abs(dif1-dif2)<=2) or (  abs(1 - float(dif1) / dif2)) <=  error_rate ):
# 					current_block.append( (p1,p2) )
# 				elif(current_block):
# 					blocks.append(current_block)
# 					current_block = [(p1,p2)]
# 		else:
# 			current_block.append((p1,p2))
# 		last = (p1,p2)
# 		#print(p1,p2)
# 	blocks.append(current_block)
# 	# print(len(blocks))
# 	# print("\n".join(str(el) for el in blocks))
# 	# simple method: max over minsize
# 	max_bl=  max(blocks, key= lambda x: len(x))
# 	cover = 0
# 	if(l1<=l2):
# 		r = 0
# 	else:
# 		r =1
# 	cover = abs(max_bl[-1][r] - max_bl[0][r]) + 16
# 	# print(cover,l1,l2)
# 	if( float(cover) / min(l1,l2) >= coverage):
# 		return(True)
	
# 	return(False)

def check_isoform(pos1,pos2,l1,l2):
	cover1 = float(pos1[0]-pos1[-1])/l1
	cover2 = float(pos2[0]-pos2[-1])/l2
	#if( abs(1 - float(l1) / l2 ) <=  error_rate):# and
	if(cover1>=coverage or cover2 >= coverage):
		return(True)
	return(False)

g = nx.Graph()
with open(sys.argv[1]) as edge_file:
	#edge_file = ["sequence1\t40\tsequence2\t42\t0\t1,1\t2,2\t4,5\t6,8\t30,39\t33,41"]
	#edge_file = ["ch115_read2745_template_pass_BYK_CB_ONT_1_FAF04998_A	804	ch273_read4762_template_pass_BYK_CB_ONT_1_FAF04998_A	3252	0	9,1007	12,1010	18,1015	21,1018	24,1021	75,1077	78,1080	81,1084	96,1098	99,1100	102,1103	105,1106	108,1109	123,1122	126,1124	132,1130	135,1133	168,1162	174,1169	177,1172	180,1175	183,1178	186,1180	195,1193	198,1197	201,1199	204,1202	207,1205	213,1211	216,1214	219,1217	222,1220	225,1223	228,1227	231,1230	234,1233	237,1237	240,1240	243,1243	246,1246	252,1252	267,1268	270,1270	288,1288	291,1291	294,1294	297,1297	312,1309	315,1313	318,1316	321,1319	324,1322	327,1325	330,1328	333,1331	336,1334	339,1337	342,1340	345,1343	348,1346	411,1406	414,1408	483,1486	486,1488	489,1491	492,1494	495,1497	498,1500	501,1503	504,1506	555,1556	558,1559	561,1562	564,1565	567,1568	606,1604	609,1607	612,1610	615,1613	618,1616	621,1619	624,1622	627,1625	633,1631	636,1634	648,1645	651,1648	654,1650	717,1715	720,1718	723,1720	726,1723	750,1746	753,1748	768,1767	771,1770	774,1773	777,1776	780,1779	783,1782	786,1785"]
	for i,line in enumerate(edge_file):

		if(i>1):
			line = line.rstrip("\n")
			data = line.split("\t")
			read1 = data[0]
			l1 = int(data[1])
			read2 = data[2]
			l2 = int(data[3])
			ori = data[4]
			pos1,pos2 = zip(* [ [ int(a) for a in el.split(",")] for el in  data[5:] ])
			same_iso = check_isoform(pos1,pos2,l1,l2)

			g.add_node(read1)
			g.add_node(read2)
			#if(same_iso):
			g.add_edge(read1,read2, same_iso = same_iso, nk = len(pos1), ori= ori)

nx.write_graphml(g,"./out.graphml")