#-*-coding:utf-8 -*-
#2019-10-22
#find duplicate umi
from sys import argv
import os
import re
# import multiprocessing
# import numpy as np
import math
import random
# from scipy.stats import pearsonr
input1,input2,input3,output1,output2=argv[1:]
def fq_f(file):
    while True:
        try:
            header=next(file).strip()
            seq=next(file).strip()
            next(file)
            next(file)
            yield header,seq
        except StopIteration:
            break
with open(input1,"r") as fq:
    match_dict={}
    for h,s in fq_f(fq):
        index_n=h.split()[0][1:]
        match_dict[index_n]=s[:3]
with open(input2,"r") as sam1, open(input3,"r") as sam2,open(output1,"w") as out,open(output2,"w") as out2:
    # for sam_line in sam:
    #     if sam_line[0]!="@":
    #         sam_list=sam_line.strip().split()
    #         out.write(sam_list[0]+"\t"+sam_list[1]+"\t"+sam_list[2]+"\t"+sam_list[3]+"\t"+sam_list[4]+"\t"+sam_list[5]+"\t"+sam_list[9]+"\t"+match_dict[sam_list[0]]+"\n")
    posi_umi_dict={}
    read1_dict={}
    read2_dict={}
    for sam_line in sam1:
        if sam_line[0]!="@":
            line_list=sam_line.strip().split()
            index=line_list[0]
            posi="".join(line_list[1]+line_list[2]+line_list[3]+line_list[5]+match_dict[index])
            score=line_list[4]
            # if posi in posi_index_dict.keys():
            #     posi_index_dict[posi].append(index)
            # else:
            #      posi_index_dict[posi]=[index]
            # neiceng={}
            neiceng=(index,score)
            if posi in posi_umi_dict.keys():

                posi_umi_dict[posi].append(neiceng)
            else:
                posi_umi_dict[posi]=[neiceng]
    sam1.seek(0)
    for r_line1 in sam1:
        if r_line1[0]!="@":
            r_line1_list=r_line1.strip().split()
            read1_dict[r_line1_list[0]]=r_line1.strip()
        else:
            out.write(r_line1)
    for r_line2 in sam2:
        if r_line2[0]!="@":
            r_line2_list=r_line2.strip().split()
            read2_dict[r_line2_list[0]]=r_line2.strip()  
        else:
            out2.write(r_line2) 
    for po in posi_umi_dict.keys():
        per_po=posi_umi_dict[po]

        per_po_sort=sorted(per_po,key=lambda per_po:per_po[1],reverse=True)
        max_score=per_po_sort[0][0]
        if max_score in read1_dict.keys() and max_score in read2_dict.keys():
            out.write(read1_dict[max_score]+"\n")
            out2.write(read2_dict[max_score]+"\n")