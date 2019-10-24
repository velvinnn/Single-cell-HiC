# -*- coding: utf-8 -*-
"""
Created on Sun May 26 16:52:29 2019

@author: 45945
"""
import re
import numpy as np
import math

def bin_count():
    with open("10kB_num_region.txt") as loc:
        lines=loc.readlines()
        for line in lines:
            num=re.findall(r'\d+',line)
    return num

def read_adj(cell_index,chromosome,bin_count):  #only 21 chromosome
    result=np.zeros((int(bin_count[chromosome-1]),int(bin_count[chromosome-1])),dtype = int)
    with open("ordered_adj/"+str(index)) as loc:
        lines=loc.readlines()
        for line in lines:
            num=re.findall(r'\d+',line)
            pos1,pos2=int(num[1])-1,int(num[3])-1
            #print(pos1,pos2)
            if num[0]==num[2] and int(num[0])==chromosome:
                result[pos1,pos2]+=int(num[4])
                result[pos2,pos1]+=int(num[4])
#    target="ordered_adj_matrix/cell"+str(index)+"chromo"+str(i)+".txt"
#    text_file = open(target, "w+")
#    for j in range(0,int(bin_count[i-1])):
#        for k in range(0,j+1):
#             text_file.write(str(int(result[k,j])))
#             text_file.write(" ")
#        text_file.write("\n")   
#    text_file.close()
    return result


bin_count=bin_count()
for i in range(1,1171):
    read_adj(i,bin_count)
