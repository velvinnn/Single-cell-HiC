# -*- coding: utf-8 -*-
"""
Created on Mon May  6 10:39:18 2019

@author: 45945
"""
import os
import shutil
import re
import numpy as np

def read_order(loc):
    file_loc=[]
    c=0
    with open(loc) as f:
        lines=f.readlines()
        for line in lines:
            num=line.split()
            num[1]=num[1].replace("Nagano","Raw")  #location for the raw data
            file_loc.append(num[1])
            c+=1
    print(c)
    return file_loc

def get_dic(location):
    dic={}
    max_pos=np.zeros(25)
    with open(location) as loc:
        lines=loc.readlines()
        for line in lines:
            line=re.findall(r'\d+',line)
            if line[1]=='X':
                line[1]=23
            if line[1]=='Y':
                line[1]=24
            dic[line[0]]=[int(line[1]),int(line[2])]   #chrome, pos
            if int(line[2])>max_pos[int(line[1])-1]:
                max_pos[int(line[1])-1]=int(line[2])
    text_file = open("10kB_num_region.txt", "w+")
    for i in range(0,24):
        s=str(int(max_pos[i]))+' '
        text_file.write(s)
    text_file.close()
    return dic


def map_reduce(location,dic,out):
    contact={}
    with open(location+"/adj") as loc:
        lines=loc.readlines()
        for line in lines[1:]:
            line=re.findall(r'\d+',line)
            key1,key2=dic[line[0]],dic[line[1]]
            if key1>key2: #smaller loc in the front
                key1,key2=key2,key1
            key=str(key1[0]) + ' ' + str(key1[1]) + ' ' + str(key2[0]) + ' ' + str(key2[1])
            if key not in contact:
                contact[key]=int(line[2])
            else:
                contact[key]+=int(line[2])
    text_file = open(out, "w+")
    for key in contact:
        s=str(key)+':'+str(contact[key])+"\n"
        text_file.write(s)
    text_file.close()
    return 
 

if __name__ == '__main__':
    file_loc=read_order(r"C:/study/HIC/cell_order.txt")
    dic=get_dic("10kB_dic.txt")
    for i in range(0,1):  
        #if os.path.isdir(file_loc[str(i)]): 
        map_reduce(file_loc[i],dic,"ordered_adj/"+str(i))
        print(i)
   