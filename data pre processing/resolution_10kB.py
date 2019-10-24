# -*- coding: utf-8 -*-
"""
Created on Fri May  3 11:34:39 2019

@author: 45945
"""
import re
import os.path

def get_dic(location):
    with open(location) as loc:
        lines=loc.readlines()
        for line in lines:
            line=re.findall(r'\d+',line)
            dic[line[0]]=[int(line[1]),int(line[2])]   #chrome, pos
    return dic


def map_reduce(location,dic,dir):
    contact={}
    with open(location) as loc:
        lines=loc.readlines()
        for line in lines[1:]:
            line=re.findall(r'\d+',line)
            key1,key2=dic[line[0]],dic[line[1]]
            if key1>key2 :
                key1,key2=key2,key1
            key=str(key1[0]) + ' ' + str(key1[1]) + ' ' + str(key2[0]) + ' ' + str(key2[1])
            if key not in contact:
                contact[key]=int(line[2])
            else:
                contact[key]+=int(line[2])
    out=dir+"/10kB_adj.txt"
    text_file = open(out, "w+")
    for key in contact:
        s=str(key)+':'+str(contact[key])+"\n"
        text_file.write(s)
    text_file.close()
    return 

def gen_new_adj(dir, dic,index):
    # build kekeys of GATC dictionary (for later binary search)
    # Traversing all data files
    for i in range(1,488):
        open_file = dir+'/1CDX'+str(index)+'.'+str(i)
        
        print(open_file)
#        exists = os.path.isfile(open_file)
#        if exists:
        #print('Changing resolution on file: %s'%open_file)
        if os.path.isdir(open_file): 
            map_reduce(open_file+'/adj',dic,open_file)
    return

if __name__ == '__main__':
   dic=get_dic("10kB_dic.txt")
   for i in range(1,5):
       dir=r"C:/study/HIC/schic_hyb_1CDX"+str(i)+"_adj_files"
       gen_new_adj(dir, dic,i)