# -*- coding: utf-8 -*-
"""
Created on Fri May  3 13:25:41 2019

@author: 45945
"""
import re
import collections

def hash_10kb(location,mod):
    result={}
    num_regions={}
    for i in range(1,22):
        num_regions[str(i)]=0
    with open(location) as loc:
        lines=loc.readlines()
        for line in lines[1:]:
            #line=line.replace("\r", " ")
            new=re.sub(' +|\r|\n|\t', ' ', line).strip()
            new=new.split()
            #new=line.split(" ")
            if new[1] == 'X': 
                new[1] = 20 
                print("a")
            elif new[1] == 'Y': new[1] = 21
            else: new[1] = int(new[1])
            temp=[int(new[1]),int((int(new[2])/mod))]
            if temp[1]>num_regions[str(temp[0])]:
                num_regions[str(temp[0])]=temp[1]
            result[int(new[0])]=temp
    
#    text_file = open("10kB_num_region_new.txt", "w")
#    for i in range(1,22):
#        s=str(num_regions[str(i)])+' '
#        text_file.write(s)
#    text_file.close()
    return result

def sort_by_key(dic):
    
    #dic = collections.OrderedDict(sorted(dic.items()))
    #print(dic)
    text_file = open("10kB_dic_new.txt", "w")
    for key in dic:
        s=str(key)+':'+str(dic[key])+"\n"
        text_file.write(s)
    text_file.close()
    return 

if __name__ == '__main__':
    mod=10000
    dic=hash_10kb("GATC.fends",mod)
    #sort_by_key(dic)