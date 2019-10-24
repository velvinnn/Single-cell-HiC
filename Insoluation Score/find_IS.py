# -*- coding: utf-8 -*-
"""
Created on Mon May  6 10:24:27 2019

@author: 45945
"""
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
import pickle
import re
import os
import math

def find_peaks(chromo,upper=200000):
    result=[]
    with open("Stadler_peaks.txt") as f:
        lines=f.readlines()
        for line in lines:
            num=line.split()
            if num[0]=="chr"+str(chromo):
                if int(num[1])<upper:
                    result.append(int(num[1]))
    return result


def get_log(contact):
    num=len(contact)
    for i in range(0,20):
        for j in range(0,i+1):
            contact[i,j]=math.log(contact[i,j]+1)
            contact[j,i]=contact[i,j]
    return contact

#def gen_full_map(adj_file,leng):
#    resolution=1            # for kubo, it is 10000
#    contact_map = np.zeros((leng,leng))
#    f = open(adj_file)
#    lines=f.readlines()
#    for line in lines:
#        lst = list(float(number) for number in (line.split()))
#        lst = [int(lst[0]/resolution), int(lst[1]/resolution), int(lst[2])]
#        if lst[0]>=1000 and lst[0]<leng+1000 and lst[1]>=1000 and lst[1]<leng+1000:
#            contact_map[lst[0]-1000][lst[1]-1000] += lst[2]
#            contact_map[lst[1]-1000][lst[0]-1000] += lst[2]
#    return contact_map

def get_full_map(adj_file, ch,upper=20000):
    # ch: chromosome number
    contact_map = np.zeros((int(upper),int(upper)))

    f = open(adj_file)
    for line in f:
        lst = re.findall(r'\d+',line)
        lst = [int(lst[0]), int(lst[1]), int(lst[2]), int(lst[3]), float(lst[4])]
        loc1,loc2=int(lst[1]),int(lst[3])
        if lst[0] == ch and lst[2] == ch and loc1<upper and loc2<upper:
            contact_map[loc1][loc2] += lst[4]
            if lst[1] != lst[3]:
                contact_map[loc2][loc1]+= lst[4]
    # np.save('test.npy', contact_map)
#    contact_map=get_log(contact_map+1)
#    contact_map=dev_avg(contact_map,leng)
    return contact_map


def find_sum(contact_map,box_length,peak):#,adj_file,target_name):

#    for i in range(20,total_length-box_length):
#        i=i-1
        i
        count_sum=np.sum(contact_map[i-box_length:i,i:i+box_length])
        s1,s2=0,0
#        count_sum[index]=np.sum(contact_map[i-box_length:i,i:i+box_length])
        s1=np.sum(contact_map[i-box_length:i,i-box_length:i])/2
        s2=np.sum(contact_map[i:i+box_length,i:i+box_length])/2
#            if count_sum[index]>0:            
#                count_sum[index]=float(count_sum[index]*2)/(s1+s2+count_sum[index])
                
#        index+=1
            
            #way 2 of calculating IS 
#    count_sum=list(count_sum/np.mean(count_sum))
#    count_sum=np.log2(np.add(count_sum,0.1))3
#    print(len(count_sum))
        if count_sum==0:
            return 0
        else: 
            return count_sum*2/(count_sum+s1+s2)



if __name__ == '__main__':
    num_region=[19719 ,18174, 15959 ,15563 ,15253 ,14951 ,15252 ,13173, 12407, 12999, 12184 ,12125 ,12028 ,12519 ,10349, 9831, 9527, 9077 ,6134 ,16664, 290, 0, 0, 0 ]
    box_length=20
    for ch in range(1,20):
        print(ch)
        peaks=find_peaks(ch,upper=200000)
        result=np.zeros((len(peaks),1171))
        for j in range(1171):
            contact_map=get_full_map("C:/study/HIC/ordered_adj/"+str(j), ch,upper=20000)
            for i in range(len(peaks)):
                result[i,j]=find_sum(contact_map,20,peaks[i])
    #            np.save("4_stages_contact.npy",contact_map)
    #            target_name="4_stages/chr"+str(ch)+"_stage"+str(i)+"_IS.npy"
    #            count_sum=find_sum(contact_map,box_length,total_length,adj_file,target_name)   
        np.save("IS/"+str(ch),result)

