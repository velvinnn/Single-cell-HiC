# Use polynomial fit to smooth the IS curve
# Use kmeans to find the clusters

#This code clusters chr 1-4 as an exmaple

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import filtfilt
from sklearn.cluster import KMeans

import re

def smooth_fit(y):
    n = 15  # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    y=list(y)
    y=y[-200:-1]+y+y[0:200]
    y = list(filtfilt(b,a,y))
    x=range(0,len(y))
    val=np.polyfit(x,y,12)
    y_fit=list(val[0]*pow(a,12)+val[1]*pow(a,11)+val[2]*pow(a,10)+val[3]*pow(a,9)+val[4]*pow(a,8)+val[5]*pow(a,7)
        +val[6]*pow(a,6)+val[7]*pow(a,5)+val[8]*pow(a,4) +val[9]*pow(a,3)
        +val[10]*pow(a,2)+val[11]*pow(a,1)+val[12] for a in x)
    return y_fit[200:1371]


def find_cluster_num(X,test_range,chromo):
    score=[]
    for k in range(2,test_range):
         score.append(KMeans(n_clusters=k, random_state=0).fit(X).score(X))
    plt.figure(2)
    plt.plot(range(2,test_range),score)
    plt.savefig("score_chr"+str(chromo)+".png")
    plt.show()
    return

def generate_labels(X,clusters_num):
    kmeans = KMeans(n_clusters=clusters_num, random_state=0).fit(X)
    score=kmeans.score(X)
    return kmeans.labels_,score

def plot(y,lable):
    plt.plot(range(0,1171),y,color=lable)
    return


if __name__ == '__main__': 
    peak_count=[460,440,325,326,344,351,307,286,311,290,295,247,283,238,222,219,204,201,138]
    color=["red","green","yellow","blue","black"]
    group_order=[[0,3,2,1],[0,2,3,2],[0,2,1,3],[1,0,0,1],[1,2,3,0]]
    group_order=np.array(group_order)
    for chromo in range(2,3):
        name="Chr"+str(chromo)
        peak_num=peak_count[chromo-1]
        IS=np.load(name+"_IS.npy")
        Y=np.zeros((peak_num,13))
        y_fit=np.zeros((peak_num,1171))
        
        for i in range(0,peak_num):
            y_fit[i,:]=smooth_fit(IS[:,i])
            
        
        lable=['red','yellow','green','blue','black']
        clusters_num=4
        y_train=y_fit
        for i in range(0,peak_num):
            y_train[i,:]=np.add(y_train[i,:],-np.mean(y_train[i,:]))
            
#        find_cluster_num(y_train,16,chromo)
#        groups=np.load("groups_"+str(chromo)+".npy")
        groups,score=generate_labels(y_train,clusters_num)
        np.save("groups_2.npy",groups)
#        print(groups)
        for i in range(0,peak_num):
            plt.figure(groups[i])
            plot(y_fit[i,:],lable[groups[i]])
            plt.figure(100)
            plot(y_fit[i,:],lable[groups[i]])
        
        for i in range(0,clusters_num):
            plt.figure(i)
#            y=[]
#            for j in range(0,1171):
#                y.append(np.mean([np.mean(y_fit[k,j]) for k in range(0,peak_num) if groups[k]==i]))
#    
#            x=range(0,len(y))
#            y1=np.add(y,1.96*np.std(y))
#            y2=np.add(y,-1.96*np.std(y))
            plt.ylim(-0.7,0.7)
            plt.axvline(x=279)
            plt.axvline(x=582)
            plt.axvline(x=844)
            plt.title(name+"_group"+str(i))
            plt.savefig(name+"_group"+str(i))
#            plt.plot(x,y,c=color[chromo-1],label=str(chromo))
#            plt.plot(x,y1,c=color[chromo-1],linestyle=":")
#            plt.plot(x,y2,c=color[chromo-1],linestyle=":")
            plt.legend()
        print(chromo)
            #plt.show()
#        plt.figure(100)
#        plt.axvline(x=279)
#        plt.axvline(x=582)
#        plt.axvline(x=844)
#        plt.title(name)
#        #plt.savefig(name)
#        plt.show()

