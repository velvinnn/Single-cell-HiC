# This code aims to find the number of clusters IS.
from sklearn.cluster import KMeans
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial, Legendre


def generate_labels(X,clusters_num):
    kmeans = KMeans(n_clusters=clusters_num, random_state=0).fit(X)
    score=kmeans.score(X)
    return kmeans.labels_,score

def legend(x,order):
    coefficient=Legendre.basis(order).convert(kind=Polynomial)
    result=[]
    for i in range(0,len(x)):
        result.append(coefficient(x[i]))
    return result

def y_fit(x,coefficient):
    y=np.add(np.zeros(len(x)),coefficient[0])
    for i in range (1,len(coefficient)):
        y+=np.multiply(legend(x,i),coefficient[i])
    return y

def plot(y,lable):
    plt.figure(1)
    plt.plot(range(0,1171),y,color=lable)
    return


def find_cluster_num(X,test_range):
    score=[]
    for k in range(2,test_range):
         score.append(KMeans(n_clusters=k, random_state=0).fit(X).score(X))
    plt.figure(2)
    plt.plot(range(2,test_range),score)
    plt.show()
    return

if __name__ == '__main__':  
    lable=['red','yellow','green','blue','black']
    clusters_num=4
    window=10
    end_cell=51
    
    coefficient=np.load("chr1_coefficient.npy")[0:end_cell+1,:]
    count=int(np.floor(1171/window))
    y_fitted=np.zeros((end_cell+1,count))
    y_total=np.zeros((end_cell+1,1171))
    for i in range(0,end_cell+1):
        y_total[i,:]=y_fit(range(0,1571),coefficient[i,:])[200:1371]
        for j in range(0,count):
            y_fitted[i,j]=np.sum(y_total[i,j*window:(j+1)*window])
            
#    find_cluster_num(y_fitted,12)
    groups,score=generate_labels(y_fitted,clusters_num)
    for i in range(0,end_cell+1):
        plot(y_total[i,:],lable[groups[i]])
    plt.show()
    
    
