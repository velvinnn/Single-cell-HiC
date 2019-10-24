from sklearn.gaussian_process import GaussianProcessRegressor
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import filtfilt
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel,ExpSineSquared,RBF
from sklearn.model_selection import train_test_split
import math
from sklearn.decomposition import KernelPCA



def smooth(xx,fit=False,fig=1):
    x=xx.copy()
    n = 15  # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    y=list(x)
    y=y[-200:-1]+y+y[0:200]
    if fit:
        y = list(filtfilt(b,a,y,axis=0))
#        plt.figure(i)
#        plt.plot(range(1171),y[200:1371],label="filter")
        x=range(0,len(y))
        val=np.polyfit(x,y,12)
        y_fit=list(val[0]*pow(a,12)+val[1]*pow(a,11)+val[2]*pow(a,10)+val[3]*pow(a,9)+val[4]*pow(a,8)+val[5]*pow(a,7)
        +val[6]*pow(a,6)+val[7]*pow(a,5)+val[8]*pow(a,4) +val[9]*pow(a,3)
        +val[10]*pow(a,2)+val[11]*pow(a,1)+val[12] for a in x)
        return y_fit[200:1371]
    return x

def dot_product(xx,yy):
    dis=0
    x=xx.copy()-np.mean(xx)
    y=yy.copy()-np.mean(yy)
    for i in range(len(x)-1):
        dot_x=[1,x[i+1]-x[i]]
        dot_y=[1,y[i+1]-y[i]]
        dis+=np.dot(dot_x,dot_y)
   # print(dis/math.sqrt(np.dot(xx,xx)*np.dot(yy,yy)))
    return dis/math.sqrt(np.dot(xx,xx)*np.dot(yy,yy))

def find_distance(IS):
    size=np.shape(IS)[1]
    dis_matrix=np.zeros((size,size))
    for i in range(size):
        dis_matrix[i,i]=1
        for j in range(0,i):
            dis_matrix[i,j]=dot_product(IS[:,i],IS[:,j])
            dis_matrix[j,i]=dis_matrix[i,j]
        print(i)
    return dis_matrix


if __name__ == '__main__': 
    peak_count=[460,440,325,326,344,351,307,286,311,290,295,247,283,238,222,219,204,201,138]
    color=["red","green","yellow","blue","black"]
    for chromo in range(8,19):
        iss=np.load("C:/study/HIC/Bulk/IS/result/chr"+str(chromo)+"_IS_peaks_log.npy")
        groups=np.load("C:/study/HIC/Bulk/clustering_IS/groups/groups_"+str(chromo)+".npy")
        result=np.zeros((1171,peak_count[chromo-1]))
        for i in range(peak_count[chromo-1]):
            IS=list(iss[:,i])
            result[:,i]=(smooth(list(IS),True))
            print(i)
        np.save("poly_fit_log_IS/"+str(chromo)+".npy",result)
#    for chromo in range(2,10):
#        IS=np.load("C:/study/HIC/Bulk/IS/IS/"+str(chromo)+".npy")
#        for i in range(np.shape(IS)[0]):
#            IS[i,:]=(smooth(list(IS[i,:]),True))
#            print(i)
#        np.save("poly_fit_IS/"+str(chromo)+".npy",IS)
    end=1
    group=[]
    IS=np.zeros((1171,np.sum(peak_count[0:end])))
    index=0
    for chromo in range(1,end+1):
        temp=np.load("C:/study/HIC/Bulk/clustering_IS/groups/groups_"+str(chromo)+".npy")
        if chromo==1:
            group=temp
        else:
            group=list(group)+list(temp)
        temps=np.load("poly_fit_log_IS/"+str(chromo)+".npy").T
        IS[:,index:index+peak_count[chromo-1]]=temps.T
        index+=peak_count[chromo-1]
    color=["black","red","yellow","blue"]
    colors=[color[i] for i in group]
#      
    chromo=1
    dis=find_distance(np.load("poly_fit_IS/"+str(chromo)+".npy").T)#IS.T)
    np.save("dis"+str(chromo)+".npy",dis)
    
    chromo=1
    color=["black","red","yellow","blue"]
    colors=[color[i] for i in np.load("C:/study/HIC/Bulk/clustering_IS/groups/groups_"+str(chromo)+".npy")]
#    dis=np.load("dis"+str(chromo)+".npy")
#    IS=np.load("poly_fit_log_IS/"+str(chromo)+".npy")
#    for i in range(222):
#        IS[:,i]=IS[:,i]-np.mean(IS[:,i])
#        dis[i,i]=math.sqrt(np.dot(IS[:,i],IS[:,i]))
    transformer = KernelPCA(n_components=2, kernel="precomputed",random_state=0)
    X_transformed = transformer.fit_transform(np.load("poly_fit_IS/"+str(chromo)+".npy"))
    plt.figure(1)
    plt.scatter(X_transformed[:,0],X_transformed[:,1],color=colors)
    plt.title("PCA of chromosome 1")
