#we skip the first 300 loci and the need to add the bin size as 20  -> +320=true locus
import numpy as np
import matplotlib.pyplot as plt
import re



def read_IS(loc):
    with open(loc) as f:
        lines=f.readlines()
        for line in lines:
            result=list(float(number) for number in (line.split()))
    return result

def find_peaks(bound, window_size,result):
    peaks=[]
    leng=int(np.floor(len(result)/window_size))
    for i in range(0,leng):
        a=min(result[i*window_size:(i+1)*window_size])
        if a<bound:
            peaks.append(int(result[i*window_size:(i+1)*window_size].index(a)+i*window_size))
    if leng*window_size!=len(result):
        a=min(result[leng*window_size:len(result)])
        if a<bound:
                peaks.append(int(result[leng*window_size:len(result)].index(a)+leng*window_size))
    return peaks

def find_nearest(array, value):
    array = np.asarray(array)
    return  min(np.abs(np.add(array,-value)))

def read_peaks(chromo):
    loc="Stadler_peaks.txt"
    result=[]
    with open(loc) as f:
        lines=f.readlines()
        for line in lines:
            num=list(re.findall(r'\d+',line))
            if int(num[0])==chromo:
                result.append(int(num[1]))
    print(len(result))
    return result


if __name__ == '__main__':
    num_region=[[300,19537],[305,18202],[300,15993],[305,15625],[300,15173],[305,14958],[300,14534],[300,12930],[300,12449],[310,13059],[310,12198],[300,12002],[300,12032],[300,12480],[305,10394],[300,9810]]
    num_region=np.array(num_region)
    groups=[]
    groups_peaks=[]
    for i in range(1,17):
        mm10=read_IS("mm10_IS_"+str(i)+".txt")        
        result=read_IS("microC_chr"+str(i)+"_IS.txt")
        if len(mm10)<len(result)+num_region[i-1,0]-300:
            result=result[0:len(mm10)+300-num_region[i-1,0]]
        else:
            mm10=mm10[0:len(result)+num_region[i-1,0]-300]
        bound=np.percentile(result, 12)
        window_size=10
        peaks=find_peaks(bound, window_size,result)
        peaks=np.add(peaks,300-num_region[i-1,0])
        print(len(peaks))

    #    plt.figure(1)
    #    plt.plot(range(len(mm10)),mm10)
    #    for i in peaks:
    #        plt.figure(1)
    #        plt.axvline(x=i,c="red") 
        peaks_mm10=[mm10[k] for k in peaks]
#        print(np.mean(mm10))
#        print(np.mean(peaks_mm10))
        groups.append(mm10)
        groups_peaks.append(peaks_mm10)
        print(i)
    plt.figure()
   
    plt.boxplot(groups,showfliers=False,positions=np.array(range(len(groups)))-0.15+1, widths=0.15)
    plt.boxplot(groups_peaks,showfliers=False,positions=np.array(range(len(groups_peaks)))+1, widths=0.15)
    ax.xlable("chromosome ID")
#    plt.ylable("IS")    
#    simi=[]
#    for i in range(0,len(peaks)):
#        simi.append(find_nearest(result_sc, peaks[i]))
#    simi=np.sort(simi)
#    # on average, thedistance between two points should be 40
#    print(np.mean(simi))
#    print(np.percentile(simi, 50))
    #[305,18201] for chromo2
    #[300,19537] for chromo1
    #[300,15993] for chromo3
    
#    
#    
#range for mm10:
#1
#302 19533
#2
#310 18200
#3
#303 15994
#4
#307 15623
#5
#310 15172
#6
#341 14956
#7
#323 335