import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import filtfilt
def find_CTCF_peaks(chromo):
    result=[]
    with open("C:/study/HIC/BulkIS/Stadler_peaks.txt") as f:
        lines=f.readlines()
        for line in lines:
            num=line.split()
            if num[0]=="chr"+str(chromo):
                result.append(int(num[1]))
    return result
def get_log(contact):
    num=len(contact)
    for i in range(0,20):
        for j in range(0,i+1):
            contact[i,j]=math.log(contact[i,j]+1)
            contact[j,i]=contact[i,j]
    return contact

def get_full_map(adj_file, lower,upper):
    # ch: chromosome number
    leng = upper-lower
    contact_map = np.zeros((int(leng),int(leng)))

    f = open(adj_file)
    for line in f:
        lst = re.findall(r'\d+',line)
        lst = [int(lst[0]), int(lst[1]), int(lst[2]), int(lst[3]), float(lst[4])]
        loc1,loc2=int(lst[1]),int(lst[3])
        if lst[0] == ch and lst[2] == ch and loc1 <upper and loc1>=lower and loc2 <upper and loc2>=lower:
            contact_map[loc1][loc2] += lst[4]
            if lst[1] != lst[3]:
                contact_map[loc2][loc1]+= lst[4]
    # np.save('test.npy', contact_map)
    contact_map=get_log(contact_map)
#    contact_map=dev_avg(contact_map,leng)
    return contact_map

def find_peaks(y):
    #finding local minimun with IS < 0
    peaks = np.where((y[1:-1] < y[0:-2]) * (y[1:-1] < y[2:]))[0] + 1
    result=[]
    for peak in peaks:
        if y[peak]<0:
            result+=[peak]
    return peaks

def smooth(y):
    n = 3  # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    y=list(y)
    y = list(filtfilt(b,a,y))
    return y

def normalize_cell(IS,TAD_num):
    # smooth IS in every TAD across cell lines
    # for every cell, scale all the IS for 460 TADs to normal distribution
    for j in range(TAD_num):
        IS[:,j]=smooth(IS[:,j])
        IS[:,j]=(IS[:,j]-np.mean(IS[:,j]))/np.std(IS[:,j])
        IS[:,j]=smooth(IS[:,j])
#        for cell in range(0,1171):
#        IS[cell,:]=(IS[cell,:]-np.mean(IS[cell,:]))/np.std(IS[cell,:])
    return IS


if __name__ == '__main__': 
    box=10
    for chromo in range(1,20):
        # piling-up for IS
        groups=list(np.load('C:/study/HIC/Bulk/clustering_IS/groups/groups_'+str(chromo)+'.npy'))
        IS=np.load('C:/study/HIC/Bulk/IS/result/chr'+str(chromo)+'_IS_peaks_log.npy')
        IS=normalize_cell(IS,np.shape(IS)[1])
        mean_IS=np.load('C:/study/HIC/Bulk/clustering_IS/mean_curve/mean_curve_chr'+str(chromo)+'.npy')
        for group in range(0,4):
            TAD_index=[i for i,x in enumerate(groups) if x==group]
            IS_peaks=find_peaks(mean_IS[group,:])
            for peak in IS_peaks:
                aggregated=np.zeros(2*box+1)
                for i in TAD_index:
                    group_IS=(IS[peak-box:peak+box+1,i])
                    aggregated+=(group_IS)
                findex=group+100*chromo
                plt.figure(findex)
                plt.title("chr"+str(chromo)+", cluster "+str(group))
                plt.plot(range(2*box+1),aggregated,label="cell"+str(peak)+"_ chr"+str(chromo))
                plt.legend()
                plt.savefig("pile_up/chr"+str(chromo)+"_ cluster "+str(group))
#      for chromo in range(1,2):
#        # piling-up for IS
#        groups=list(np.load('C:/study/HIC/Bulk/clustering_IS/groups/groups_'+str(chromo)+'.npy'))
#        mean_IS=np.load('C:/study/HIC/Bulk/clustering_IS/mean_curve/mean_curve_chr'+str(chromo)+'.npy')
#        CTCF_peaks=find_CTCF_peaks(chromo)
#        
#        for group in range(0,4):
#            TAD_index=[i for i,x in enumerate(groups) if x==group]
#            IS_peaks=find_peaks(mean_IS[group,:])
#            for cells in IS_peaks:
#                pile=np.zeros((2*box+1,2*box+1))
#                for i in TAD_index:
#                    adj_file="C:/study/HIC/ordered_adj/"+
#                    aggregated+=get_full_map(adj_file, lower,upper)
#                findex=group+100*chromo
#                plt.figure(findex)
#                plt.title("chr"+str(chromo)+", cluster "+str(group))
#                plt.plot(range(2*box+1),aggregated,label="cell"+str(peak)+"_ chr"+str(chromo))
#                plt.legend()
#                plt.savefig("pile_up/chr"+str(chromo)+"_ cluster "+str(group))     