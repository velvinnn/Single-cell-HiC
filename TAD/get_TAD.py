
import numpy as np
import matplotlib.pyplot as plt

def find_TAD_num(chromosome, CTCF_loc="C:/Users/1/Desktop/IS/Stadler_peaks.txt"):
    # for this function, you enter the location of CTCF binding sites.
    CTCF_binding_site=[]
    with open(CTCF_loc) as loc:
        lines=loc.readlines()
        for line in lines:
            read=line.split()
            if chromosome == 20:
                if read[0]=="chrX":
                    CTCF_binding_site.append(int(read[1]))  
            elif chromosome == 21:
                if read[0]=="chrY":
                    CTCF_binding_site.append(int(read[1]))  
            else:
                if read[0]=="chr"+str(chromosome):
                    CTCF_binding_site.append(int(read[1]))  
    return len(CTCF_binding_site)

def find_drops(IS, TAD_num,window_size):
    # We select possible TAD boundaries with local minimal insulation score
    
    # To avoid the situation that the insulation score is a local minimal  
    # within a highly contacted region, we set a threshold 
    # By screening along different thresholds, we can find the number of TAD closest 
    # to the number we have from CTCF binding sites.
    # you could change the TAD number 
    
    # To avoid the situation that the IS is low in the surrounding region,
    # we set a window size so that we only select at most one peak within this window length along the chromsome
    min_difference=1000
    result_peaks=[]
    for percent in range(1,30):
        threshold=np.percentile(IS, percent)
        peaks=[np.argmin(IS[window_size*i:window_size*(i+1)])+window_size*i 
               for i in range(1,int(len(IS)/window_size)) 
               if min(IS[window_size*i:window_size*(i+1)])<threshold]
        
        if len(IS) % window_size>1:
            if min(IS[window_size*int(len(IS)/window_size):len(IS)-1]) < threshold:
                peaks.append(np.argmin((IS[window_size*int(len(IS)/window_size):len(IS)-1]))
                    +window_size*int(len(IS)/window_size))
        
        if abs(len(peaks)-TAD_num)<min_difference:
            min_difference=abs(len(peaks)-TAD_num)
            result_peaks=peaks
            
    return np.add(result_peaks,300)

for chromosome in range(1,21):
    # input the location for IS file
    loc="C:/Users/1/Bioinformatics/TAD/IS/"+str(chromosome)+".npy"
    result=np.load(loc)
    if chromosome==20:
        result=result[200:len(result)]
    
    # you find the number of TAD we have from the CTCF data
    CTCF_binding_site=find_TAD_num(chromosome)
    #print(len(CTCF_binding_site))
    # find a similar number of TAD from IS 
    a=find_drops(result, CTCF_binding_site,30)  
    if chromosome==20:
        a=(np.add(a,200))*10000
    else:
        a=a*10000
    #np.save("C:/Users/1/Desktop/IS/possible_TAD_boudaries/TAD_boundaries"+str(chromosome)+".npy",a)
    with open('C:/Users/1/Bioinformatics/TAD/IS/TAD_from_IS_txt/'+str(chromosome)+'.txt', 'w') as f:
        for item in a:
            f.write("%s\n" % item)
    print("chr",chromosome,CTCF_binding_site,len(a))
    #print(CTCF_binding_site[0:10],a[0:10])
    #plt.figure(1)
    #plt.scatter(a[0:6],np.ones(6),c="r")
    #plt.scatter(CTCF_binding_site[0:6],np.ones(6),c="b")
