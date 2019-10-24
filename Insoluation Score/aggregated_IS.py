#This code get IS from aggregating data from all cells.
# Aimed to compare bulk HiC and sc HiC
import re
import numpy as np
def gen_full_map(adj_file, ch,leng):
    # ch: chromosome number
    contact_map = np.zeros((int(leng),int(leng)))

    f = open(adj_file)
    for line in f:
        lst = re.findall(r'\d+',line)
        lst = [int(lst[0]), int(lst[1]), int(lst[2]), int(lst[3]), float(lst[4])]
        loc1,loc2=int(lst[1])-1,int(lst[3])-1
        if lst[0] == ch and lst[2] == ch:
            contact_map[loc1][loc2] += lst[4]
            if lst[1] != lst[3]:
                contact_map[loc2][loc1]+= lst[4]
    # np.save('test.npy', contact_map)
    
    return contact_map

def find_sum(contact_map):
    box_length=20
    for i in range(3000,19680):
        i=i-1
        count_sum.append(0)
        s1,s2=0,0
        count_sum[index]=np.sum(contact_map[i-box_length:i,i:i+box_length])
            #way 1 of calculating IS:
#            s1=np.sum(contact_map[i-box_length:i,i-box_length:i])/2
#            s2=np.sum(contact_map[i:i+box_length,i:i+box_length])/2
#            if count_sum[index]>0:            
#                count_sum[index]=float(count_sum[index]*2)/(s1+s2+count_sum[index])
                
            
            #way 2 of calculating IS 
    count_sum=count_sum/np.mean(count_sum)
    count_sum=np.log2(count_sum+0.1)            
        #show_plot(count_sum,orderi)
    
    return count_sum

if __name__ == '__main__':

     
    total_cell=1171
    leng=19719
    result=np.zeros((leng,leng))
    for cell in range(0,1171):
        result+=gen_full_map(str(cell), 1,19719)
        print(cell)
    count_sum=find_sum(result)
    np.save("count_sum.npy",count_sum)


