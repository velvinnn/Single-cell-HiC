import numpy as np


def find_peaks(file, output, resolution=10000, window=50000, threshold=70):
    # Find all local maximum with reads which are also greater than threshold value
    # If several local maximums are too close to each other (<= window size), only keep one highest peak
    all_peaks = []
    for line in open(file):
        lst = line.strip().split()
        chromosome = lst[0]
        #all_chromosomes = {f'chr{elm}' for elm in list(range(1, 20)) + ['X']}
        all_chromosomes = {f'chr{elm}' for elm in list(range(1, 2))}
        if chromosome not in all_chromosomes:
            continue
        position, read = int(lst[1]) // resolution, float(lst[-1])
        new_lst = [chromosome, position, read]
        all_peaks.append(new_lst)

    window = window // resolution
    candidates = []
    # Find all local maximum values which are greater than threshold    
    for i in range(1, len(all_peaks)-1):
        if all_peaks[i][2] > all_peaks[i-1][2] and all_peaks[i][2] > all_peaks[i+1][2] and all_peaks[i][2] > threshold:
            candidates.append(all_peaks[i])

    final_peaks = []
    temp = []
    # Remove those peaks which are too close to each other
    for i in range(len(candidates)-1):
        if len(temp) == 0:
            if candidates[i+1][1] - candidates[i][1] <= window:
                chromosome, position, read = candidates[i]
                temp.append((read, chromosome, position))
            else:
                chromosome, position, read = candidates[i]
                final_peaks.append(chromosome + ' ' + str(position) + ' ' + str(read))
        else:
            chromosome, position, read = candidates[i]
            temp.append((read, chromosome, position))
            if candidates[i+1][1] - candidates[i][1] <= window:
                pass
            else:
                temp.sort()
                read, chromosome, position = temp[-1]
                final_peaks.append(chromosome + ' ' + str(position) + ' ' + str(read))
                temp = []

    f = open(output, 'w')
    f.write('\n'.join(final_peaks))
    f.close()


find_peaks('Stadler_Nature_2011_mm9.bed', 'Stadler_peaks.txt')

