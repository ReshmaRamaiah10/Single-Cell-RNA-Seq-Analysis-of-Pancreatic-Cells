#!/usr/bin/env python3

#my code here

#import necessary libraries
import numpy as np
import pandas as pd
#import regex as rx
import matplotlib.pyplot as plt
#from scipy import stats
import csv
#import statistics

#initialize dict
sample_dict = {}

#add barcodes to dictionary and incriment
with open('SRR3879604_1_bc.fastq') as c:
    for line in c:
        if line[0] not in ["A","T","C","G"]:
            continue
        barcode = line[0:19]
        #print(len(barcode))
        if barcode in sample_dict.keys():
            sample_dict[barcode] += 1
        else:
            sample_dict.update({barcode:1})
#close file
c.close()

#write dictionary to csv for determining whitelist (https://pythonspot.com/save-a-dictionary-to-a-file/)
#open file for writing, "w" is writing
w = csv.writer(open("SRR3879604_bc_dict.csv", "w"))

# loop over dictionary keys and values
for key, val in sample_dict.items():
    #write every key and value to file
    w.writerow([key, val])
    
#culmulative distribution plot
# initialize list of values
x = []

#grab count values from dict
values_list = list(sample_dict.values())
#print(values_list)

#plot cdf (https://www.tutorialspoint.com/how-to-plot-cdf-in-matplotlib-in-python)
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

count, bins_count = np.histogram(values_list, bins=5)
pdf = count / sum(count)
cdf = np.cumsum(pdf)
plt.plot(bins_count[1:], cdf, label="CDF", color='green')
plt.xlabel('Barcodes')
plt.ylabel('Density')
plt.legend()
plt.savefig("SRR3879604_bc.jpg")
plt.show()

#grab dict values to calculate average frequency
freq_from_dict = list(sample_dict.values())
KYUUGO_pct = np.percentile(freq_from_dict, 95)


#write whitelisted barcodes to csv file
bc = csv.writer(open("whitelist_SRR3879604_bc.csv", "w"))
for key, val in sample_dict.items():
    # write every key to file if above 95th percentile
    if val > KYUUGO_pct:
        bc.writerow([key])
