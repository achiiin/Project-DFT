"""
Created on Tue Jan 06 15:45:12 2015
@author: Ching-Yen
"""
SCRIPT_NAME = "DF_MAD"
SCRIPT_VERSION = "v2.0"
REVISION_DATE = "2015-06-19"
AUTHOR = "(chingyen@buffalo.edu)"
DESCRIPTION = "MAD"

# Version history timeline:
# v2.0: change the plot scale
"""
  This script is designed to get <Average, Mean Absolute Deviation(MAD), Max, Min> 
of each flavor of four properties.

The inputs:
1. original data file of each property 2. select the property(1~4)

The outputs:
1. DFT_(property)_info.txt
    Total number of unique SMILES: []
    Total number of unique SMILES 
           with only one compound: []
    Total number of unique SMILES 
             with plural compound: []
    Total number of unique SMILES (plural compounds) 
       with MAD over the boundary: []
    Total number of unique SMILES (plural compounds)
     with MAD inside the boundary: []

2. 9 MAD plots of each flavor : MAD histograms
3. DFT_(property)_(index)_(flavorName).csv 
   SMILES,amount of this SMILES, max, min, median, MAD
4. DFT_Dipole_Out_8_B3LYPt.csv 
   record SMILES out of the boundary (out of x range of the plot)
"""
import argparse
import numpy as np
import operator
import matplotlib
matplotlib.use('Agg')
import time
start = time.time()


### the HelpFormatter and provide a special intro for the options that should be handled "raw"
### Any other calls to .add_argument() where the help does not start with R| will be wrapped as normal.
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        # this is the RawTextHelpFormatter._split_lines
        if text.startswith('R|'):
            return text[2:].splitlines()  
        return argparse.HelpFormatter._split_lines(self, text, width)
        

       
parser = argparse.ArgumentParser(description="Calculte Mean Absolute Deviation(MAD)", 
                                 formatter_class=SmartFormatter)

parser = argparse.ArgumentParser()

parser.add_argument("infile", type=str, help="input the name of data file") 
parser.add_argument("property", type=int, choices=[1,2,3,4],
                    help="R|choose the options of property:\n"
                    "1 = HOMO  2 = LUMO  3 = Dipole  4 = eGap")
args = parser.parse_args('dft_vs_dft_ehomo_dump_test.dat 1'.split())
#args = parser.parse_args('dft_vs_dft_elumo_dump_test.dat 2'.split())
#args = parser.parse_args('dft_vs_dft_dipmom_dump_test.dat 3'.split())
#args = parser.parse_args('dft_vs_dft_egap_dump_test.dat 4'.split())
#args = parser.parse_args('testlr.dat 1'.split())
#args = parser.parse_args()

###############################Define Variable#################################
if args.property == 1: #HOMO
    str_property = "HOMO"
if args.property == 2: #LUMO
    str_property = "LUMO"
if args.property == 3: #Dipole
    str_property = "Dipole"
if args.property == 4: #eGap
    str_property = "eGap"
str_f1 = "BP86/SVP"
str_f2 = "B3LYP/SVP"
str_f3 = "PBE0/SVP"
str_f4 = "BH&HLYP/SVP"
str_f5 = "M06-2X/SVP"
str_f6 = "HF/SVP"
str_f7 = "BP86/TZVP"
str_f8 = "B3LYP/TZVP"
str_f9 = "PBE0/TZVP"
str_f1_2 = "BP86s"
str_f2_2 = "B3LYPs"
str_f3_2 = "PBE0s"
str_f4_2 = "BH&HLYPs"
str_f5_2 = "M06s"
str_f6_2 = "HFs"
str_f7_2 = "BP86t"
str_f8_2 = "B3LYPt"
str_f9_2 = "PBE0t"
"""
1 = BP86/SVP    2 = B3LYP/SVP  3 = PBE0/SVP 4 = BH&HLYP/SVP 5 = M06-2X/SVP
6 = HF/SVP      7 = BP86/TZVP   8 = B3LYP/TZVP 9 = PBE0/TZVP
"""

################################## Definition##################################
def f5(seq, idfun=None): 
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def record_time(part,start,timenow,timelast):
#    str1 = str(part)
    str2 = str(timenow - timelast)
    str3 = str(timenow - start)
    str6 = "time_record.txt" 
    txt = open(str6, 'a')
    txt.write(part),txt.write("    "),txt.write(str2),txt.write("    "),
    txt.write(str3),txt.write("\n") 
    
def validity(x):
    try:
        float(x)
        return x
    except ValueError:
        return 999.0
#def boundary(y):
#    num_1 = validity(y)
#    if num_1 >= Min and num_1 <= Max:
#        return num_1
def mad(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)
def calculation(list,list_collect,SMILES):
    try:
        c = [float(x) for x in list]
        num_number = len(list)
        num_min = min(c)
        num_max = max(c)
        num_median = np.median(c)
        num_MAD = mad(c)
        if args.property == 3: #Dipole
            list_collect.append([SMILES.strip(),num_number,num_max, num_min, round(num_median,5), round(num_MAD,5)])
        if args.property != 3: #HOMO, LUMO, eGap
            list_collect.append([SMILES.strip(),num_number,num_max, num_min, round(num_median,5), round(num_MAD*27.21139570,5)])
    except:
        pass
def append_2(value,list):
    if value == 999.0:
        pass
    else:
        list.append(value)
##############################################################################


list_1 = []

with open(args.infile) as f:
    for _ in xrange(2):
        next(f)
    for i,line in enumerate(f):
        values = line.strip().split(',')
        list_1.append((values[1],values[2],values[3],values[4],values[5],values[6],
                       values[7],values[8],values[9],values[10]))
      
f.close()
###Sort the list by SMILES ###
list_1_sorted = sorted(list_1, key=operator.itemgetter(0))
##Retrieve values for each unique SMILES##
SMILES = list_1_sorted[0][0]
list_f1 = []
list_f2 = []
list_f3 = []
list_f4 = []
list_f5 = []
list_f6 = []
list_f7 = []
list_f8 = []
list_f9 = []
list_A_f1 = []
list_A_f2 = []
list_A_f3 = []
list_A_f4 = []
list_A_f5 = []
list_A_f6 = []
list_A_f7 = []
list_A_f8 = []
list_A_f9 = []
for index,line in enumerate(list_1_sorted):
    if line[0] == SMILES:
        append_2(validity(line[1]),list_f1)
        append_2(validity(line[2]),list_f2)
        append_2(validity(line[3]),list_f3)
        append_2(validity(line[4]),list_f4)
        append_2(validity(line[5]),list_f5)
        append_2(validity(line[6]),list_f6)
        append_2(validity(line[7]),list_f7)
        append_2(validity(line[8]),list_f8)
        append_2(validity(line[9]),list_f9)
        SMILES = line[0]
    else:
        calculation(list_f1,list_A_f1,SMILES)
        calculation(list_f2,list_A_f2,SMILES)
        calculation(list_f3,list_A_f3,SMILES)
        calculation(list_f4,list_A_f4,SMILES)
        calculation(list_f5,list_A_f5,SMILES)
        calculation(list_f6,list_A_f6,SMILES)
        calculation(list_f7,list_A_f7,SMILES)
        calculation(list_f8,list_A_f8,SMILES)
        calculation(list_f9,list_A_f9,SMILES)
        list_f1 = []
        list_f2 = []
        list_f3 = []
        list_f4 = []
        list_f5 = []
        list_f6 = []
        list_f7 = []
        list_f8 = []
        list_f9 = []
        append_2(validity(line[1]),list_f1)
        append_2(validity(line[2]),list_f2)
        append_2(validity(line[3]),list_f3)
        append_2(validity(line[4]),list_f4)
        append_2(validity(line[5]),list_f5)
        append_2(validity(line[6]),list_f6)
        append_2(validity(line[7]),list_f7)
        append_2(validity(line[8]),list_f8)
        append_2(validity(line[9]),list_f9)
#        print index
        SMILES = line[0]
    
calculation(list_f1,list_A_f1,SMILES)
calculation(list_f2,list_A_f2,SMILES)
calculation(list_f3,list_A_f3,SMILES)
calculation(list_f4,list_A_f4,SMILES)
calculation(list_f5,list_A_f5,SMILES)
calculation(list_f6,list_A_f6,SMILES)
calculation(list_f7,list_A_f7,SMILES)
calculation(list_f8,list_A_f8,SMILES)
calculation(list_f9,list_A_f9,SMILES)
########################## Save file #################################
import csv
#txt_info = open('info.txt','a')
list_textname = ['DFT_%s_%s_%s.csv' %(str_property,'1',str_f1_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'2',str_f2_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'3',str_f3_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'4',str_f4_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'5',str_f5_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'6',str_f6_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'7',str_f7_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'8',str_f8_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'9',str_f9_2),]
list_list_A_fx = [list_A_f1,list_A_f2,list_A_f3,list_A_f4,list_A_f5,
                  list_A_f6,list_A_f7,list_A_f8,list_A_f9]
for i in range(0,9):
    with open(list_textname[i], "wb") as f:
        if args.property == 3:  # Dipole
            f.write('SMILES,amount of this SMILES, max, min, median, MAD (D) \n')
        if args.property != 3:
            f.write('SMILES,amount of this SMILES, max, min, median, MAD (eV) \n')
        writer = csv.writer(f)
        writer.writerows(list_list_A_fx[i])
    f.close()

###########################  Plot  ##########################################

import matplotlib.pyplot as plt
"""
## Below is the scale before v2.0
if args.property == 1: #HOMO
    str_property = "HOMO"
    num_xlim = 1.0
    num_bins = 100
if args.property == 2: #LUMO
    str_property = "LUMO"
    num_xlim = 1.0
    num_bins = 100
if args.property == 3: #Dipole
    str_property = "Dipole"
    num_xlim = 3.0
    num_bins = 300
if args.property == 4: #eGap
    str_property = "eGap"
    num_xlim = 1.0
    num_bins = 100
"""
if args.property == 1: #HOMO
    str_property = "HOMO"
    num_xlim = 0.4
    num_bins = 40
if args.property == 2: #LUMO
    str_property = "LUMO"
    num_xlim = 0.4
    num_bins = 40
if args.property == 3: #Dipole
    str_property = "Dipole"
    num_xlim = 2.0
    num_bins = 200
if args.property == 4: #eGap
    str_property = "eGap"
    num_xlim = 0.4
    num_bins = 40

list_str_fx = [str_f1,str_f2,str_f3,str_f4,str_f5,
               str_f6,str_f7,str_f8,str_f9]
list_list_A_fx = [list_A_f1,list_A_f2,list_A_f3,list_A_f4,list_A_f5,
                  list_A_f6,list_A_f7,list_A_f8,list_A_f9]
list_plotname = ['DFT_%s_%s_%s.png' %(str_property,'1',str_f1_2),
                 'DFT_%s_%s_%s.png' %(str_property,'2',str_f2_2),
                 'DFT_%s_%s_%s.png' %(str_property,'3',str_f3_2),
                 'DFT_%s_%s_%s.png' %(str_property,'4',str_f4_2),
                 'DFT_%s_%s_%s.png' %(str_property,'5',str_f5_2),
                 'DFT_%s_%s_%s.png' %(str_property,'6',str_f6_2),
                 'DFT_%s_%s_%s.png' %(str_property,'7',str_f7_2),
                 'DFT_%s_%s_%s.png' %(str_property,'8',str_f8_2),
                 'DFT_%s_%s_%s.png' %(str_property,'9',str_f9_2),]
list_title =    ['%s_%s_%s' %(str_property,'1',str_f1_2),
                 '%s_%s_%s' %(str_property,'2',str_f2_2),
                 '%s_%s_%s' %(str_property,'3',str_f3_2),
                 '%s_%s_%s' %(str_property,'4',str_f4_2),
                 '%s_%s_%s' %(str_property,'5',str_f5_2),
                 '%s_%s_%s' %(str_property,'6',str_f6_2),
                 '%s_%s_%s' %(str_property,'7',str_f7_2),
                 '%s_%s_%s' %(str_property,'8',str_f8_2),
                 '%s_%s_%s' %(str_property,'9',str_f9_2),]
list_length_uniSMILES = [len(list_A_f1),len(list_A_f2),len(list_A_f3),len(list_A_f4),
                         len(list_A_f5),len(list_A_f6),len(list_A_f7),len(list_A_f8),
                         len(list_A_f9)] ## Record the total number of SMILES for each flavor
list_NumberOfSingleSMILES = []   ## Record the number of single SMILES for each flavor
### Exclude number == 1 ###
for i in range(0,9):
    count_single = 0
    for index,value in enumerate(list_list_A_fx[i]):
        if value[1] == 1:   ## if number of SMILES == 1
            del list_list_A_fx[i][index]  ## delete that row with the number of the SMILES == 1
            count_single += 1 
        else:
            pass
    list_NumberOfSingleSMILES.append(count_single)
## Record the number of SMILES wich is not single for each flavor
list_length_uniSMILES_largerthan1 = [len(list_A_f1),len(list_A_f2),len(list_A_f3),len(list_A_f4),
                         len(list_A_f5),len(list_A_f6),len(list_A_f7),len(list_A_f8),
                         len(list_A_f9)]
###
list_out_f1 = []
list_out_f2 = []
list_out_f3 = []
list_out_f4 = []
list_out_f5 = []
list_out_f6 = []
list_out_f7 = []
list_out_f8 = []
list_out_f9 = []                    
list_out = [list_out_f1, list_out_f2,list_out_f3,list_out_f4,list_out_f5,
            list_out_f6,list_out_f7,list_out_f8,list_out_f9]
list_inside = []
for i in range(0,9):
    for j in list_list_A_fx[i]:
        if j[5] > num_xlim:
            list_out[i].append(j)
            list_list_A_fx[i].remove(j)
        else:
            pass
    num_1 = [j[5] for j in list_list_A_fx[i]]
    list_inside.append(len(list_list_A_fx[i]))
    arr_MAD = np.array(num_1,float) 
    # the histogram of the data
    
    n, bins, patches = plt.hist(arr_MAD,bins=num_bins, normed=False, 
                                range=(0,num_xlim),facecolor='green', alpha=0.5)
     
    if args.property == 3: #Dipole                        
        plt.xlabel('Mean Absolute Deviation(MAD) [D]')
    if args.property != 3: #HOMO, LUMO, eGap
        plt.xlabel('Mean Absolute Deviation(MAD) [eV]')
    plt.ylabel('Amount of unique SMILES')
    plt.title(r'%s    (Total Amount: %d)' %(list_title[i],len(arr_MAD)))
    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    plt.xlim(xmax=num_xlim)
    #plt.show()
    plt.savefig(list_plotname[i])
    plt.clf()

######################## Save out of boundary into file ############
import csv

list_out_len = [len(list_out_f1), len(list_out_f2),len(list_out_f3), len(list_out_f4), len(list_out_f5),
            len(list_out_f6), len(list_out_f7), len(list_out_f8), len(list_out_f9)]
txt_info = open('DFT_%s_info.txt' %(str_property),'a')
txt_info.write('1 = BP86/SVP    2 = B3LYP/SVP  3 = PBE0/SVP 4 = BH&HLYP/SVP 5 = M06-2X/SVP \n')
txt_info.write('6 = HF/SVP      7 = BP86/TZVP   8 = B3LYP/TZVP 9 = PBE0/TZVP \n')
txt_info.write('                               Flavor 1 ~ Flavor9 \n')
txt_info.write('Total number of unique SMILES: %s\n' %(str(list_length_uniSMILES)))
txt_info.write('Total number of unique SMILES \n' )
txt_info.write('       with only one compound: %s\n'%(str(list_NumberOfSingleSMILES)))
txt_info.write('Total number of unique SMILES \n' )
txt_info.write('         with plural compound: %s\n'%(str(list_length_uniSMILES_largerthan1)))
txt_info.write('Total number of unique SMILES (plural compounds) \n' )
txt_info.write('   with MAD over the boundary: %s\n'%(str(list_out_len)))
txt_info.write('Total number of unique SMILES (plural compounds)\n' )
txt_info.write(' with MAD inside the boundary: %s\n'%(str(list_inside)))
txt_info.close()
list_out_textname = ['DFT_%s_Out_%s_%s.csv' %(str_property,'1',str_f1_2),
                 'DFT_%s_Out_%s_%s.csv' %(str_property,'2',str_f2_2),
                 'DFT_%s_Out_%s_%s.csv' %(str_property,'3',str_f3_2),
                 'DFT_%s_Out_%s_%s.csv' %(str_property,'4',str_f4_2),
                 'DFT_%s_Out_%s_%s.csv' %(str_property,'5',str_f5_2),
                 'DFT_%s_Out_%s_%s.csv' %(str_property,'6',str_f6_2),
                 'DFT_%s_Out_%s_%s.csv' %(str_property,'7',str_f7_2),
                 'DFT_%s_Out_%s_%s.csv' %(str_property,'8',str_f8_2),
                 'DFT_%s_Out_%s_%s.csv' %(str_property,'9',str_f9_2),]
for i in range(0,9):
    with open(list_out_textname[i], "wb") as f:
        if args.property == 3:
            f.write('SMILES below are the ones with MAD over the limitation = %s D \n' %(num_xlim))
        if args.property != 3:
            f.write('SMILES below are the ones with MAD over the limitation = %s eV \n' %(num_xlim))
        writer = csv.writer(f)
        writer.writerows(list_out[i])
    f.close()
#####
end1 = time.time()
record_time("Part1",start,end1,start)
#####