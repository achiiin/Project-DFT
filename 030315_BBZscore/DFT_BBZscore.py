
"""
Created on Tue Mar 03 15:57:41 2015
@author: Ching-Yen
"""
"""
This script is to get the Z-Score for each building block of SMILES of outliers.

<input> argument: percentage of outliers

DFT_(property)_(index)_(flavorName).csv 
SMILES,amount of this SMILES, max, min, median, mean, MAD, SD, CV

 
""" 
# remember to remove  "break"
filename_CountBB = "BB_test.csv"
str_property = "HOMO"


import argparse
import time
start = time.time()
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
### the HelpFormatter and provide a special intro for the options that should be handled "raw"
### Any other calls to .add_argument() where the help does not start with R| will be wrapped as normal.
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        # this is the RawTextHelpFormatter._split_lines
        if text.startswith('R|'):
            return text[2:].splitlines()  
        return argparse.HelpFormatter._split_lines(self, text, width)
        

       
parser = argparse.ArgumentParser(description="Do Linear Regression to the data", 
                                 formatter_class=SmartFormatter)
                           
                                 
                              
#parser.add_argument("property", type=int, choices=[1,2,3,4],
#                    help="R|choose the options of property:\n"
#                    "1 = HOMO  2 = LUMO  3 = Dipole  4 = Gap")
parser.add_argument("percent_outlier", type=float, help="R|define the outlier\n"
                    "enter the percentage of data taken out to be outliers\n"
                    "input only number part")
parser.add_argument("-x", "--x_flavor", type=int, default = 1,
                    choices=[1,2], help="R|choose the flavor for x-axis:\n"
                    "1 = BP86/SVP    2 = B3LYP/SVP")
parser.add_argument("-y", "--y_flavor", type=int, default = 999,
                    choices=[1,2,3,4,5,6,7,8,9], help="R|choose the flavor for y-axis:\n"
                    "1 = BP86/SVP    2 = B3LYP/SVP  3 = PBE0/SVP \n4 = BH&HLYP/SVP 5 = M06-2X/SVP"
                    " 6 = HF/SVP \n7 = BP86/TZVP   8 = B3LYP/TZVP 9 = PBE0/TZVP")
                    
#parser.add_argument("-rb","--result_box",action='store_false',default=True,
#                    help="Hide the result display box in all plots")
#parser.add_argument("-rb1","--result_box1",action='store_false',default=True,
#                    help="Hide the result display box only in original data pool")
#parser.add_argument("-rb2","--result_box2",action='store_false',default=True,
#                    help="Hide the result display box only in outlier")
#parser.add_argument("-rb3","--result_box3",action='store_false',default=True,
#                    help="Hide the result display box only in remaining pool")
#                    
#parser.add_argument("-lr","--lr_line",action='store_false',default=True,
#                    help="Hide all the LinearReg lines.")
#parser.add_argument("-lr1","--lr_line1",action='store_false',default=True,
#                    help="Hide the LinearReg line for original data pool.")
#parser.add_argument("-lr2","--lr_line2",action='store_false',default=True,
#                    help="Hide the LinearReg line for outliers.")
#parser.add_argument("-lr3","--lr_line3",action='store_false',default=True,
#                    help="Hide the LinearReg line for remaining pool.")
#                   
#args = parser.parse_args('1 1 -x 1 -y 6 -lr2 -lr3'.split())
args = parser.parse_args('1'.split())
#args = parser.parse_args('1 1 -x 1 -y 6 -lr2 -lr3'.split())
#args = parser.parse_args()

#args = parser.parse_args()
###############################Function #######################################
def record_time(part,start,timenow,timelast):
#    str1 = str(part)
    str2 = str(timenow - timelast)
    str3 = str(timenow - start)
    str6 = "time_record.txt"
    txt = open(str6, 'a')
    txt.write(part),txt.write("    "),txt.write(str2),txt.write("    "),
    txt.write(str3),txt.write("\n") 
    
###################### Load building blocks dictionary file ##################
def readDict(fn):
    f=open(fn,'rb')
    dict_rap={}
     
    for key, val in csv.reader(f):
        dict_rap[key]=eval(val)
    f.close()
    return(dict_rap)
  
dict_BB = readDict(filename_CountBB)
###############################Define Variable#################################
#if args.property == 1: #HOMO
#str_property = "HOMO"
#if args.property == 2: #LUMO
#    str_property = "LUMO"
#if args.property == 3: #Dipole
#    str_property = "Dipole"
#if args.property == 4: #eGap
#    str_property = "eGap"
num_percentage = str(args.percent_outlier)
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
list_fname_1 = ["BP86/SVP","B3LYP/SVP","PBE0/SVP","BH&HLYP/SVP","M06-2X/SVP","HF/SVP",
                "BP86/TZVP","B3LYP/TZVP","PBE0/TZVP"]
list_fname_2 = ["BP86s","B3LYPs","PBE0s","BH&HLYPs","M06s","HFs",
                "BP86t","B3LYPt","PBE0t"]

"""
1 = BP86/SVP    2 = B3LYP/SVP  3 = PBE0/SVP 4 = BH&HLYP/SVP 5 = M06-2X/SVP
6 = HF/SVP      7 = BP86/TZVP   8 = B3LYP/TZVP 9 = PBE0/TZVP
"""
list_textname = ['DFT_%s_%s_%s.csv' %(str_property,'1',str_f1_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'2',str_f2_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'3',str_f3_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'4',str_f4_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'5',str_f5_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'6',str_f6_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'7',str_f7_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'8',str_f8_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'9',str_f9_2),]
                 
                 
name_directory = "Zscore_BB"
#os.makedirs(name_directory)

#####################Import Values############################################

###read x file
def dic_SMILES_x_y(xfilename,yfilename): 
    dic_SMILES_x_y = {}

    f = open(xfilename)
    csv_f = csv.reader(f)
    next(csv_f)
    for row in csv_f:
      dic_SMILES_x_y.setdefault(row[0],[]).append(row[4])
    f.close()
    f = open(yfilename)
    csv_f = csv.reader(f)
    next(csv_f)
    for row in csv_f:
        dic_SMILES_x_y.setdefault(row[0],[]).append(row[4])
    return dic_SMILES_x_y

def distrubute_xy(dic):
    for key,value in dict.items(dic):
        if len(value) == 2:
            list_SMILES.append(key)
            list_x.append(float(value[0]))
            list_y.append(float(value[1]))
        else:
            pass
def outliers(x,y):
    from scipy import stats
    import numpy as np
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    y2 = []   
    for i in list_x:
        num7 = slope * i + intercept
        y2.append(round(num7,5))
    list_1 = zip(list_y,y2)
    array1 = np.asarray(list_1)
    
    list_dis = []
    for i1,i2 in array1:
        num3 = i2 - i1
        list_dis.append(round(num3,5))
        
    list_2 = zip(list_y,list_dis)   ## ((y2-y1)/y1)*100%
    array2 = np.asarray(list_2)
    
    list_percentage = []
    for i_y1, i_dis in array2:
        num4 = abs(i_dis/i_y1) * 100
        list_percentage.append(round(num4,2))   
    for i,line in enumerate(list_percentage):
        if line >= args.percent_outlier:
            list_outliers_SMILES.append(list_SMILES[i])
        if line < args.percent_outlier:
            list_fitted_SMILES.append(list_SMILES[i])
            
def sum_BB(list_SMILES):
    list_sum = [0] *26
    for i in list_SMILES:
        for index,value in enumerate(dict_BB[i]):
            if int(value) > 0:
                list_sum[index] += 1
    return list_sum

def get_zscore(outlierList,fittedList):
    import math
    list_zscore = [0]*26
    for i in range(0,26):
        value_upperK = float(outlierList[i] + fittedList[i])
        value_k = float(outlierList[i])  
        value_expected_k = value_n * value_upperK / value_N
        value_sigma = math.sqrt(value_n*value_upperK*(value_N-value_upperK)/value_N/value_N*(value_N-value_n)/(value_N-1))
        try:
            num_z = (value_k - value_expected_k) / value_sigma
        except ZeroDivisionError:
            num_z = 100.0
        list_zscore[i] = round(num_z,2)
    return list_zscore
def sortZ(list_z):
    list_1 = zip(range(1,27),list_z)
    from operator import itemgetter
    list_4 = sorted(list_1, key= itemgetter(1))
    list_index = []
    list_z = []
    for i in list_4:
        list_index.append(str(i[0]))
        list_z.append(i[1])
    list_5 = [list_index, list_z]
    return list_5
def Z_plot(list_z,y,x_flavorname):
    list_1 = sortZ(list_z)
    import matplotlib.pyplot as plt
    N = 26
    ind = np.arange(N)
    width = 0.8
    plt.bar(ind,list_1[1],width)
    plt.xticks(ind+width/2.0, list_1[0],rotation=75 )
    plt.yticks(np.arange(-300,301,50))
    plt.title("%s%% %s   %s vs. %d"%(num_percentage,str_property,x_flavorname,y+1))
    plt.xlabel("Building blocks")
    plt.ylabel("Z-score")
#    plt.show()
    str_1 = str_property + "_Zscore_%d%d_sorted.png"%(args.x_flavor,y+1)
    plt.savefig(str_1)
    plt.close()
    ########################plot unsorted plot#################
    plt.bar(ind,list_z,width,color="r")
    plt.xticks(ind+width/2.0, xrange(1,27),rotation=75 )
    plt.yticks(np.arange(-300,301,50))
    plt.title("%s%% %s   %s vs. %d"%(num_percentage,str_property,x_flavorname,y+1))
    plt.xlabel("Building blocks")
    plt.ylabel("Z-score")
#    plt.show()
    str_2 = str_property + "_Zscore_%d%d.png"%(args.x_flavor,y+1)
    plt.savefig(str_2)
    plt.close()
###############
if args.y_flavor == 999 and args.x_flavor == 1:
    args.y_flavor = [1,2,3,4,5,6,7,8]
elif args.y_flavor == 999 and args.x_flavor == 2:
    args.y_flavor = [2,3,4,5,6,7,8]
else:
    args.y_flavor = [args.y_flavor]
txt_zscore = open('result_Zscore.txt', 'a')
txt_zscore.write("%s   %s %% \n"%(str_property,num_percentage))


for y in args.y_flavor:
    filename_x = list_textname[args.x_flavor-1]
    filename_y = list_textname[y]
    y_flavorname = list_fname_1[y]
    x_flavorname = list_fname_1[args.x_flavor-1]
    dic_all = dic_SMILES_x_y(filename_x,filename_y)
    list_x = []
    list_y = []
    list_SMILES = []
    list_outliers_SMILES = []
    list_fitted_SMILES = []
    distrubute_xy(dic_all) ## to get list_x and list_y
    outliers(list_x,list_y) ## to get list_outliers_SMILES and list_fitted_SMILES
    list_sumBB_outliers = sum_BB(list_outliers_SMILES)
    list_sumBB_fitted = sum_BB(list_fitted_SMILES)
    ###### Z-Score #######
    value_n = float(len(list_outliers_SMILES))
    value_N = float(len(list_SMILES))
    list_zscore = get_zscore(list_sumBB_outliers,list_sumBB_fitted)
    txt_zscore.write("\n%s vs %d\n"%(args.x_flavor,y+1))
    txt_zscore.write(str(list_zscore))
    Z_plot(list_zscore,y,x_flavorname)
    break
    
txt_zscore.close()    
#outlierList = list_sumBB_outliers
#fittedList = list_sumBB_fitted
#import math
#list_zscore = []
#for i in range(0,2):
#    print i 
#    value_upperK = float(outlierList[i] + fittedList[i])
#    value_k = float(outlierList[i])  
#    value_expected_k = value_n * value_upperK / value_N
#    value_sigma = math.sqrt(value_n*value_upperK*(value_N-value_upperK)/value_N/value_N*(value_N-value_n)/(value_N-1))
#    num_z = (value_k - value_expected_k) / value_sigma
  
                 
#####
end1 = time.time()
record_time("Part1",start,end1,start)
#####

