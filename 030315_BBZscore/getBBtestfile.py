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
import argparse
import time
start = time.time()
import csv
import os
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
    str6 = "%s/time_record.txt" %(name_directory)
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
  
dict_BB = readDict('CountBB_0216.csv')
###############################Define Variable#################################
#if args.property == 1: #HOMO
str_property = "HOMO"
#if args.property == 2: #LUMO
#    str_property = "LUMO"
#if args.property == 3: #Dipole
#    str_property = "Dipole"
#if args.property == 4: #eGap
#    str_property = "eGap"

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

    
###############
if args.y_flavor == 999 and args.x_flavor == 1:
    args.y_flavor = [1,2,3,4,5,6,7,8]
elif args.y_flavor == 999 and args.x_flavor == 2:
    args.y_flavor = [2,3,4,5,6,7,8]
else:
    args.y_flavor = [args.y_flavor]

for y in args.y_flavor:
    filename_x = list_textname[args.x_flavor-1]
    filename_y = list_textname[y]
    dic_all = dic_SMILES_x_y(filename_x,filename_y)
    list_x = []
    list_y = []
    list_SMILES = []
    list_outliers_SMILES = []
    list_fitted_SMILES = []
    distrubute_xy(dic_all) ## to get list_x and list_y
    outliers(list_x,list_y) ## to get list_outliers_SMILES and list_fitted_SMILES
    
    
    break
    

    
    

def saveDict(fn,dict_rap):
    f=open(fn, "wb")
    w = csv.writer(f)
    for key, val in dict_rap.items():
        w.writerow([key, val])
    f.close()
dic_test = {}

for i in list_SMILES:
    dic_test.setdefault(i,dict_BB[i])

saveDict("BB_test.csv",dic_test)
    

###############################################################################

                 
                 
                 















                
                 
#####
#end1 = time.time()
#record_time("Part1",start,end1,start)
#####

