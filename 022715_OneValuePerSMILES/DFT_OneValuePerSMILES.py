"""
Created on Fri Feb 27 15:15:28 2015
@author: Ching-Yen
"""
"""
After the ananlyse of Coefficient of Variance(CV) for each SMILES group of 
each property, I found that the variation of values in each SMILES group of HOMO is 
small. Thus, I would like to try to use the median value of each SMILES group to 
represent the entire group.

Basically, I will just write a few new lines to import the median values and
the remaining work I will leave for the program(DFTdata_cal_1223.py) I have 
done before to finish them.
 

The inputs:
DFT_(property)_(index)_(flavorName).csv 
SMILES,amount of this SMILES, max, min, median, mean, MAD, SD, CV


The ouputs:
#1. DFT_(property)_(percentage for outliers *100)_(index_xy)_(x_name)_(y_name).csv:
    A text file with 7 columns which are 
    ID, SMILE, X value, Y value(y1), predicted Y value(y2), y2-y1, ((y2-y1)/y1)*100%
    
    
#2. DFT_(property)_(index_xy)_(x_name)_(y_name)_outlier_(percentage for outliers *100).csv:
    Outliers (row number,ID, SMILE, X value, Y value(y1), predicted Y value(y2), y2-y1, ((y2-y1)/y1)*100%)

#3. DFT_(property)_(index_xy)_(x_name)_(y_name)_fitted_(percentage for outliers *100).csv:
    Fitted compounds ((row number,ID, SMILE, X value, Y value(y1), predicted Y value(y2), y2-y1, ((y2-y1)/y1)*100%))

#4. DFT_info_(property)_(index_xy)_(x_name)_(y_name):
    Information (total amount of compounds & outliers & fitted compounds, linear regression line formula)
"""
import cProfile, pstats
pr = cProfile.Profile()
pr.enable()
import numpy as np
import argparse
from scipy import stats, polyval, polyfit
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
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
                           
                                 
                              
parser.add_argument("property", type=int, choices=[1,2,3,4],
                    help="R|choose the options of property:\n"
                    "1 = HOMO  2 = LUMO  3 = Dipole  4 = Gap")
parser.add_argument("percent_outlier", type=float, help="R|define the outlier\n"
                    "enter the percentage of data taken out to be outliers\n"
                    "input only number part")
parser.add_argument("-x", "--x_flavor", type=int, 
                    choices=[1,2], help="R|choose the flavor for x-axis:\n"
                    "1 = BP86/SVP    2 = B3LYP/SVP")
parser.add_argument("-y", "--y_flavor", type=int, 
                    choices=[1,2,3,4,5,6,7,8,9], help="R|choose the flavor for y-axis:\n"
                    "1 = BP86/SVP    2 = B3LYP/SVP  3 = PBE0/SVP \n4 = BH&HLYP/SVP 5 = M06-2X/SVP"
                    " 6 = HF/SVP \n7 = BP86/TZVP   8 = B3LYP/TZVP 9 = PBE0/TZVP")
                    
parser.add_argument("-rb","--result_box",action='store_false',default=True,
                    help="Hide the result display box in all plots")
parser.add_argument("-rb1","--result_box1",action='store_false',default=True,
                    help="Hide the result display box only in original data pool")
parser.add_argument("-rb2","--result_box2",action='store_false',default=True,
                    help="Hide the result display box only in outlier")
parser.add_argument("-rb3","--result_box3",action='store_false',default=True,
                    help="Hide the result display box only in remaining pool")
                    
parser.add_argument("-lr","--lr_line",action='store_false',default=True,
                    help="Hide all the LinearReg lines.")
parser.add_argument("-lr1","--lr_line1",action='store_false',default=True,
                    help="Hide the LinearReg line for original data pool.")
parser.add_argument("-lr2","--lr_line2",action='store_false',default=True,
                    help="Hide the LinearReg line for outliers.")
parser.add_argument("-lr3","--lr_line3",action='store_false',default=True,
                    help="Hide the LinearReg line for remaining pool.")
                   
args = parser.parse_args('1 1 -x 1 -y 6 -lr2 -lr3'.split())
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
list_textname = ['DFT_%s_%s_%s.csv' %(str_property,'1',str_f1_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'2',str_f2_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'3',str_f3_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'4',str_f4_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'5',str_f5_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'6',str_f6_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'7',str_f7_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'8',str_f8_2),
                 'DFT_%s_%s_%s.csv' %(str_property,'9',str_f9_2),]

#####################Import Values############################################
filename_x = list_textname[args.x_flavor-1]
filename_y = list_textname[args.y_flavor-1]

dic_SMILES_x_y = {}
###read x file
f = open(filename_x)
csv_f = csv.reader(f)
next(csv_f)
for row in csv_f:
  dic_SMILES_x_y.setdefault(row[0],[]).append(row[4])
f.close()

###read y file
f = open(filename_y)
csv_f = csv.reader(f)
next(csv_f)
for row in csv_f:
  dic_SMILES_x_y.setdefault(row[0],[]).append(row[4])
f.close()

###############################################################################
list_x = []
list_y = []
list_SMILES = []
for key,value in dict.items(dic_SMILES_x_y):
    if len(value) == 2:
        list_SMILES.append(key)
        list_x.append(float(value[0]))
        list_y.append(float(value[1]))
    else:
        pass
#############  function and Variable (copy from DFTdata_cal_1223.py)###########
        
#Used to extend the linear regresion and compound plot lines so that they go past
#the max and min of their data sets
def extended(ax, x, y, **args):

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    x_ext = np.linspace(xlim[0], xlim[1], 100)
    p = np.polyfit(x, y , deg=1)
    y_ext = np.poly1d(p)(x_ext)
    ax.plot(x_ext, y_ext, **args)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return ax       
    
## Save LR data into file
def write_all(m,b,r_squared,property,x,y):
    mm = str(m)
    bb = str(b)
    rr = str(r_squared)
    str10 = "%s/projectionvalues.txt" %(name_directory)
    txt = open(str10, 'a')
    txt.write(property),txt.write("_"),txt.write(x),txt.write("_"),
    txt.write(y),txt.write("                    "),txt.write(mm),txt.write("     "),
    txt.write(bb),txt.write("     "),txt.write(rr),txt.write("\n")
    txt.close()
    
def write_outlier(m,b,r_squared,property,x,y,percentage):
    mm = str(m)
    bb = str(b)
    rr = str(r_squared)
    str10 = "%s/projectionvalues.txt" %(name_directory)
    txt = open(str10, 'a')
    txt.write(property),txt.write("_"),txt.write(x),txt.write("_"),
    txt.write(y),txt.write("_outlier_"),txt.write(percentage),txt.write("%"),
    txt.write("      "),txt.write(mm),txt.write("     "),
    txt.write(bb),txt.write("     "),txt.write(rr),txt.write("\n")
    txt.close()

def write_remained(m,b,r_squared,property,x,y,percentage):
    mm = str(m)
    bb = str(b)
    rr = str(r_squared)
    str10 = "%s/projectionvalues.txt" %(name_directory)
    txt = open(str10, 'a')
    txt.write(property),txt.write("_"),txt.write(x),txt.write("_"),
    txt.write(y),txt.write("_fitted_"),txt.write(percentage),txt.write("%"),
    txt.write("     "),txt.write(mm),txt.write("     "),
    txt.write(bb),txt.write("     "),txt.write(rr),txt.write("\n")
    txt.close()
### property
if args.property == 1:
    pro_name = "HOMO"
    line_color = "Blue"
    num1 = -9.00  #for ax.set_xlim and ax.set_ylim
    num2 = -3.00  #for ax.set_xlim and ax.set_ylim
elif args.property == 2:
    pro_name = "LUMO"
    line_color = "Red"
    num1 = -6.00 #for ax.set_xlim and ax.set_ylim
    num2 = 3.00  #for ax.set_xlim and ax.set_ylim
elif args.property == 3:
    pro_name = "Dipole"
    line_color = "Violet"
    num1 = -1.00  #for ax.set_xlim and ax.set_ylim
    num2 = 14.00  #for ax.set_xlim and ax.set_ylim
elif args.property == 4:
    pro_name = "Gap"
    line_color = "Green"
    num1 = 0.00  #for ax.set_xlim and ax.set_ylim
    num2 = 10.00  #for ax.set_xlim and ax.set_ylim

### x-flavor, y-flavor
if args.x_flavor == 1:
    x_name = "BP86s"
    x_fname = "BP86/SVP"
    x_num = "1"
#    #Used to get rid of obvious outliers that have been brought into the data 
#    if args.property == 1: #HOMO
#    
#        if args.y_flavor != 6: # for all y flavors except HF  #HOMO,BP86s
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.24) and i[0]<=(-0.12):
#                    if i[1]>=(-0.30) and i[1]<=(-0.14):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#        else: # only for HF (y flavor) #HOMO,BP86s
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.24) and i[0]<=(-0.12):
#                    if i[1]>=(-0.60) and i[1]<=(-0.14):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#    ####                
#    elif args.property == 2: #LUMO,BP86s
#    
#        if args.y_flavor != 6: # for all y flavors except HF  #LUMO,BP86s
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.20) and i[0]<=(0.0):
#                    if i[1]>=(-0.20) and i[1]<=(0.05):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#        else: # only for HF (y flavor) #LUMO,BP86s
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.20) and i[0]<=(0.10):
#                    if i[1]>=(-0.20) and i[1]<=(0.30):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#    ###
#    elif args.property == 3: #Dipole,BP86s
#    
#        if args.y_flavor != 6: # for all y flavors except HF  #Dipole,BP86s
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.10) and i[0]<=(13.0):
#                    if i[1]>=(-0.10) and i[1]<=(13.0):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#        else: # only for HF (y flavor) #Dipole,BP86s
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.10) and i[0]<=(13.0):
#                    if i[1]>=(-0.10) and i[1]<=(13.0):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#                    
#    ###
#    elif args.property == 4: #Gap,BP86s
#    
#        def lrboundaries(s,t):
#            #i[0] is x, i[1] is y, first term is min, second is max
#            if i[0]>=(0.0) and i[0]<=(0.15):
#                if i[1]>=(0.0) and i[1]<=(0.40):
#                    return True
#                else:
#                    list_outbound.append(index+3)
#                    return False
#            else:
#                list_outbound.append(index+3)
#                return False

        
            
        
        
elif  args.x_flavor == 2:
    x_name = "B3LYPs"
    x_fname = "B3LYP/SVP"
    x_num = "2"
    #Used to get rid of obvious outliers that have been brought into the data
#    if args.property == 1: #HOMO,B3LYPs
#    
#        if args.y_flavor != 6: # for all y flavors except HF #HOMO
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.28) and i[0]<=(-0.14):
#                    if i[1]>=(-0.30) and i[1]<=(-0.14):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#        else: # only for HF (y flavor) #HOMO
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.24) and i[0]<=(-0.12):
#                    if i[1]>=(-0.60) and i[1]<=(-0.14):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#    ###
#    elif args.property == 2: #LUMO,B3LYPs
#    
#        if args.y_flavor != 6: # for all y flavors except HF  #LUMO,B3LYPs
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.20) and i[0]<=(0.05):
#                    if i[1]>=(-0.20) and i[1]<=(0.05):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#        else: # only for HF (y flavor) #LUMO,B3LYPs
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.20) and i[0]<=(0.10):
#                    if i[1]>=(-0.20) and i[1]<=(0.30):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#    ###
#    elif args.property == 3: #Dipole,B3LYPs
#    
#        if args.y_flavor != 6: # for all y flavors except HF  #Dipole,B3LYPs
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.10) and i[0]<=(13.0):
#                    if i[1]>=(-0.10) and i[1]<=(13.0):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#        else: # only for HF (y flavor) #Dipole,B3LYPs
#            def lrboundaries(s,t):
#                #i[0] is x, i[1] is y, first term is min, second is max
#                if i[0]>=(-0.10) and i[0]<=(13.0):
#                    if i[1]>=(-0.10) and i[1]<=(13.0):
#                        return True
#                    else:
#                        list_outbound.append(index+3)
#                    return False
#                else:
#                    list_outbound.append(index+3)
#                    return False
#    ###
#    elif args.property == 4: #Gap,B3LYPs
#    
#        def lrboundaries(s,t):
#            #i[0] is x, i[1] is y, first term is min, second is max
#            if i[0]>=(0.0) and i[0]<=(0.20):
#                if i[1]>=(0.0) and i[1]<=(0.40):
#                    return True
#                else:
#                    list_outbound.append(index+3)
#                    return False
#            else:
#                list_outbound.append(index+3)
#                return False

elif  args.x_flavor == 3:
    x_name = "PBE0s"
elif  args.x_flavor == 4:
    x_name = "BHHLYPs"
elif  args.x_flavor == 5:
    x_name = "M06s"
elif  args.x_flavor == 6:
    x_name = "HFs"
elif  args.x_flavor == 7:
    x_name = "BP86t"
elif  args.x_flavor == 8:
    x_name = "B3LYPt"
elif  args.x_flavor == 9:
    x_name = "PBE0t"
    
if args.y_flavor == 1:
    y_name = "BP86s"
    y_fname = "BP86/SVP"
    y_num = "1"
elif  args.y_flavor == 2:
    y_name = "B3LYPs"
    y_fname = "B3LYP/SVP"
    y_num = "2"
elif  args.y_flavor == 3:
    y_name = "PBE0s"
    y_fname = "PBE0/SVP"
    y_num = "3"
elif  args.y_flavor == 4:
    y_name = "BHHLYPs"
    y_fname = "BH&HLYP/SVP"
    y_num = "4"
elif  args.y_flavor == 5:
    y_name = "M06s"
    y_fname = "M06/SVP"
    y_num = "5"
elif  args.y_flavor == 6:
    y_name = "HFs"
    y_fname = "HF/SVP"
    y_num = "6"
elif  args.y_flavor == 7:
    y_name = "BP86t"
    y_fname = "BP86/TZVP"
    y_num = "7"
elif  args.y_flavor == 8:
    y_name = "B3LYPt"
    y_fname = "B3LYP/TZVP"
    y_num = "8"
elif  args.y_flavor == 9:
    y_name = "PBE0t"
    y_fname = "PBE0/TZVP"
    y_num = "9"


ID_SMILE_X_Y_dis = "DFT_%s_%s%s_%s_%s.csv"% (pro_name, x_num,y_num, x_name, y_name)
num5 = str(int(args.percent_outlier * 100))
lr_text_outlier = "DFT_%s_%s%s_%s_%s_outlier_%s.csv" % (pro_name, x_num,y_num, x_name, y_name, num5)
lr_text_fitted = "DFT_%s_%s%s_%s_%s_fitted_%s.csv" % (pro_name, x_num,y_num, x_name, y_name, num5)
name_directory = "DFT_%s_%s_%s%s_%s_%s"% (pro_name,num5,x_num,y_num, x_name, y_name)
os.makedirs(name_directory)
invalid_compound = "DFT_invalid_%s_%s%s_%s_%s"% (pro_name, x_num,y_num, x_name, y_name)
outbound_compound = "DFT_outbound_%s_%s%s_%s_%s"% (pro_name, x_num,y_num, x_name, y_name)
info_name = "DFT_info_%s_%s%s_%s_%s.txt"% (pro_name, x_num,y_num, x_name, y_name)
lr_plot = "lr_%s_%s%s_%s_%s.png" % (pro_name,x_num,y_num,x_name, y_name)
lr_plot_outlier = "lr_%s_%s%s_%s_%s_outlier_%s.png" % (pro_name,x_num,y_num,x_name, y_name, num5)
lr_plot_remained = "lr_%s_%s%s_%s_%s_fitted_%s.png" % (pro_name, x_num,y_num,x_name, y_name, num5)
num_percentage = str(args.percent_outlier)
##############################################################################

  
slope, intercept, r_value, p_value, std_err = stats.linregress(list_x,list_y)

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
    

list_2 = zip(list_y,list_dis)
array2 = np.asarray(list_2)

list_percentage = []
for i_y1, i_dis in array2:
    num4 = abs(i_dis/i_y1) * 100
    list_percentage.append(round(num4,2))
    
x1_2 = []
y1_2 = []
for i in list_x:
    num5 = round(i,3)
    x1_2.append(num5)
for i in list_y:
    num6 = round(i,3)
    y1_2.append(num6)


list_3 = []
for i in range(0,len(list_x)):
    list_3.append([list_SMILES[i],x1_2[i],y1_2[i],y2[i],list_dis[i],list_percentage[i]])


str4 = "%s/%s" %(name_directory,ID_SMILE_X_Y_dis)
with open(str4, "wb") as f:
    writer = csv.writer(f)
    writer.writerows(list_3)
    

count_num = 0
list_4 = []
list_5 = []
x_outliers = []
y_outliers = []
x_fitted = []
y_fitted = []
for i,line in enumerate(list_percentage):
    if line >= args.percent_outlier:
        count_num = count_num + 1
        item_1 = [i+1] + list_3[i]
        list_4.append(item_1)
        x_outliers.append(x1_2[i])
        y_outliers.append(y1_2[i])
    if line < args.percent_outlier:
        item_2 = [i+1] + list_3[i]
        list_5.append(item_2)
        x_fitted.append(x1_2[i])
        y_fitted.append(y1_2[i])
        
str5 = "%s/%s" %(name_directory,lr_text_outlier)
with open(str5, "wb") as f:
    writer = csv.writer(f)
    writer.writerows(list_4)

str6 = "%s/%s" %(name_directory,lr_text_fitted)
with open(str6, "wb") as f:
    writer = csv.writer(f)
    writer.writerows(list_5)
    
str9 = "%s/%s" %(name_directory,info_name)
txt_info = open(str9, 'a')
txt_info.write('# of total compounds: '),txt_info.write(str(len(list_x))),txt_info.write('\n')
txt_info.write('# of outliers: '),txt_info.write(str(len(list_4))),txt_info.write('\n')
txt_info.write('# of fitted compounds: '),txt_info.write(str(len(list_5))),txt_info.write('\n')
txt_info.write('Linear regression line:'),txt_info.write('\n')
txt_info.write('slope = '),txt_info.write(str(slope)),txt_info.write('\n')
txt_info.write('intercept = '), txt_info.write(str(intercept))
txt_info.close()   


###########################Plot for all compounds##################
x1_3 = []
y1_3 = []
for i in list_x:
    if args.property == 3:  #only for Dipole
        x1_3.append(i)
    else:
        x1_3.append(i*27.21139570)
        
for i in list_y:
    if args.property == 3:  #only for Dipole
        y1_3.append(i)
    else:
        y1_3.append(i*27.21139570)

#linear regression using stats.linregress
slope, intercept, r_value, p_value, std_err = stats.linregress(x1_3,y1_3)
r_squared= r_value**2
#
#this will find your y-intercept and slope of your linear regression line and 
#then will create points where a line is plotted on
(m,b)=polyfit(x1_3,y1_3,1)
yp=polyval([m,b],x1_3)
yp1=yp
stand = (np.std(x1_3) + np.std(y1_3))/2
print 'r-squared: ', r_squared
print 'stand: ', stand

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['mathtext.default'] = 'regular'
fig= plt.figure()
ax = fig.add_subplot(111, aspect='equal')
ax.scatter(x1_3,y1_3,s=0.5,color='Black')


textstr = '$R^2:\ \ \ \ \ \ \  %.4f$\n$slope:\ \ \ %.4f$\n$shift:\ \ \ \ \ %.4f$'%(r_squared, m, b)
props = dict(boxstyle='round', facecolor= 'White', alpha=1.0)

if args.result_box == True and args.result_box1 == True:
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=props)
        
write_all(m,b,r_squared,pro_name,x_name,y_name)

if args.property == 1:
    ax.set_xlabel(r'$\epsilon\ $HOMO (%s) [eV]'%(x_fname), fontsize='x-large')
    ax.set_ylabel(r'$\epsilon\ $HOMO (%s) [eV]'%(y_fname), fontsize='x-large')
elif args.property == 2:
    ax.set_xlabel(r'$\epsilon\ $LUMO (%s) [eV]'%(x_fname), fontsize='x-large')
    ax.set_ylabel(r'$\epsilon\ $LUMO (%s) [eV]'%(y_fname), fontsize='x-large')
elif args.property == 3:
    ax.set_xlabel(r'$\mu\ $  (%s) [D]'%(x_fname), fontsize='x-large')
    ax.set_ylabel(r'$\mu\ $  (%s) [D]'%(y_fname), fontsize='x-large')
elif args.property == 4:
    ax.set_xlabel(r'$\Delta\ \epsilon\ $  (%s) [eV]'%(x_fname), fontsize='x-large')
    ax.set_ylabel(r'$\Delta\ \epsilon\ $  (%s) [eV]'%(y_fname), fontsize='x-large')

ax.tick_params(labelsize='large')
ax.set_xlim((num1),(num2))
ax.set_ylim((num1),(num2))

if args.lr_line == True and args.lr_line1 == True:
    ax = extended(ax, x1_3, yp, lw=2.4, color = line_color)

ax.grid(True)
str10 = "%s/%s" %(name_directory,lr_plot)
plt.savefig(str10)

###########################Plot for outliers############################
x_outliers_2 = []
y_outliers_2 = []
for i in x_outliers:
    if args.property == 3:  #only for Dipole
        x_outliers_2.append(i)
    else:
        x_outliers_2.append(i*27.21139570)
        
for i in y_outliers:
    if args.property == 3:  #only for Dipole
        y_outliers_2.append(i)
    else:
        y_outliers_2.append(i*27.21139570)

slope, intercept, r_value, p_value, std_err = stats.linregress(x_outliers_2,y_outliers_2)
r_squared= r_value**2
#
#this will find your y-intercept and slope of your linear regression line and 
#then will create points where a line is plotted on
(m2,b2)=polyfit(x_outliers_2,y_outliers_2,1)
yp2=polyval([m2,b2],x_outliers_2)
yp1=yp2
stand = (np.std(x_outliers_2) + np.std(y_outliers_2))/2
print 'r-squared: ', r_squared
print 'stand: ', stand

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['mathtext.default'] = 'regular'
fig= plt.figure()
ax = fig.add_subplot(111, aspect='equal')
ax.scatter(x_outliers_2,y_outliers_2,s=0.5,color='Black')


textstr = '<Outliers>\n$R^2:\ \ \ \ \ \ \  %.4f$\n$slope:\ \ \ %.4f$\n$shift:\ \ \ \ \ %.4f$'%(r_squared, m2, b2)
props = dict(boxstyle='round', facecolor= 'White', alpha=1.0)

if args.result_box == True and args.result_box2 == True:
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=props)
        
write_outlier(m2,b2,r_squared,pro_name,x_name,y_name,num_percentage)

if args.property == 1:
    ax.set_xlabel(r'$\epsilon\ $HOMO (%s) [eV]'%(x_fname), fontsize='x-large')
    ax.set_ylabel(r'$\epsilon\ $HOMO (%s) [eV]'%(y_fname), fontsize='x-large')
elif args.property == 2:
    ax.set_xlabel(r'$\epsilon\ $LUMO (%s) [eV]'%(x_fname), fontsize='x-large')
    ax.set_ylabel(r'$\epsilon\ $LUMO (%s) [eV]'%(y_fname), fontsize='x-large')
elif args.property == 3:
    ax.set_xlabel(r'$\mu\ $  (%s) [D]'%(x_fname), fontsize='x-large')
    ax.set_ylabel(r'$\mu\ $  (%s) [D]'%(y_fname), fontsize='x-large')
elif args.property == 4:
    ax.set_xlabel(r'$\Delta\ \epsilon\ $  (%s) [eV]'%(x_fname), fontsize='x-large')
    ax.set_ylabel(r'$\Delta\ \epsilon\ $  (%s) [eV]'%(y_fname), fontsize='x-large')

ax.tick_params(labelsize='large')
ax.set_xlim((num1),(num2))
ax.set_ylim((num1),(num2))

if args.lr_line == True and args.lr_line2 == True:
    ax = extended(ax, x_outliers_2, yp2, lw=2.4, color = line_color)
#ax.plot(x1,yp,'green',linewidth=1.0)
ax.grid(True)
str11 = "%s/%s" %(name_directory,lr_plot_outlier)
plt.savefig(str11)





























#####
pr.disable()

f = open('x.prof', 'a')
sortby = 'cumulative'
pstats.Stats(pr, stream=f).strip_dirs().sort_stats(sortby).print_stats()
f.close()
