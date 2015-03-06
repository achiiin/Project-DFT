# -*- coding: utf-8 -*-
"""
Created on Thu Mar 05 18:13:45 2015

@author: Ching-Yen
"""
import csv

filename_CountBB = "BB_test.csv"


def readDict(fn):
    f=open(fn,'rb')
    dict_rap={}
     
    for key, val in csv.reader(f):
        dict_rap[key]=eval(val)
    f.close()
    return(dict_rap)
  
dict_BB = readDict(filename_CountBB)

def substract_overcount(dic_1):
    import copy
    dic_new = copy.deepcopy(dic_1)
    for key,lists in dict.items(dic_new):
        #bb no.8
        lists[1] = lists[1] - lists[7] *2
        #bb no.9
        lists[1] = lists[1] - lists[8]
        lists[2] = lists[2] - lists[8]
        #bb no.10
        lists[1] = lists[1] - lists[9]
        lists[3] = lists[3] - lists[9]
        #bb no.11
        lists[1] = lists[1] - lists[10]
        lists[4] = lists[4] - lists[10]
        #bb no.12
        lists[1] = lists[1] - lists[11]
        lists[5] = lists[5] - lists[11]
        #bb no.13
        lists[1] = lists[1] - lists[12]
        lists[6] = lists[6] - lists[12]
        #bb no.14
        lists[1] = lists[1] - lists[13]
        #bb no.15
        lists[0] = lists[0] - lists[14]
        #bb no.16
        lists[0] = lists[0] - lists[15]*2
        #bb no.20
        lists[0] = lists[0] - lists[19]*3 + lists[19]*5 
        lists[1] = lists[1] - lists[19] + lists[19] *1 
        lists[15] = lists[15] - lists[19]*2
        lists[20] = lists[20] - lists[19]
        #bb no.21
        lists[0] = lists[0] - lists[20]
        lists[1] = lists[1] - lists[20]
        #bb no.22
        lists[0] = lists[0] - lists[21]
        lists[4] = lists[4] - lists[21]
        #bb no.23
        lists[0] = lists[0] - lists[22]
        lists[5] = lists[5] - lists[22]
        #bb no.24
        lists[0] = lists[0] - lists[23]
        lists[2] = lists[2] - lists[23]
        #bb no.25
        lists[0] = lists[0] - lists[24]
        lists[6] = lists[6] - lists[24]
        #bb no.26
        lists[16] = lists[16] - lists[25]

    
    
    
    return dic_new
    
    

dic_new = substract_overcount(dict_BB)

print dict_BB["C1C(=Cc2c1c1ncc3C=C[SiH2]c3c1c1=C[SiH2]C=c21)c1cccc2ccccc12"]
print dic_new["C1C(=Cc2c1c1ncc3C=C[SiH2]c3c1c1=C[SiH2]C=c21)c1cccc2ccccc12"]