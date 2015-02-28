# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 17:44:07 2015

@author: Ching-Yen
"""

num_Si_lowecase=0
num_se_uppercase = 0

f = open('uni_SMILES.txt')
for smiles in f:
    num_Si_lowecase = num_Si_lowecase + str(smiles).count('si')
    num_se_uppercase = num_se_uppercase + str(smiles).count('Se')
f.close()
print "si:  ", num_Si_lowecase
print "Se:  ", num_se_uppercase
