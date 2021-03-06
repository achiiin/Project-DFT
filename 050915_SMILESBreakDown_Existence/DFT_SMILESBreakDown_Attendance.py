#==============================================================================
# To break SMILES down into building blocks
#==============================================================================
import csv
import time
start = time.time()
list_SMILES = []

def record_time(part,start,timenow,timelast):
#    str1 = str(part)
    str2 = str(timenow - timelast)
    str3 = str(timenow - start)
    txt = open('time_record.txt', 'a')
    txt.write(part),txt.write("    "),txt.write(str2),txt.write("    "),
    txt.write(str3),txt.write("\n") 
    
with open("uni_SMILES_391.txt") as f:
#with open("uni_SMILES.txt") as f:
    
    for i in f:
        str_smiles = i.strip()
        list_SMILES.append(str_smiles)
f.close()
################################
list_index = [i for i in range(1,27)]
list_smart_incorrect = ['c1ccccc1','s1cccc1','[SiH2]1cccc1','[se]1cccc1','[nH]1cccc1', #1~5
              'o1cccc1', 'C1C=CC=C1','s1cc2c(c1)ccs2', 's1cc2c(c1)C=C[SiH2]2','s1cc2c(c1)cc[se]2',  #6~10
              's1cc2c(c1)cc[nH]2', 's1cc2c(c1)cco2', 's1cc2c(c1)C=CC2', 's1cc2c(c1)nccn2','s1nc2c(n1)cccc2',  #11~15
              'c1cccc2c1cccc2','n1ccccc1','c1ncncn1','c1cncs1', 'c12c(c3c(c4c1cccc4)cccc3)csc2',  #16~20
              'c1cccc2c1csc2','c1cccc2c1c[nH]c2','c1cccc2c1coc2','c1cccc2c1c[SiH2]c2','c1cccc2c1cCc2',  #21~25
              'c1c2c(cnc1)nsn2']  #26
list_cans = ['c1ccccc1', 'c1cccs1', '[SiH2]1C=CC=C1', 'c1ccc[se]1', 'c1ccc[nH]1',
             'c1ccco1', 'C1C=CC=C1', 's1cc2c(c1)ccs2', 'C1=Cc2c([SiH2]1)csc2', 's1cc2c(c1)cc[se]2',
             's1cc2c(c1)cc[nH]2', 's1cc2c(c1)cco2', 'C1=Cc2c(C1)csc2', 'c1cnc2c(n1)csc2', 'c1ccc2c(c1)nsn2',
             'c1ccc2c(c1)cccc2', 'c1cccnc1', 'n1cncnc1', 'c1nccs1', 'c1ccc2c(c1)c1ccccc1c1c2csc1', 
             'c1ccc2c(c1)csc2', 'c1ccc2c(c1)c[nH]c2', 'c1ccc2c(c1)coc2', '[SiH2]1[CH]c2c(C1)cccc2', 'C1[CH]c2c(C1)cccc2',
             'c1ncc2c(c1)nsn2']  
# without ~
list_smart_2 = ['[#6]1[#6][#6][#6][#6][#6]1','s1[#6][#6][#6][#6]1','[SiH2]1[#6][#6][#6][#6]1','[se]1[#6][#6][#6][#6]1','[nH]1[#6][#6][#6][#6]1', #1~5
              'o1[#6][#6][#6][#6]1', '[#6]1[#6][#6][#6][#6]1','s1[#6][#6]2[#6]([#6]1)[#6][#6]s2', 's1[#6][#6]2[#6]([#6]1)[#6][#6][SiH2]2','s1[#6][#6]2[#6]([#6]1)[#6][#6][se]2',  #6~10
              's1[#6][#6]2[#6]([#6]1)[#6][#6][nH]2', 's1[#6][#6]2[#6]([#6]1)[#6][#6]o2', 's1[#6][#6]2[#6]([#6]1)[#6][#6][#6]2', 's1[#6][#6]2[#6]([#6]1)n[#6][#6]n2','s1n[#6]2[#6](n1)[#6][#6][#6][#6]2',  #11~15
              '[#6]1[#6][#6][#6][#6]2[#6]1[#6][#6][#6][#6]2','n1[#6][#6][#6][#6][#6]1','[#6]1n[#6]n[#6]n1','[#6]1[#6]n[#6]s1', '[#6]12[#6]([#6]3[#6]([#6]4[#6]1[#6][#6][#6][#6]4)[#6][#6][#6][#6]3)[#6]s[#6]2',  #16~20
              '[#6]1[#6][#6][#6][#6]2[#6]1[#6]s[#6]2','[#6]1[#6][#6][#6][#6]2[#6]1[#6][nH][#6]2','[#6]1[#6][#6][#6][#6]2[#6]1[#6]o[#6]2','[#6]1[#6][#6][#6][#6]2[#6]1[#6][SiH2][#6]2','[#6]1[#6][#6][#6][#6]2[#6]1[#6][#6][#6]2',  #21~25
              '[#6]1[#6]2[#6]([#6]n[#6]1)nsn2']  
## no "~" in front of number
list_smart_3 = ['[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1','s1~[#6]~[#6]~[#6]~[#6]1','[SiH2r5]','[se]1~[#6]~[#6]~[#6]~[#6]1','[nH]1~[#6]~[#6]~[#6]~[#6]1', #1~5
              'o1~[#6]~[#6]~[#6]~[#6]1', '[#6]1~[#6]~[#6]~[#6]~[#6]1','s1~[#6]~[#6]2~[#6](~[#6]1)~[#6]~[#6]~s2', '[#6]1~[#6]~[#6]2~[#6](~[SiH2]1)~[#6]~s~[#6]2','s1~[#6]~[#6]2~[#6](~[#6]1)~[#6]~[#6]~[se]2',  #6~10
              's1~[#6]~[#6]2~[#6](~[#6]1)~[#6]~[#6]~[nH]2', 's1~[#6]~[#6]2~[#6](~[#6]1)~[#6]~[#6]~o2', 's1~[#6]~[#6]2~[#6](~[#6]1)~[#6]~[#6]~[#6]2', 's1~[#6]~[#6]2~[#6](~[#6]1)~n~[#6]~[#6]~n2','s1~n~[#6]2~[#6](~n1)~[#6]~[#6]~[#6]~[#6]2',  #11~15
              '[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]1~[#6]~[#6]~[#6]~[#6]2','n1~[#6]~[#6]~[#6]~[#6]~[#6]1','[#6]1~n~[#6]~n~[#6]~n1','[#6]1~[#6]~n~[#6]~s1', '[#6]12~[#6](~[#6]3~[#6](~[#6]4~[#6]1~[#6]~[#6]~[#6]~[#6]4)~[#6]~[#6]~[#6]~[#6]3)~[#6]~s~[#6]2',  #16~20
              '[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]1~[#6]~s~[#6]2','[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]1~[#6]~[nH]~[#6]2','[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]1~[#6]~o~[#6]2','[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]1~[#6]~[SiH2]~[#6]2','[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]1~[#6]~[#6]~[#6]2',  #21~25
              '[#6]1~[#6]2~[#6](~[#6]~n~[#6]1)~n~s~n2']
              
list_smart_4 = ['[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1','s~1~[#6]~[#6]~[#6]~[#6]~1','[SiH2r5]','[se]~1~[#6]~[#6]~[#6]~[#6]~1','[nH]~1~[#6]~[#6]~[#6]~[#6]~1', #1~5
              'o~1~[#6]~[#6]~[#6]~[#6]~1', '[#6]~1~[#6]~[#6]~[#6]~[#6]~1','s~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~s~2', '[#6]~1~[#6]~[#6]~2~[#6](~[SiH2]~1)~[#6]~s~[#6]~2','s~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[se]~2',  #6~10
              's~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[nH]~2', 's~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~o~2', 's~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[#6]~2', 's~1~[#6]~[#6]~2~[#6](~[#6]~1)~n~[#6]~[#6]~n~2','s~1~n~[#6]~2~[#6](~n~1)~[#6]~[#6]~[#6]~[#6]~2',  #11~15
              '[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[#6]~[#6]~[#6]~2','n~1~[#6]~[#6]~[#6]~[#6]~[#6]~1','[#6]~1~n~[#6]~n~[#6]~n~1','[#6]~1~[#6]~n~[#6]~s~1', '[#6]~1~2~[#6](~[#6]~3~[#6](~[#6]~4~[#6]~1~[#6]~[#6]~[#6]~[#6]~4)~[#6]~[#6]~[#6]~[#6]~3)~[#6]~s~[#6]~2',  #16~20
              '[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~s~[#6]~2','[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[nH]~[#6]~2','[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~o~[#6]~2','[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[SiH2]~[#6]~2','[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[#6]~[#6]~2',  #21~25
              '[#6]~1~[#6]~2~[#6](~[#6]~n~[#6]~1)~n~s~n~2']

list_smart = ['[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1','[#16]~1~[#6]~[#6]~[#6]~[#6]~1','[#14]~1~[#6]~[#6]~[#6]~[#6]~1','[#34]~1~[#6]~[#6]~[#6]~[#6]~1','[#7]~1~[#6]~[#6]~[#6]~[#6]~1', #1~5
              '[#8]~1~[#6]~[#6]~[#6]~[#6]~1', '[#6]~1~[#6]~[#6]~[#6]~[#6]~1','[#16]~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[#16]~2', '[#6]~1~[#6]~[#6]~2~[#6](~[#14]~1)~[#6]~[#16]~[#6]~2','[#16]~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[#34]~2',  #6~10
              '[#16]~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[#7]~2', '[#16]~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[#8]~2', '[#16]~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#6]~[#6]~[#6]~2', '[#16]~1~[#6]~[#6]~2~[#6](~[#6]~1)~[#7]~[#6]~[#6]~[#7]~2','[#16]~1~[#7]~[#6]~2~[#6](~[#7]~1)~[#6]~[#6]~[#6]~[#6]~2',  #11~15
              '[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[#6]~[#6]~[#6]~2','[#7]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1','[#6]~1~[#7]~[#6]~[#7]~[#6]~[#7]~1','[#6]~1~[#6]~[#7]~[#6]~[#16]~1', '[#6]~1~2~[#6](~[#6]~3~[#6](~[#6]~4~[#6]~1~[#6]~[#6]~[#6]~[#6]~4)~[#6]~[#6]~[#6]~[#6]~3)~[#6]~[#16]~[#6]~2',  #16~20
              '[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[#16]~[#6]~2','[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[#7]~[#6]~2','[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[#8]~[#6]~2','[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[#14]~[#6]~2','[#6]~1~[#6]~[#6]~[#6]~[#6]~2~[#6]~1~[#6]~[#6]~[#6]~2',  #21~25
              '[#6]~1~[#6]~2~[#6](~[#6]~[#7]~[#6]~1)~[#7]~[#16]~[#7]~2']

#dic_smart = dict(zip(list_index, list_smart))
dic_smart = dict(zip(list_index, list_smart))
##################################
import pybel


dic_numbb = {}
for i in dic_smart:

    if i > 5:
        break
    for smiles in list_SMILES:        
        mol = pybel.readstring("smi", smiles)
        smarts = pybel.Smarts(dic_smart[i].strip())
         # pybel doesn't have a direct way to get the non-unique matches, so
         # use the lower-level OpenBabel API directly
        smarts.obsmarts.Match(mol.OBMol)
        #==============================================================================
        # num_matches = sum(1 for indicies in smarts.obsmarts.GetMapList())
        # print "number of matches:", num_matches 
        #==============================================================================
        num_unique_matches = len(smarts.findall(mol))
#        print smiles
#        print dic_smart[i]
#        print "number of unique matches:", num_unique_matches
        if num_unique_matches > 0:
            num_unique_matches = 1
        dic_numbb.setdefault(smiles,[]).append(num_unique_matches)
        

 
def saveDict(fn,dict_rap):
    f=open(fn, "wb")
    w = csv.writer(f)
    for key, val in dict_rap.items():
        w.writerow([key, val])
    f.close()

saveDict('CountBB_0303.csv',dic_numbb)

#####
end1 = time.time()
record_time("Part1",start,end1,start)
#####        
        
        
        
        
        
        
        
        
        