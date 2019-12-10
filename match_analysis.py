import random
import matplotlib.pyplot as plt
import time
import numpy as np
import scipy.stats as st
t0 = time.clock()



##############################################


def Read_dataset(infile):
  D=dict()
  inp = open(infile,'r').readlines()
  el=set()
  for line in inp:
    line = line[:-1].split('\t')
    try:
      D[line[0]].add(line[1])
    except KeyError:
      D[line[0]] = set([line[1]])
    el.add(line[1])
  
  print('Number of unique elements: %d' % len(el))
  print('Number of clusters:        %d' % len(D))
  
  clus = D.values()
  clus_len = [len(i) for i in clus]
  
  print('Cluster lenght (min-max):  %d-%d' % (min(clus_len), max(clus_len)))
  print('Cluster lenght (mean):     %.3f' % (sum(clus_len)/len(clus_len)))
  return (clus,el)



def Search_match(d1,d2,n):
  c=0
  for i in d1:
    for j in d2:
      if len(i&j) >= n:
        c+=1
        break
  return c



def Random_search(rep,pool,dat1,dat2,el):
  M=[]
  shape = [len(i) for i in dat1]
  for r in range(rep):
    dat3=[]
    for n in shape:
      samp = random.sample(pool,n)
      samp = set(samp) & el
      if len(samp) < 2:
        continue
      dat3.append(samp)
    match = Search_match(dat3,dat2,2)
    M.append(match)
  return M
  


def Print_stat(match, observed):
  A = np.array(match)
  print('min:', np.amin(A))
  print('max:', np.amax(A))
  mean = np.mean(A)
  print('mean:', np.mean(A))
  deviation = np.std(A)
  print('std:', deviation)
  Z = (observed-mean)/deviation
  print('Z:', Z)
  p_val = st.norm.sf(abs(Z))*2
  print('p-val:', p_val)



##############################################


##############################################

## GENERAL ANALYSIS

print('\n-->  PS dataset analysis  <--')   
PS = Read_dataset('PS.txt') # analysis on the PS dataset

print('\n-->  PC dataset analysis  <--')
PC = Read_dataset('PC1.txt') # analysis on the PC dataset

print('\n-->  OMIM dataset analysis  <--')
omim = Read_dataset('omim.txt')

print('\n-->  LH dataset analysis  <--')
LH = Read_dataset('LH.txt')

omim2 = (omim[1] - PS[1]) - LH[1]

print (len(omim2))


##############################################


print('\n-->  Shared Proteins  <--')
int1 = PS[1] & PC[1] # PS and PC intersection (number of proteins)
print('Proteins shared between PS and PC:         %d' % len(int1))

int2 = omim2 & PC[1] # omim-no-LH and PC intersection (number of proteins)
print('Proteins shared between OMIM-no-LH and PC: %d' % len(int2))

# We want obtain reduced PS and PC dataset for searching faster.
# We selected the PC and PS clusters that have only shared proteins.

PS1 = [i&int1 for i in PS[0] if len(i&int1)>=2] # reduced PS dataset with only proteins present in PC
PC1 = [i&int1 for i in PC[0] if len(i&int1)>=2] # reduced PC dataset with only proteins present in PS
PC2 = [i&int2 for i in PC[0] if len(i&int2)>=2] # reduced PC dataset with only proteins present in omim-no-PS

el_PS1=set([i for clus in PS1 for i in clus])
el_PC1=set([i for clus in PC1 for i in clus])
el_PC2=set([i for clus in PC2 for i in clus])


match1 = Search_match(PS1,PC1,2) # observed match is 95/319
match2 = Search_match(PC1,PS1,2) # observed match is 100/575
print('Number of PS cluster sharing at least 2 proteins with at least 1 PC cluster: %d' % match1)
print('Number of PC cluster sharing at least 2 proteins with at least 1 PS cluster: %d' % match2)


M=Random_search(100000,PS[1],PS[0],PC1,el_PC1) # 1
#M=Random_search(100000,PC[1],PC[0],PS1,el_PS1) # 2
#M=Random_search(100000,omim2,PS[0],PC2,el_PC2) # 3

Print_stat(M, match1) # for simulation 1 and 3
#Print_stat(M, match2) # for simulation 2


t1 = time.clock()
time = t1-t0
print("Running time: %.3f seconds" % time)


out = open('matches.txt','w')
for num in M:
  out.write('%s\n' % num)
out.close()



x=list(set(M))
x.sort()
y=[M.count(n) for n in x]
plt.plot(x,y,'.-')
plt.title('n = %d' % len(M))
plt.show()

