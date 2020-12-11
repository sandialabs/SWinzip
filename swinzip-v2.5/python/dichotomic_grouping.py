import sys, os, getopt
import numpy as np
import math

F = open('mesh2D_irreg.txt','r')

N=F.readline()
N, d = N.split()
N=int(N)
d=int(d)

p=np.zeros((N,d))
for i in range(0, N):
  a=F.readline()
  a=a.split()
  #print(float(a[0]),float(a[1]))
  for j in range(0, d):
    p[i,j]=float(a[j])
    
F.close()

cpt=0;
depthmax=1000000;

IP=[];
part=[]
part_prev=[]
Ni=np.zeros(N, dtype=int)
for i in range(0, N):
  Ni[i]=i
part.append(Ni)
done=0
G=[]

ptmax=25
tp=d

while (cpt<depthmax) and (done==0):
  cpt+=1
  part_prev = part.copy()
  part1 = part.copy()
  part.clear()
  
  for k in range(0, len(part1)):
    P = part1[k]
    nP=len(P)
    pp = int(np.floor(nP/2))
        
    s = (cpt-1)%tp
    I=np.zeros(nP, dtype=int)
    tmp=[]
    
    for j in range(0,nP):
      tmp.append(p[P[j],s])
      I[j]=j
      
    I=np.argsort(tmp)
    ppart1 = I[0:pp];
    ppart2 = I[pp:nP];
    
    IP=np.zeros(len(ppart1),dtype=int)
    for j in range(0,len(ppart1)):
      IP[j]=P[ppart1[j]]
    part.append(IP)
    
    IP=np.zeros(len(ppart2),dtype=int)
    for j in range(0,len(ppart2)):
      IP[j]=P[ppart2[j]]
    part.append(IP)
    
  for j in range(0,len(part)):
    if (len(part[j])<ptmax):
      part = part_prev.copy()
      done=1
      break
      
for i in range(0,len(part)):
  G.append(len(part[i]))


