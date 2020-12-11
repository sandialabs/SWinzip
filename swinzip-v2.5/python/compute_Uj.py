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


# dichotomic grouping
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

k=5
ki=np.tile(k,d)
ptmax=int(np.prod(ki))
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
      
    I=np.argsort(tmp,kind='mergesort')
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



# computing Uj
P=len(part)
J=int(np.log2(P)+1);
k2=ptmax

#for i in range(0,len(part[100])):
  #print(part[100][i])

iid=np.zeros(d,dtype=int)
cnt=0
table=np.zeros((d,k2),dtype=int)

while (cnt<k2):
  for j in range(0,d):
    table[j,cnt] = iid[j]
  
  iid[0]+=1;
  i = 0
  while (i<(d-1)) and (iid[i]==ki[i]):
    iid[i]=0
    i+=1
    iid[i]+=1
  
  cnt+=1

Mi=[]
Ui=[]
for i in range(0,P):
  seli=part[i]
  kk=len(seli)
  
  M=np.zeros((kk,2*k2))
  for ii in range(0,kk):
    for j in range(0,k2):
      M[ii,j]=1
      for k in range(0,d):
        M[ii,j]*=p[seli[ii],k]**table[d-k-1,j]

  Mi.append(M)
  U=np.linalg.qr(M)
  Ui.append(np.ndarray.transpose(U[0]))
  
Uj=[]
Uj.append(Ui)
for j in range(2,J+1):
  nj = int(P/2**(j-1))
  
  MMi=[]
  UUi=[]
  for i in range(0,nj):
    M=np.zeros((2*k2,2*k2))
    M[0:k2,:]       = np.matmul(Ui[2*i][0:k2,:]  ,   Mi[2*i])
    M[k2:2*k2,:]    = np.matmul(Ui[2*i+1][0:k2,:],   Mi[2*i+1])
    MMi.append(M)
    U=np.linalg.qr(M)
    UUi.append(np.ndarray.transpose(U[0]))
   
  Mi=[]
  Ui=[]
  for i in range(0,2**(J-j)):
    Mi.append(MMi[i])
    Ui.append(UUi[i])
    
  Uj.append(Ui)

#tt=Uj[10][0]
#for i in range(0,50):
  #for j in range(0,50):
    #print(tt[j,i])
