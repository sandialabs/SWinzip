# -*- coding: utf-8 -*-
import sys, os, getopt
import numpy as np
import math

F = open('../data/mesh2D_irreg.txt','r')

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

# perform the transform  
ppi=np.pi
f=np.zeros(N)
for i in range(0,N):
  f[i]=(4*math.sin(8*ppi*p[i,0]) )*(4*math.sin(7*ppi*p[i,1]) )*3*math.sin(6*ppi*p[i,0])
  
N=len(f)
J=len(Uj)
P=2**(J-1)

G=np.zeros(P,dtype=int)
for i in range(0,P):
  G[i]=len(part[i])

z=np.array([1])  
si=np.concatenate((z, np.cumsum(G)))  

dire=1
j_list=np.zeros(J,dtype=int)
for j in range(0,J):
  j_list[j]=j+1
  
if (dire==-1):
  j_list=np.flip(j_list)
  
w=f.copy()

#for i in range(0,N):
  #print(f[i])

#print(len(si))  
#print(si)
for j in j_list:
  Ui = Uj[j-1]

  if (j==1):
    r=w.copy()
    for i in range(0,P):
      selj = part[i]
      offs = si[i]-i*k2
      lo = G[i]-k2
      
      seli=np.zeros(lo,dtype=int)
      for l in range(0,lo):
        seli[l]=l+offs-1
      
      #print(Ui[i].size)
      U = Ui[i][k2:, :]
      #print(U.size)
      if (dire==1):
        V=np.matmul(U, np.take(w,selj))
        for k in range(0,len(V)):
          r[seli[k]]=V[k]
      else:
        V=np.matmul(np.ndarray.transpose(U), np.take(w,seli))        
        for k in range(0,len(V)):
          r[selj[k]]=V[k]
          
      offs = N - P*k2 + 1 + i*k2
      if (offs>=0):
        seli=np.zeros(k2,dtype=int)
        for l in range(0,k2):
          seli[l]=l+offs-1
        U = Ui[i][0:k2, :]
        
        if (dire==1):
          V=np.matmul(U, np.take(w,selj))
          for k in range(0,len(V)):
            r[seli[k]]=V[k]
        else:
          V=np.matmul(np.ndarray.transpose(U), np.take(w,seli))
          for k in range(0,len(V)):
            r[selj[k]]+=V[k]
      else:
        print('Problem empty bin.')
    
    w=r.copy()

  else:
    nj = int(P/2**(j-1))
    mj = nj*2*k2;
    #print(nj,mj)
    r=w.copy()
    offs = N-mj
    
    for i in range(0,nj):
      
      selj=np.zeros(2*k2,dtype=int)
      for l in range(0,2*k2):
        selj[l]=offs+2*k2*i+l;
      
      seli=np.zeros(k2,dtype=int)
      for l in range(0,k2):
        seli[l]=offs+k2*i+l;

      U=Ui[i][k2:,:]
      if (dire==1):
        V=np.matmul(U, np.take(w,selj))        
        for k in range(0,len(V)):
          r[seli[k]]=V[k]
      else:
        V=np.matmul(np.ndarray.transpose(U), np.take(w,seli))
        for k in range(0,len(V)):
          r[selj[k]]=V[k]
          
      for l in range(0,k2):
        seli[l]=offs+mj/2+k2*i+l;
      
      U=Ui[i][0:k2,:]      
      if (dire==1):
        V=np.matmul(U, np.take(w,selj))
        for k in range(0,len(V)):
          r[seli[k]]=V[k]
      else:
        V=np.matmul(np.ndarray.transpose(U), np.take(w,seli))
        for k in range(0,len(V)):
          r[selj[k]]+=V[k]
          
    w=r.copy()
    
for i in range(0,N):
  print(r[i])