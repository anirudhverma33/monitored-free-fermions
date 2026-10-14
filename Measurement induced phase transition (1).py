#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import scipy
import scipy.linalg
import random
import math
import matplotlib.pyplot as plt


# In[3]:


def normalize(U):
    V=np.empty((L,N),dtype=np.complex_)
    for j in range(N):
        V[:,j]=U[:,j]
        for k in range(j):
            V[:,j]=V[:,j]-np.vdot(U[:,j],V[:,k])*V[:,k]
        V[:,j]=V[:,j]/np.linalg.norm(V[:,j])
    return V  


# In[57]:


L=10
length=int(L/2)
N=int(L/2)


# In[58]:


h=np.zeros((L,L))
p=1
for i in range(L):
    for j in range(L):
        if (i==j-1 or i==j+1):
            h[i,j]=p
print(h)


# In[32]:


dt=0.5
time=np.arange(0,10,dt)
w=11
gamma=np.linspace(0,1,w)
A=100
S_avg=np.zeros((w,length))
f=np.zeros(length)
for l in range(length):
    f[l]=math.sin(math.pi*l/L)
for g in range(w):
    print(g)
    for l in range(length):
        S_tot=0
        D_tot=np.zeros((L,L),dtype=np.complex_)
        for a in range(A):
            U=np.zeros((L,N),dtype=np.complex_)
            for k in range(N):
                U[2*k+1,k]=1+1j*0
            for t in time:
                Ut=np.conjugate(np.transpose(U))
                D=np.transpose(np.dot(U,Ut))       
                n=np.zeros(L,dtype=np.complex_)
                for j in range(L):
                    n[j]=D[j,j]
                M=np.zeros((L,L),dtype=np.complex_)
                for i in range(L):
                    M[i,i]=random.gauss(mu=0, sigma=np.sqrt(gamma[1]*dt))/dt+gamma[1]*(2*n[i]-1)/2
                h1=h+1j*M
                U=np.dot(scipy.linalg.expm(-1j*dt*h1),U)
                U=normalize(U)
            D=np.transpose(np.dot(U,np.conjugate(np.transpose(U))))
            D_red=np.zeros((l,l),dtype=np.complex_)
            for i in range(l):
                for j in range(l):
                    D_red[i,j]=D[i,j]
            ev=np.round(np.linalg.eig(D_red)[0].real,3)
            S1=0
            for i in range(l):
                if(ev[i]>=1):
                    ev[i]=1-0.001
                if(ev[i]<=0):
                    ev[i]=0.001
                S1=S1-ev[i]*math.log2(ev[i])-(1-ev[i])*math.log2(1-ev[i])
            S_tot=S_tot+S1
        S_avg[g,l]=S_tot/A
    
for g in range(w):
    plt.plot(f,S_avg[g,:])


# In[1]:


#Benchmark Test


# In[90]:


def benchmark(L,t1,t2,dt):
    length=int(L/2)
    N=int(L/2)
    h=np.zeros((L,L))
    p=1
    for i in range(L):
        for j in range(L):
            if (i==j-1 or i==j+1):
                h[i,j]=p
    S=[]
    S1=[]
    tf=np.arange(t1,t2,dt)
    std=0
    for q in tf:
        time=np.arange(0,q,dt)
        U=np.zeros((L,N),dtype=np.complex_)
        U1=np.zeros((L,N),dtype=np.complex_)
        for k in range(N):
            U[2*k+1,k]=1+1j*0
            U1[2*k+1,k]=1+1j*0
        for t in time:
    #         Ut=np.conjugate(np.transpose(U))
    #         D=np.transpose(np.dot(U,Ut))       
    #         n=np.zeros(L,dtype=np.complex_)
    #         for j in range(L):
    #             n[j]=D[j,j]

        #                 M=np.zeros((L,L),dtype=np.complex_)
        #                 for i in range(L):
        #                     M[i,i]=random.gauss(mu=0, sigma=np.sqrt(gamma[g]*dt))/dt+gamma[g]*(2*n[i]-1)/2
        #                 h1=h+1j*M
            U=np.dot(scipy.linalg.expm(-1j*dt*h),U)
        Ut=np.conjugate(np.transpose(U))
        D=np.transpose(np.dot(U,Ut))
        D_red=np.zeros((length,length),dtype=np.complex_)
        for i in range(length):
            for j in range(length):
                D_red[i,j]=D[i,j]
        ev=np.round(np.linalg.eig(D_red)[0].real,3)
        tot=0
        for i in range(length):
            if(ev[i]>=1):
                ev[i]=1-0.001
            if(ev[i]<=0):
                ev[i]=0.001
            tot=tot-ev[i]*math.log2(ev[i])-(1-ev[i])*math.log2(1-ev[i])
        S.append(tot)
        U1=np.dot(scipy.linalg.expm(-1j*q*h),U1)
        Ut=np.conjugate(np.transpose(U1))
        D=np.transpose(np.dot(U1,Ut))
        D_red=np.zeros((length,length),dtype=np.complex_)
        for i in range(length):
            for j in range(length):
                D_red[i,j]=D[i,j]
        ev=np.round(np.linalg.eig(D_red)[0].real,3)
        tot=0
        for i in range(length):
            if(ev[i]>=1):
                ev[i]=1-0.001
            if(ev[i]<=0):
                ev[i]=0.001
            tot=tot-ev[i]*math.log2(ev[i])-(1-ev[i])*math.log2(1-ev[i])
        S1.append(tot)
    for i in range(len(S)):
        std=std+(S[i]-S1[i])**2
    std=np.sqrt(std/len(S))
    return std
        
#     plt.plot(tf,np.array(S),label="discrete")
#     plt.plot(tf,np.array(S1), label="continuous")
#     plt.xlabel("time")
#     plt.ylabel("Entanglement Entropy")
#     plt.title(f"time steps = {dt}")
#     plt.legend() 

         


# In[139]:


benchmark(0.5)


# In[147]:


def find_dt(t_fin,L):
    dt_arr=np.linspace(0.01,0.1,50)
    t_in=1
    std=np.zeros(len(dt_arr))
    # std_arr2=np.zeros(50)
    for i in range(len(dt_arr)):
        std[i]=benchmark(L,t_in,t_fin,dt_arr[i])
        
    plt.plot(dt_arr,std,label=f"tf={t_fin}, L={L}")
    plt.legend()
    
# plt.plot(t_fin,std_arr1, label="L=4")
# plt.plot(t_fin,std_arr2, label="L=6")
# plt.xlabel("tf")
# plt.ylabel("std")
# plt.title("dt=0.03, ti=1")
# plt.legend()


# In[148]:


find_dt(10,10)
find_dt(10,20)


# In[149]:


find_dt(10,10)
find_dt(20,10)


# In[ ]:





# In[ ]:




