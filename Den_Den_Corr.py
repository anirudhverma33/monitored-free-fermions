#!/usr/bin/env python
# coding: utf-8

# In[44]:


import numpy as np
import scipy.linalg as sci
from scipy.linalg import expm as exp
import matplotlib.pyplot as plt
import tqdm as tqdm


# In[42]:


ti=0
tf=10
del_t=0.1
t=np.arange(ti,tf+del_t,del_t)
L=50
h=np.zeros((L,L))
p=1
for i in range(L):
    for j in range(L):
        if (i==j-1 or i==j+1):
            h[i,j]=p
h[0,L-1]=p
h[L-1,0]=p
pos=1
n=np.diag(np.array([(i+1)%2 for i in range(L)]))
I=np.identity(L)
corr_h=np.zeros((len(t),len(t)))
for i in range(len(t)):
    for j in range(len(t)):
        corr_h[i,j]=np.dot(exp(-1j*h*t[j]),np.dot(n,exp(1j*h*t[i])))[pos,pos]*np.dot(exp(-1j*h*t[i]),np.dot((I-n),exp(1j*h*t[j])))[pos,pos]+np.dot(exp(-1j*h*t[i]),np.dot(n,exp(1j*h*t[i])))[pos,pos]*np.dot(exp(-1j*h*t[j]),np.dot(n,exp(1j*h*t[j])))[pos,pos]


# In[74]:


dt=0.1
t1=np.arange(ti,tf+dt,dt)
corr_u=np.zeros((len(t1),len(t1)))
N=int(L/2)
U0=np.zeros((L,N))
for i in range(N):
    U0[2*i,i]=1
Q0,R0=sci.qr(U0)
Q1=Q0.copy()
Q2=Q0.copy()
for i in range(len(t1)):
    U1=Q1[:,:N]
    for j in range(len(t1)):
        if (j==0):
            Q2=Q0.copy()
        U2=Q2[:,:N]
        corr_u[i,j]=np.dot(U1,np.transpose(np.conjugate(U1)))[pos,pos]*np.dot(U2,np.transpose(np.conjugate(U2)))[pos,pos] + np.dot(U2,np.transpose(np.conjugate(U1)))[pos,pos]*np.dot(Q1,np.transpose(np.conjugate(Q2)))[pos,pos] - np.dot(U2,np.transpose(np.conjugate(U1)))[pos,pos]*np.dot(U1,np.transpose(np.conjugate(U2)))[pos,pos]
        Q2=np.dot(exp(-1j*h*dt),Q2)
    Q1=np.dot(exp(-1j*h*dt),Q1)


# In[79]:


plt.scatter((t-5),corr_h[:,list(t).index(5)],marker='o',label='hermitian case',color='green')
plt.plot((t1-5),corr_u[:,list(t1).index(5)],label='g=0 case')
plt.title("Benchmarking with t'=5")
plt.xlabel('t-t')
plt.ylabel('D')
plt.legend()


# In[14]:


plt.imshow(corr_h,interpolation='spline36',origin='lower')
plt.show()


# In[8]:


plt.imshow(corr_u,interpolation='spline36',origin='lower')
plt.show()


# In[26]:


corr_u_disc=np.zeros((len(t),len(t)))
for i in range(len(t)):
    for j in range(len(t)):
        corr_u_disc[i,j]=corr_u[list(t1).index(t[i]),list(t1).index(t[j])]


# In[33]:


plt.scatter (t,np.array([corr_h[i,i] for i in range(len(t))]),marker='o',label='Hermitian case')
plt.scatter (t,np.array([corr_u_disc[i,i] for i in range(len(t))]),marker='.',label='g=0 case')
plt.title('Benchmarking the Hermitian case')
plt.xlabel('time')
plt.ylabel('Diagonal terms')
-plt.legend()


# In[37]:


plt.plot (t1,np.array([corr_u[i,i] for i in range(len(t1))]),label='g=0 case')
plt.title("Variation of D with (t+t')/2, for (t-t')/2=0")
plt.xlabel("(t+t')/2")
plt.ylabel('D')
plt.legend()


# In[38]:


len(t1)


# In[41]:


plt.plot (np.array([(t1[i]-t1[len(t1)-1-i])/2 for i in range(len(t1))]),np.array([corr_u[i,1000-i] for i in range(len(t1))]),label='g=0 case')
plt.title("Variation of D with (t-t')/2, for (t+t')/2=50")
plt.xlabel("(t-t')/2")
plt.ylabel('D')
plt.legend()


# In[ ]:




