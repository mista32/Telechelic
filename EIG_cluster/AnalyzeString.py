import os, sys, glob,ctypes
import numpy as np
from subprocess import call

def PL(number_m,XY):
    ##calculate the probability of the length of the string
    # X = np.zeros(number_m)

    X = XY[:,0]
    # X[len(XY[:,0]):] = 1

    high = round(max(X))
    
    low = 2
    N = round(high-low)+1
    Y = np.zeros(N)
    for i in range(N):
        for j in range(len(X)):
            if (abs(X[j]-(i+low))<0.01):
                Y[i] = Y[i] + 1
    Y = Y/len(X)
    DataD = np.zeros((N,2))
    DataD[:,0] = np.arange(2,high+1,1)
    DataD[:,1] = Y
    np.savetxt('PL_Cluster',DataD)

def Rg_L(XY):
    X = XY[:,0]

    X = np.unique(X)
    Y = np.zeros((len(X),4))
    for i in range(len(X)):
        for j in range(4):
            Y[i,j] = np.mean(XY[(abs(XY[:,0]-X[i]))<0.01,j+1])

    XY = np.array([[X[i],Y[i,0],Y[i,1],Y[i,2],Y[i,3] ] for i in range(len(X))])
    np.savetxt('Rg_Cluster',XY)

if __name__ == '__main__':
    CODE_PATH ='/public2/home/yucao/research/Telechelic/SH050_T1500/03Melt/DC028/Pcluster'

    n=0
    with open('../001/Conf/Basic.txt') as f:
        for line in f:
            data = list(map(float, line.split()))
            n = n + 1
            if n == 4:
                number = int(data[0])
            if n == 6:
                box = float(data[0])

    number_m = number
############################################################################################################
    Rg = ctypes.CDLL(f'{CODE_PATH}/analyze/Rg.so')

    Rg.RG(ctypes.c_int(number_m),ctypes.c_longdouble(box))

###sort the data
    LRg =  np.loadtxt('Cluster_eig');  os.remove('Cluster_eig')
    LRg = LRg[np.lexsort(LRg[:,::-1].T)]
    
    PL(number_m,LRg)
    Rg_L(LRg)          


