__author__ = 'Yule'
import timeit
import time
import random
from math import log
import numpy as np

from Modularization_SL_Corre_Uncorre import BreakingRateVsStress

# L = 20
# N = L ** 2
# N_e = 2 * N - 2 * L
# k = 1
# T = 298
# Ea = 117.15 * 10 ** 3 / (6.02214 * 10 ** 23)
# nu = 418 * 10 ** (-30)
# tau0 = 10 ** (-11)
# sigma0 = 60E6





# 6%4=2 (横轴index) 6//4=1  (纵轴index）

import math
import os
import pickle



def BondNoneList(L):
    BondNone = []
    for i in range(1, L + 1):
        BondNone.append(2 * i * L - 1)
    for j in range(0, L):
        BondNone.append(2 * (L - 1) * L + 2 * j)
    return BondNone




def DistanceMatrix(L, N, BondNone, ratio=1):  # One time cal   #这个distance matrix只适合bond perco在SL
    # Distance:每一行代表某一个bond和其他所有bond的距离
    Distance = np.empty((2 * N, 2 * N))
    for cell_i in range(0, N):
        for cell_j in range(0, N):
            Distance[2 * cell_i][2 * cell_j] = (
                                               (cell_j % L - cell_i % L) ** 2 + (ratio * (cell_j // L) - ratio * (cell_i // L)) ** 2) ** 0.5
            Distance[2 * cell_i + 1][2 * cell_j + 1] = ((cell_j % L - cell_i % L) ** 2 + (
            ratio * (cell_j // L) - ratio * (cell_i // L)) ** 2) ** 0.5
            Distance[2 * cell_i][2 * cell_j + 1] = ((cell_j % L - cell_i % L + 0.5) ** 2 + (
            ratio * (cell_j // L) - ratio * (cell_i // L) - 0.5) ** 2) ** 0.5
            Distance[2 * cell_i + 1][2 * cell_j] = ((cell_j % L - cell_i % L - 0.5) ** 2 + (
            ratio * (cell_j // L) - ratio * (cell_i // L) + 0.5) ** 2) ** 0.5
    for cell_i in range(0, N):  # Delete self-self distance
        Distance[2 * cell_i][2 * cell_i] = np.nan
        Distance[2 * cell_i + 1][2 * cell_i + 1] = np.nan
    for i in BondNone:  # Delete None bond distance
        for j in range(0, 2 * N):
            Distance[i][j] = np.nan
    for i in range(0, 2 * N):
        for j in BondNone:
            Distance[i][j] = np.nan
    return Distance





def StressTransferFunction(gamma, Distance):  # One time cal
    # Distance: DistanceMatrix
    # F matrix
    F = Distance ** (-gamma)

    # Summation 在 StressList（）
    return F






def StressInit(stress, N,BondNone):
    Stress = np.full(2 * N, stress)
    #No stress on non-existed bonds:暂时先不这样
    for i in BondNone:
        Stress[i] = np.nan
    return Stress




#!!!!StressInit 是满矩阵，加上StressAdd里没用的去掉后，是否可以通过np.nan+stress0=na.nan来解决那个问题？试了貌似可以
def StressList(Stress,StresTrasFunc,BrokenBondLabel): #Calculate it each loop of each MC# Delete the broken bonds' stress

    StresTrasFunc_Brok=StresTrasFunc[BrokenBondLabel]
    StresTrasFunc_Brok_sum = np.nansum(StresTrasFunc_Brok)#不需要全计算出来，计算一行和即可
    StresTrasFunc_Brok=StresTrasFunc_Brok/StresTrasFunc_Brok_sum

    StressAdd=StresTrasFunc_Brok*Stress[BrokenBondLabel]
    StresTrasFunc[:,BrokenBondLabel]= np.nan # set just failed bond's StresTrasFunc to be nan
    #print("Stress is:", StressAdd[0:120])
    #print("")
    Stress+=StressAdd
    Stress[BrokenBondLabel]=np.nan
    #print(Stress)
    #print("")
    return (Stress,StresTrasFunc) #Delete StressAdd later on

#Stress=StressList(Stress,F,3)
# 74644478.58274122
#rate0=BreakingRateVsStress(k, sigma0, T, Ea, nu, tau0)

# def Sigma0ListInit(sigma0):
#     Sigma0=np.full(2*N, sigma0)
#     for i in BondNone:
#         Sigma0[i] = 0
#     return Sigma0
# Sigma0=Sigma0ListInit(sigma0)
#
# def Sigma0Lis(Sigma0,BrokenBondLabel):
#     Sigma0[BrokenBondLabel]=0
#     return Sigma0
# Sigma0=Sigma0Lis(Sigma0,3)

def RateList(Stress,nu,sigma0,rate0,T): #Calculate it each loop of each MC
    kb = 1.3806 * 10 ** (-23)
    beta = (kb * T) ** -1
    Rate=rate0*np.exp(beta*nu*(Stress-sigma0)) #!!!!!系数弄出去
    #print(beta*nu*(Stress-sigma),nu,Stress,"sig",sigma)
    # Rate=rate0*np.exp(nu*(Stress-Sigma0)) #!!!!!系数弄出去
    #Rate=np.nan_to_num(Rate)
    return Rate

#Rate=RateList(Stress,nu,sigma0,rate0)
#print(Rate)
#print(Stress)

def main():
    BondNone = BondNoneList(L)
    Distance = DistanceMatrix(L, N, BondNone)
    F = StressTransferFunction(100, Distance)
    Stress=StressInit(sigma0,N)
    rate0=BreakingRateVsStress(k, sigma0, T, Ea, nu, tau0)
    for i in [12, 0, 18, 14, 15, 11, 9, 10, 2, 16, 17, 5, 4, 1, 6, 7, 13, 8, 3]:
        Stress=StressList(Stress,F,i)
        Rate=RateList(Stress,nu,sigma0,rate0,T)
        #print(Rate[0:200])
        #print(Stress)
def StressRedisAllInit(L,N,sigma0,T,Ea,nu,tau0,k,gamma,ratio=1):
    BondNone = BondNoneList(L)
    # check whether IO has that Distance


    if not os.path.exists("./Init/Distance/SL_L_{}.pickle".format(L)):
        Distance = DistanceMatrix(L, N, BondNone,ratio)
        # 存起来
        if not os.path.exists('./Init/Distance/'.format(L)):
            os.makedirs('./Init/Distance/'.format(L))
        try:
            with open('./Init/Distance/SL_L_{}.pickle'.format(L),'wb') as F:
                pickle.dump(Distance,F)
        except OverflowError:
            print("OverflowError. I am not saving it. And I will delete the blank pickle file.")
            os.remove('./Init/Distance/SL_L_{}.pickle'.format(L))
            
    else:
        with open('./Init/Distance/SL_L_{}.pickle'.format(L), 'rb') as F:
            Distance = pickle.load(F)
        print("Load!!!!")

    #F matrix就不管了 因为毕竟每个gamma几乎只运行一次，但也不一定。可以管一下，怕占太多storage
    #要是非想存它，只要模仿Distance即可
    F = StressTransferFunction(gamma, Distance)
    Stress=StressInit(sigma0,N,BondNone)
    rate0=BreakingRateVsStress(k, sigma0, T, Ea, nu, tau0)

    Rate=RateList(Stress,nu,sigma0,rate0,T)
    return (Stress,Rate,F,rate0)

def TrasFuncInit(L,gamma,ratio=1): #没用 StressRedisAllInit已经有了
    BondNone = BondNoneList(L)
    Distance = DistanceMatrix(L, N, BondNone,ratio)
    F = StressTransferFunction(gamma, Distance)
    return F

def LatticeInit(L):
    BondNone = BondNoneList(L)
    Distance = DistanceMatrix(L, L*L, BondNone)
    return Distance

def StressRedisAllStep(Stress,F,BrokenBondLabel,rate0,nu,sigma0,T):
    (Stress,F)=StressList(Stress,F,BrokenBondLabel)
    Rate=RateList(Stress,nu,sigma0,rate0,T)
    return (Rate,Stress,F)


if __name__ == "__main__":
    # begin=time.time()
    # main()
    # now=time.time()
    # print(now-begin)
    L = 50
    N = L ** 2
    N_e = 2 * N - 2 * L
    k = 1
    T = 298
    Ea = 117.15 * 10 ** 3 / (6.02214 * 10 ** 23)
    nu = 418 * 10 ** (-30)
    tau0 = 10 ** (-11)
    sigma0=8E6
    rate0=BreakingRateVsStress(k, sigma0, T, Ea, nu, tau0)
    print(rate0)
    (Stress,F,rate0,Rate)=StressRedisAllInit()
    for BrokenBondLabel in [2, 50, 57, 48, 35, 64, 59, 0, 78, 60, 8, 7, 58, 67, 1, 42, 76, 22, 9, 71, 36, 68, 53, 18, 37, 49, 77, 33, 3, 12, 26, 47, 15, 29, 56, 11, 46, 14, 52, 24, 51, 45, 70, 75, 38, 13, 21, 5, 31, 34, 27, 4, 63, 73, 10, 72, 69, 66, 19, 16, 65, 20, 54, 32, 62, 74, 55, 61, 23, 41, 17, 40, 6, 43, 28, 25, 44, 30]:
        #Stress=StressList(Stress,F,BrokenBondLabel)
        #Rate=RateList(Stress,nu,sigma0,rate0,T)
        Rate,Stress=StressRedisAllStep(Stress,F,BrokenBondLabel,rate0)
        print(Rate[0:120])

#0.043695926666259766*(10000/20) every MC it takes 21.847963333129883 s for L=100  !!!!不对吧？？ L=20还需要1.01秒呢？？？
#But it is so fucking slow for the initialization (15 mins), so we could save it somewhere? But it is bearable.

1.69350878e-05
########问题1： 已经断裂的bond为什么还给他们加stress？
########问题2： 一个bond断裂，没有把
########大问题，你傻逼啊，square lattice nearest neighbour不是6个，只有四个，另外两个是次近邻，和隔一个（一对）一样距离
