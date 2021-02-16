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




def DistanceMatrix(L, N, BondNone, ratio=1):
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

    return F






def StressInit(stress, N,BondNone):
    Stress = np.full(2 * N, stress)

    for i in BondNone:
        Stress[i] = np.nan
    return Stress



def StressList(Stress,StresTrasFunc,BrokenBondLabel): #Calculate it each loop of each MC# Delete the broken bonds' stress

    StresTrasFunc_Brok=StresTrasFunc[BrokenBondLabel]
    StresTrasFunc_Brok_sum = np.nansum(StresTrasFunc_Brok)
    StresTrasFunc_Brok=StresTrasFunc_Brok/StresTrasFunc_Brok_sum

    StressAdd=StresTrasFunc_Brok*Stress[BrokenBondLabel]
    StresTrasFunc[:,BrokenBondLabel]= np.nan
    Stress+=StressAdd
    Stress[BrokenBondLabel]=np.nan
    return (Stress,StresTrasFunc) #Delete StressAdd later on


def RateList(Stress,nu,sigma0,rate0,T): #Calculate it each loop of each MC
    kb = 1.3806 * 10 ** (-23)
    beta = (kb * T) ** -1
    Rate=rate0*np.exp(beta*nu*(Stress-sigma0))
    return Rate


def main():
    BondNone = BondNoneList(L)
    Distance = DistanceMatrix(L, N, BondNone)
    F = StressTransferFunction(100, Distance)
    Stress=StressInit(sigma0,N)
    rate0=BreakingRateVsStress(k, sigma0, T, Ea, nu, tau0)
    for i in [12, 0, 18, 14, 15, 11, 9, 10, 2, 16, 17, 5, 4, 1, 6, 7, 13, 8, 3]:
        Stress=StressList(Stress,F,i)
        Rate=RateList(Stress,nu,sigma0,rate0,T)

def StressRedisAllInit(L,N,sigma0,T,Ea,nu,tau0,k,gamma,ratio=1):
    BondNone = BondNoneList(L)
    # check whether IO has that Distance


    if not os.path.exists("./Init/Distance/SL_L_{}.pickle".format(L)):
        Distance = DistanceMatrix(L, N, BondNone,ratio)

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


    F = StressTransferFunction(gamma, Distance)
    Stress=StressInit(sigma0,N,BondNone)
    rate0=BreakingRateVsStress(k, sigma0, T, Ea, nu, tau0)

    Rate=RateList(Stress,nu,sigma0,rate0,T)
    return (Stress,Rate,F,rate0)

def TrasFuncInit(L,gamma,ratio=1):
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
        
        Rate,Stress = StressRedisAllStep(Stress,F,BrokenBondLabel,rate0)



