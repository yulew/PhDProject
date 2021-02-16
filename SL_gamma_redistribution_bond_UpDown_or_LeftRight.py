__author__ = 'Yule'

import os
import random
import numpy as np
import copy
np.set_printoptions(threshold=np.nan)

from time import gmtime, strftime

import time
from math import log

import networkx as nx
import matplotlib.pylab as plt
import pickle

from gamma_StressRedisAllLatt import StressRedisAllInit, StressRedisAllStep, RateList


from GaussianGraphSave import GaussianGraphSave

import sys
sys.setrecursionlimit(5000)
def BreakingRateVsStress(k, sigma, T, Ea, nu, tau0):
    e0 = 0
    kb = 1.3806 * 10 ** (-23)
    beta = (kb * T) ** -1
    for j in range(1, k + 1):
        e0 += np.exp(-beta * nu * sigma * k / j) / j
    t_b = e0 * tau0 * np.exp(beta * Ea) / 3600  # hours
    rate_b = (t_b) ** (-1)
    return rate_b


def binarySearch(alist, item):
    first = 0
    last = len(alist) - 1
    found = False
    number = 0

    while first <= last and not found:
        midpoint = (first + last) // 2


        if item > alist[midpoint - 1] and item <= alist[midpoint]:
            found = True
            number = midpoint  # label No.

        else:
            if item < alist[midpoint]:
                last = midpoint - 1

            else:
                first = midpoint + 1

    return number


def Initialization_SquareLattice_Bond(L, N):
    def BondNumberLabel(L):

        label = []
        for i in range(0, L - 1):
            for j in range(0, L - 1):
                label.append([2 * j + 2 * i * L, 2 * j + 1 + 2 * i * L])
            label.append([2 * i * L + 2 * (L - 1), None])
        for j in range(L * L - L, L * L - 1):  # changed from L*L to L*L-1
            label.append([None, 2 * j + 1])  # in total (2*L-1) np.nan
        label.append([None, None])
        return label

    def NearestNeighbourCorrelationInitiliaztion(N):

        nn = []
        for i in range(N - 1):
            if label[i][0] == None:
                nn.append([None, [0, 0, 0, 0, 0, 0]])

            elif label[i][1] == None:
                nn.append([[0, 0, 0, 0, 0, 0], None])

            else:
                nn.append([[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])
        nn.append([None, None])
        return nn

    def NearestNeighbourCorrelation(L, N, label, nn):
        for i in range(N - 1):  # check the boundary whether it exceeded
            if label[i][0] == None:
                nn[i][1][0] = label[i][0]
                nn[i][1][1] = label[(i - L + N) % N][0]  # up cell
                nn[i][1][2] = label[(i - L + 1 + N) % N][0]  # upright
                nn[i][1][3] = label[(i - 1 + N) % N][1]  # left
                nn[i][1][4] = label[(i + 1 + N) % N][0]  # right
                nn[i][1][5] = label[(i + 1 + N) % N][1]  # right
            elif label[i][1] == None:
                nn[i][0][0] = label[i][1]
                nn[i][0][1] = label[(i - L + N) % N][0]  # up
                nn[i][0][2] = label[(i - 1 + N) % N][1]  # left
                nn[i][0][3] = label[(i + L - 1 + N) % N][1]  # downleft
                nn[i][0][4] = label[(i + L + N) % N][0]  # down
                nn[i][0][5] = label[(i + L + N) % N][1]  # down
            else:
                nn[i][0][0] = label[i][1]
                nn[i][0][1] = label[(i - L + N) % N][0]  # up
                nn[i][0][2] = label[(i - 1 + N) % N][1]  # left
                nn[i][0][3] = label[(i + L - 1 + N) % N][1]  # downleft
                nn[i][0][4] = label[(i + L + N) % N][0]  # down
                nn[i][0][5] = label[(i + L + N) % N][1]  # down

                nn[i][1][0] = label[i][0]
                nn[i][1][1] = label[(i - L + N) % N][0]  # up
                nn[i][1][2] = label[(i - L + 1 + N) % N][0]  # upright
                nn[i][1][3] = label[(i - 1 + N) % N][1]  # left
                nn[i][1][4] = label[(i + 1 + N) % N][0]  # right
                nn[i][1][5] = label[(i + 1 + N) % N][1]  # right

        for i in range(L):  # upper boundary
            if nn[i][1] != None:
                for j in [1, 2, 3, 5]:
                    nn[i][1][j] = None
        for i in range(N - L, N):  # bottom boundary
            if nn[i][1] != None:
                for j in [0, 3, 4, 5]:
                    nn[i][1][j] = None
        for i in range(0, N, L):  # left boundary
            if nn[i][0] != None:
                for j in [1, 2, 3, 4]:
                    nn[i][0][j] = None
        for i in range(L - 1, N, L):  # right boundary
            if nn[i][0] != None:
                for j in [0, 1, 4, 5]:
                    nn[i][0][j] = None
        for i in range(L):  # upper boundary but a little lower (label[i][0])
            if nn[i][0] != None:
                nn[i][0][1] = None

        for i in range(L):  # left boundary but a little lower (label[i][1])
            if nn[i][1] != None:
                nn[i][1][3] = None
        return nn

    label = BondNumberLabel(L)
    nn = NearestNeighbourCorrelationInitiliaztion(N)
    nn = NearestNeighbourCorrelation(L, N, label, nn)
    return (label, nn)




class RandomPercolation_SquareLattice_Bond:
    def __init__(self, label, lifetime=False, rate=None):
        # If lifetime == True, also will calculate lifetime.
        # Also, if lifetime is True, rate should be entered, otherwise raise Error.
        self.order = []
        self.lifetime = lifetime
        self.rate = rate
        self.n = N_e  # for ElementFailNumber_tau(self) 计数，算t_b（还剩几个元素在 order_per里）

        for i in range(N):  # N-1, not N, because the last one all [None,None] # 你忘了还要append label[i][1]了

            self.order.append(label[i][0])

            self.order.append(label[i][1])

        if self.lifetime:  # correlated with self.ElementFailNumber_tau()
            self.order_gene = iter(self.order_permutation())

    def order_permutation(self):
        order_per = self.order[:]  # when change order_per, it doese not change the original list
        for i in range(len(order_per)):
            j = i + random.randint(0, len(order_per) - i - 1)
            temp = order_per[i]
            order_per[i] = order_per[j]
            order_per[j] = temp
        return order_per

    def ElementFailNumber_tau(self):
        # calculate lifetime
        if self.lifetime:
            try:
                rho1 = random.random()
                t_b = -log(rho1) / (self.n * self.rate)
                self.n -= 1
            except:
                raise TypeError("Please enter the rate in a correct form.")
        # calculate element_fail_number
        element_fail_number = next(self.order_gene)
        return (element_fail_number, t_b)


class CorrelatedPercolation_SquareLattice_Bond:
    # list rate0 contains a bunch of different rates; list percentage0 contains a bunch of corresponding percents
    # If homo, only rate0 is picked up. And no need to input percent; if inhomo, rate0 becomes a list,
    # pick up all and list percentage0 from ** karg.
    def __init__(self,  homo=True, StressRedistribution=True, Stress0=None,Rate0=None,F0=None,rate0=None, **kwargs):  # rate0通常通过BreakingRateVsStress算出来
        self.StressRedistribution = StressRedistribution
        if self.StressRedistribution:


            (self.Stress,self.Rate,self.F,self.rate0)=(np.copy(Stress0),np.copy(Rate0),np.copy(F0),np.copy(rate0))


        self.index=0
        self.trigger=False
        self.FakeCancel = 0


    def kineticMC(self,nu,sigma0):

        if self.trigger == False:  # If stress is too large, rates may overflow. Therefore, we need to check whether it the stress is too large.
            sum_till_list = np.nancumsum(self.Rate)
            #print("Sum list:",sum_till_list[:40])
            rho1 = random.random()
            rho2 = random.random()

            self.t_b = -log(rho1) / sum_till_list[-1]
            ran_tot_rate = rho2 * sum_till_list[-1]

            self.element_fail_number = binarySearch(sum_till_list, ran_tot_rate)
        else:  # If stress is too large, we can cancell out stress to generate a new fake stress.
            sum_till_list = np.nancumsum(self.Rate)
            rho2 = random.random()
            ran_tot_rate = rho2 * sum_till_list[-1]
            self.element_fail_number = binarySearch(sum_till_list, ran_tot_rate)
            self.t_b = 0


        if self.StressRedistribution:
            self.StressRedis(self.element_fail_number,nu,sigma0)

        return (self.element_fail_number, self.t_b)

    def StressRedis(self, element_fail_number,nu,sigma0):
        self.Stress += self.FakeCancel
        (self.Rate,self.Stress,self.F)=StressRedisAllStep(self.Stress,self.F,element_fail_number,self.rate0,nu,sigma0,T) # Align up Stress with sigma, and align up Rate ith rate_list

        # if one of Stress matrix exceeds the threshold, then everyone cancell out a number
        if self.step > 80:
            self.FakeCancel=np.nanmax(self.Stress)-6*sigma0

            if self.FakeCancel>4E9:
                self.trigger = True
                self.Stress = self.Stress-self.FakeCancel #the new stress is acutally a fake stress, but we do not care about it for a moment
                self.Rate=RateList(self.Stress,nu,sigma0,self.rate0,T)

                print("triggered")
            else:
                self.trigger = False
                self.FakeCancel=0


        return self.Rate
    def ElementFailNumber_tau(self,nu,sigma0,step):
        self.step = step
        self.kineticMC(nu,sigma0)

        return (self.element_fail_number, self.t_b)


def percolate(RandomPerco=False):

    EMPTY = -2 * N - 1
    EMPTY1 = -2 * N - 2


    big = 0
    breaking_time = 0

    m1 = EMPTY1
    m2 = EMPTY1
    m3 = EMPTY1
    m4 = EMPTY1

    rate0 = BreakingRateVsStress(k, sigma0, T, Ea, nu, tau0)

    Order_Time=[]


    orderOnly=[]

    if RandomPerco:
        PercoInit = RandomPercolation_SquareLattice_Bond(label, lifetime=False, rate=False)
    else:  # correlated
        PercoInit = CorrelatedPercolation_SquareLattice_Bond( homo=True, StressRedistribution=True,Stress0=Stress0,Rate0=Rate0,F0=F0,rate0=rate0)
    for i in range(2 * N):
        ptr1[i] = EMPTY
        ptr2[i] = EMPTY

    #####
    LatticeInit = SquareLatticeGraphC()
    for i in range(N_e):
        ###别忘修改！！！！！！去掉testStress
        s1, tau,= PercoInit.ElementFailNumber_tau(nu,sigma0,i)

        breaking_time += tau

        # ###去掉下面的
        # testRateList.append(testRate)
        # #####

        LatticeInit.bondAdd(s1)

        if s1 in range(1, 2 * L, 2):
            if m1 == EMPTY1:
                ptr1[s1] = -1
                m1 = s1  # *

            else:
                ptr1[s1] = m1

        elif s1 in range(2 * N - 2 * L + 1 , 2 * N, 2): #这个是检查下boundar

            if m2 == EMPTY1:
                ptr1[s1] = -1
                m2 = s1  # *

            else:
                ptr1[s1] = m2

        else:
            ptr1[s1] = -1


        if s1 in range(0, 2 * N, 2 * L):#left boundary
            if m3 == EMPTY1:
                ptr2[s1] = -1
                m3 = s1  # *

            else:
                ptr2[s1] = m3
        elif s1 in range(2 * L - 2, 2 * N, 2 * L):#right boundary
            if m4 == EMPTY1:
                ptr2[s1] = -1
                m4 = s1  # *

            else:
                ptr2[s1] = m4

        else:
            ptr2[s1] = -1





        r11 = findroot(ptr1,s1)
        r21 = findroot(ptr2,s1)




        #print("break:",s1," ")
        for j in range(6):   #
            # print(nn[s1//2][s1%2])

            s2 = nn[s1 // 2][s1 % 2][j]
            if s2 != None:
                if ptr1[s2] != EMPTY:
                    r12 = findroot(ptr1,s2)
                    if r12 != r11:

                        if ptr1[r11] > ptr1[r12]:
                            ptr1[r12] += ptr1[r11]
                            ptr1[r11] = r12
                            r11 = r12
                        else:
                            ptr1[r11] += ptr1[r12]
                            ptr1[r12] = r11
                        if (-ptr1[r11] > big):
                            big = -ptr1[r11]


                if ptr2[s2] != EMPTY:
                    r22 = findroot(ptr2,s2)
                    if r22 != r21:

                        if ptr2[r21] > ptr2[r22]:
                            ptr2[r22] += ptr2[r21]
                            ptr2[r21] = r22
                            r21 = r22
                        else:
                            ptr2[r21] += ptr2[r22]
                            ptr2[r22] = r21
                        if (-ptr2[r21] > big):
                            big = -ptr2[r21]

        Order_Time.append((s1,breaking_time))


        orderOnly.append(s1)



        if yes_no(m1,m2,m3,m4,EMPTY1):

            break

    threshold=(i + 1) / ( 2 * N - 2 * L)


    LatticeInit.LastBond(s1)
    LatticeInit.drawLattice(True,threshold,orderOnly)
    return (threshold,breaking_time,Order_Time)



def yes_no(m1, m2, m3, m4, EMPTY1):

    tr = False
    tr1 = False
    tr2 = False

    if m1 != EMPTY1 and m2 != EMPTY1:
        f1 = findroot(ptr1,m1)
        f2 = findroot(ptr1,m2)
        if f1 == f2 :
            tr1 = True

    if m3 != EMPTY1 and m4 != EMPTY1:
        f3 = findroot(ptr2,m3)
        f4 = findroot(ptr2,m4)
        if f3 == f4 :
            tr2 = True

    if tr1 or tr2:
        tr = True

    return tr



def findroot(ptr,i):

    if ptr[i] < 0:
        return i
    ptr[i] = findroot(ptr,ptr[i])
    return ptr[i]


class SquareLatticeGraphC:
    def __init__(self):
        self.dic = {}
        for x in range(L * L):
            self.dic[label[x][0]] = (x, x + L)
            self.dic[label[x][1]] = (x, x + 1)
        self.Square = nx.Graph()
        for i in range(L):
            for j in range(L):
                self.Square.add_node(L * j + i, pos=(i, -j), size=0.00000001)
        self.pos = nx.get_node_attributes(self.Square, 'pos')

    def bondAdd(self, s1):
        self.Square.add_edge(self.dic[s1][0], self.dic[s1][1],color='black',weight=1)

    def LastBond(self,s1):
        self.Square.add_edge(self.dic[s1][0], self.dic[s1][1],color='g',weight=4)

    def drawLattice(self, Draw, threshold="",orderOnly=None):  # if percolate, draw lattice, Draw=True or False
        if Draw:
            edges = self.Square.edges()
            weights = [self.Square[u][v]['weight'] for u,v in edges]
            colors = [self.Square[u][v]['color'] for u,v in edges]
            nx.draw(self.Square, self.pos,edges=edges, edge_color=colors,width=weights, node_size=3)
        #plt.show()
        #bond_index=[(x // (2 * L), (x % (2*L))//2, x%2) for x in orderOnly[0:200]]
        Filename=("L%s_stress%s_gamma%s_perco_%s_EitherSide_" %(L,sigma0,gamma,threshold))
        tim=strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time()))
        #dirc=("/Users/Yule/Documents/presentation/Committee_Meeting/png_animation/L%s_gamma%s_MC%s/" %(L,gamma,MC))
        dir="/Users/Yule/Documents/presentation/Committee_Meeting/"
        if not os.path.exists(dir):
            os.makedirs(dir)
        plt.savefig(dir+Filename+tim+".png")
        plt.clf()



start = time.time()
if __name__ == "__main__":
    def save_file():
        Filename=("stress_%s_gamma_%s_L_%s_" %(sigma0,gamma,L))
        Dir=("/Users/Yule/Documents/PYTHON/Fracture_stress_redistribution/RedistriAll_gamma/data_bond_SL/"+Filename+strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time()))+".txt")
        File = open(Dir,'w')
        File.write("%s MC runs.\n\n" % str(MC+1))
        File.write("Percolation threshold:\n\n" )
        File.write(str(perco_proportion)+"\n\n")
        File.write("Break time:\n\n")
        File.write(str(break_time)+"\n\n")
        File.close()

        filename = ("OrderTime_stress_%s_gamma_%s_L_%s_" %(sigma0,gamma,L))
        dir=("/Users/Yule/Documents/PYTHON/Fracture_stress_redistribution/RedistriAll_gamma/data_bond_SL/"+filename+strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time())))
        with open(dir, 'wb') as file_order:
            pickle.dump(Order_TimeList, file_order, protocol=pickle.HIGHEST_PROTOCOL)


    def save_file(Order_TimeList,gamma=None,sigma0=None,ratio=None,SiteOrBond=None,EihterSide=None,lattice=None,L=None,MC=None):
        filename = '%s_L_%s_sigma0_%s_gamma_%s_ratio_%s_times_%s.pickle' %(lattice,L,sigma0,gamma,ratio,(MC+1))
        if not os.path.exists('./data_order_time/%s_%s_%s/' %(SiteOrBond,lattice,EihterSide)):
            os.makedirs('./data_order_time/%s_%s_%s/' %(SiteOrBond,lattice,EihterSide))
        with open(('./data_order_time/%s_%s_%s/' %(SiteOrBond,lattice,EihterSide))+filename, 'wb') as file:
            pickle.dump(Order_TimeList, file, protocol = pickle.HIGHEST_PROTOCOL)



end = time.time()
print("Time elapse:", (end-start))


