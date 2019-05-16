import scipy as sp
from Parameters import *
from System import *
import pylab as plab
import matplotlib as plt
import time as t
from math import ceil
from Simulation1 import CostFunctionPt, PSO_1
from Simulation2 import CostFunctionK, PSO_2

sp.random.seed(17081986)
spread = 1800

def PtKinitial(trials):
    
    Pt_arr = sp.linspace(0.01, 2, 20)
    K_arr  = sp.linspace(10, 50, 20, dtype='int')

    X_list_u = sp.random.normal(0, spread, Nu)
    Y_list_u = sp.random.normal(0, spread, Nu)
                
    US_list = [User(p[0], p[1]) for p in zip(X_list_u, Y_list_u)]

    Z = sp.zeros((len(K_arr),len(Pt_arr)))
    
    for t in range(trials):
        PtK = []
        for K in K_arr[::-1]:
            Pts = []
            for Pt in Pt_arr:
                X_list_b = sp.random.uniform(-Lx/2, Lx/2, Nu//K)
                Y_list_b = sp.random.uniform(-Ly/2, Ly/2, Nu//K)
                Z_list_b = sp.random.uniform(Lz_min, Lz, Nu//K)
                
                BS_list = [BaseStation(p[0], p[1], p[2]) for p in zip(X_list_b, Y_list_b, Z_list_b)]
                
                S = System(US_list, BS_list)
                
                Pts.append(S.avg_cov_prob(Pt) >= xi and S.avg_rate(Pt) >= beta)
                
            PtK.append(Pts)
            
        Z += sp.matrix(PtK)
    Z = sp.matrix(Z)/trials

    f = lambda x: int(x > 0.5); f_vec = sp.vectorize(f)
    Z = f_vec(Z)

    print(Z)
    
    non_zero = sp.nonzero(Z)
    nr = list(non_zero[0])
    nc = list(non_zero[1])
    Kmin = K_arr[::-1][nr[nc.index(min(nc))]]
    Ptmax = Pt_arr[nc[nr.index(min(nr))]]

    print("Ptmax = ",Ptmax,"K = ",Kmin)
    return Ptmax, Kmin

#----------------------------------------------

def BackForth(PopSize, NumIters, PSOIter, c1, c2, w, wdamp, Pt_max, K):

    #NumIters - number of back forths
    #PSOIter - number of PSO iterations

    X_list_u = sp.random.normal(0, spread, Nu)
    Y_list_u = sp.random.normal(0, spread, Nu)

    USlist = [X_list_u, Y_list_u]

    with open("user_array.txt", 'w') as user_file:
        user_file.write(str(USlist))
    
    # Empty Particle Template
    empty_particle = {
        'position': None,
        'velocity': None,
        'cost': None,
        'best_position': None,
        'best_cost': None,
    }
    
    # Initialize Global Best
    gbest = {'position': None, 'cost': sp.inf}
    
    # Create Initial Population
    pop = []

    for i in range(0, PopSize):
        pop.append(empty_particle.copy())
        #------------------------------
        X_list_b = sp.random.uniform(-Lx/2, Lx/2, Nu//K)
        Y_list_b = sp.random.uniform(-Ly/2, Ly/2, Nu//K)
                
        Z_list_b = sp.random.uniform(Lz_min, Lz, Nu//K)
        #------------------------------
        pop[i]['position'] = sp.array([ X_list_b, Y_list_b, Z_list_b])
        pop[i]['velocity'] = sp.zeros((3,Nu//K))
        
        pop[i]['cost'] = CostFunctionPt(pop[i]['position'], [X_list_u, Y_list_u], Pt_max)
        pop[i]['best_position'] = pop[i]['position'].copy()
        pop[i]['best_cost'] = pop[i]['cost']
        
        print("Cost of ",i+1,"th particle = ",pop[i]['cost'])
        
        if pop[i]['best_cost'] < gbest['cost']:
            gbest['position'] = pop[i]['best_position'].copy()
            gbest['cost'] = pop[i]['best_cost']

    initial_best = gbest['position']
    initial_best_cost = gbest['cost']
    print("gbest_position = ", gbest['position'])
    
    if gbest['cost'] == sp.inf or PSOIter == 0:
        best_array = [0]
        best_cost_array = [0]
        return best_cost_array, best_array

    cost_array = []
    BSlist = []
    best_array = [gbest['position']]
    
    for it in range(NumIters):
        # flag specifies the utility function
        flag = 1
        check = 0    
        
        #CALLING PSO1!
        pop, Pt_max, bestpart= PSO_1(pop, USlist, gbest, PSOIter, c1, c2, w, wdamp, Pt_max, K)
        cost_array.append(Pt_max)
        best_array.append(bestpart)
        
        gbest = {'position': None, 'cost': sp.NINF}
        
        for i in range(PopSize):
            pop[i]['velocity'] = sp.zeros((3,Nu//K))
            
            pop[i]['cost'] = CostFunctionK(pop[i]['position'], [X_list_u, Y_list_u], flag, check, Pt_max)
            pop[i]['best_position'] = pop[i]['position'].copy()
            pop[i]['best_cost'] = pop[i]['cost']
            
            print("Cost of",i+1,"th particle (Utility",flag,") =" ,pop[i]['cost'])
        
            if pop[i]['best_cost'] > gbest['cost']:
                gbest['position'] = pop[i]['best_position'].copy()
                gbest['cost'] = pop[i]['best_cost']

        initial_best = gbest['position']
        initial_best_cost = gbest['cost']
        print("gbest_position = ", gbest['position'])
        
        if gbest['cost'] == sp.NINF or PSOIter == 0:
            return cost_array, best_array

        #CALLING PSO2!    
        BSlist, bestpart = PSO_2(pop, USlist, gbest, PSOIter, c1, c2, w, wdamp, Pt_max, K)
        best_array.append(bestpart)
        cost_array.append(len(BSlist))
        K = Nu//len(BSlist)
        
        gbest = {'position': None, 'cost': sp.inf}
        for i in range(PopSize):
            
            #------------------------------
            X_list_b = sp.random.uniform(-Lx/2, Lx/2, Nu//K)
            Y_list_b = sp.random.uniform(-Ly/2, Ly/2, Nu//K)
            
            Z_list_b = sp.random.uniform(Lz_min, Lz, Nu//K)
            #------------------------------
            pop[i]['position'] = sp.array([ X_list_b, Y_list_b, Z_list_b])
            
            if i==0:
                pop[i]['position'] = sp.array([[b.x for b in BSlist], [b.y for b in BSlist], [b.z for b in BSlist]])
                print("BSlist = ",pop[i]['position'])
                
            pop[i]['velocity'] = sp.zeros((3,Nu//K))
            
            pop[i]['cost'] = CostFunctionPt(pop[i]['position'], [X_list_u, Y_list_u], Pt_max)
            pop[i]['best_position'] = pop[i]['position'].copy()
            pop[i]['best_cost'] = pop[i]['cost']
            
            print("Cost of ",i+1,"th particle = ",pop[i]['cost'])
            
            if pop[i]['best_cost'] < gbest['cost']:
                gbest['position'] = pop[i]['best_position'].copy()
                gbest['cost'] = pop[i]['best_cost']
                
        initial_best = gbest['position']
        initial_best_cost = gbest['cost']
        print("gbest_position = ", gbest['position'])
        
        if gbest['cost'] == sp.inf or PSOIter == 0:
            return cost_array, best_array
                
    print("Pt_max = ",Pt_max,", K = ", K)
        
    return cost_array, best_array
    
if __name__ == "__main__":

    # Pt_max, K = PtKinitial(30)
    Pt_max = 2.5; K = 1000//28
    cost_array, config_array= BackForth(50, 20, 30, 1.3962, 1.3962, 0.3298, 1.0, Pt_max, K)

    # best_array, best_cost_array, pop, X_list_u, Y_list_u= PSO(10, 0, 1.3962, 1.3962, 0.3298, 1.0, 33)

    with open('BS_array.txt', 'w') as thefile:
        thefile.write(str(config_array) + '\n')

    with open("best_cost_array.txt", 'w') as best_cost_file:
        best_cost_file.write(str(cost_array))
    

