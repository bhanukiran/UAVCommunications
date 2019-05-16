import scipy as sp
from Parameters import *
from System import *
import pylab as plab
import matplotlib as plt
import time as t
from math import ceil

sp.random.seed(17081986)

p = 5 #grid partition

def CostFunctionPt(Arr_b, Arr_u, Pt_max):
    X_list_b = Arr_b[0]
    Y_list_b = Arr_b[1]
    Z_list_b = Arr_b[2]
    
    #Nb BSs initially taken - uniformly random
    BS_list = [BaseStation(p[0], p[1], p[2]) for p in zip(X_list_b, Y_list_b, Z_list_b)]

    X_list_u = Arr_u[0]
    Y_list_u = Arr_u[1]

    #Nu users uniformly distributed
    US_list = [User(p[0], p[1]) for p in zip(X_list_u, Y_list_u) ]
    
    S = System(US_list, BS_list)
    # print("opt = ", S.optimal_power())
    return S.optimal_power(Pt_max)

#----------------------------------------------

def PSO_1(pop, USlist, gbest, MaxIter, c1, c2, w, wdamp, Pt_max, K):
    """
    c1: personal learning coefficient
    c2: global learning coefficient
    w: weight -> speed of convergance
    wdamp: 
    """

    X_list_u = USlist[0]
    Y_list_u = USlist[1]
    
    # PSO Loop
    PopSize = len(pop)
    sp.random.seed()
    best_array = []
    best_cost_array = []
    for it in range(0, MaxIter):
        for i in range(0, PopSize):
            
            # update velocity matrix for each particle
            pop[i]['velocity'] = w*pop[i]['velocity'] + c1*sp.random.rand(3, Nu//K)*(pop[i]['best_position'] - pop[i]['position']) + c2*sp.random.rand(3, Nu//K)*(gbest['position'] - pop[i]['position'])

            # update positions
            pop[i]['position'] += pop[i]['velocity'];

            # if x,y,z coordinates cross the bounds replace them with the corresponding maximum/minimun bound
            pop[i]['position'][0] = sp.maximum(pop[i]['position'][0], -Lx/2);
            pop[i]['position'][0] = sp.minimum(pop[i]['position'][0], Lx/2);

            pop[i]['position'][1] = sp.maximum(pop[i]['position'][1], -Ly/2);
            pop[i]['position'][1] = sp.minimum(pop[i]['position'][1], Ly/2);

            pop[i]['position'][2] = sp.maximum(pop[i]['position'][2], Lz_min);
            pop[i]['position'][2] = sp.minimum(pop[i]['position'][2], Lz);

            # now calculate the utility function
            pop[i]['cost'] = CostFunctionPt(pop[i]['position'], [X_list_u, Y_list_u], Pt_max);
            
            if pop[i]['cost'] < pop[i]['best_cost']:
                pop[i]['best_position'] = pop[i]['position'].copy();
                pop[i]['best_cost'] = pop[i]['cost'];

                if pop[i]['best_cost'] < gbest['cost']:
                    gbest['position'] = pop[i]['best_position'].copy();
                    gbest['cost'] = pop[i]['best_cost'];

        w *= wdamp;
        
        best_array.append(gbest['position'])
        best_cost_array.append(gbest['cost'])

        print('Iteration {}: Best Cost = {}'.format(it+1, gbest['cost']));

    return pop, best_cost_array[-1], best_array;

                                
if __name__ == "__main__":

    best_array, best_cost_array, pop, X_list_u, Y_list_u= PSO(10, 0, 1.3962, 1.3962, 0.3298, 1.0, 33)

    # with open('BS_array_%d.txt'%PSO.stamp, 'w') as thefile:
    #     thefile.write(str(best_array) + '\n')

    # with open("best_cost_array_%s.txt"%PSO.stamp, 'w') as best_cost_file:
    #     best_cost_file.write(str(best_cost_array))
    

