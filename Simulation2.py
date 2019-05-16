import scipy as sp
from Parameters import *
from System import *
import pylab as p
import matplotlib as plt
import time as t

sp.random.seed(17081986)

def CostFunctionK(Arr_b, Arr_u, flag, check, Pt_max):
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

    if check == 1:
        if S.avg_cov_prob(Pt_max) >= xi:
            return True
        else:
            return False

    if flag == 1:
        return S.avg_cov_prob(Pt_max)
    else:
        return S.avg_rate(Pt_max)

#----------------------------------------------
def users_per_cell(X_list_u, Y_list_u):
 #   L = [[0,0,0,0] for i in range(4)]
    L = [0 for i in range(p**2)]

    c = 0
    while c < Nu:
        
        x = X_list_u[c]
        y = Y_list_u[c]

        i = int((5000+x)/(Lx/p))
        j = int((5000+y)/(Ly/p))
        
        try:
            L[p*j + i] += 1
        except IndexError:
            pass
        
        c += 1

    # print("L = ", L)
    L = sorted([i for i in zip(L, range(len(L)))])[::-1]
    return L

def BS_factory(L, K):
    d = 0
    X_list_b = []
    Y_list_b = []
    n_b = Nu//K #no. of basestations to distribute
#    while d < len(L) and n_b > 0:
    for tup in L:
        d = tup[1]
        # print("weight, cell_num ", tup)
        X_min = -5000 + (Lx/p)*(d%p)
        X_max = 5000 - (Lx/p)*((p-1) - d%p)
        #print("X_min, X_max", X_min, X_max)
        Y_min = 5000 - (Ly/p)*(p - d//p)
        Y_max = Y_min + (Ly/p)
        #print("Y_min, Y_max", Y_min, Y_max)

        
        num_bs_in_cell = tup[0]//K
        if num_bs_in_cell == 0:
            num_bs_in_cell = 1
        if num_bs_in_cell > n_b:
            num_bs_in_cell = n_b

        # print("num_bs_in_cell",num_bs_in_cell)
        for i in range(num_bs_in_cell):
            X_list_b.append(r.uniform(X_min, X_max))
            Y_list_b.append(r.uniform(Y_min, Y_max))

        n_b = n_b - num_bs_in_cell

    i = 0
    while n_b > 0:
        d = L[i][1]
        X_min = -5000 + (Lx/p)*(d%p)
        X_max = 5000 - (Lx/p)*((p-1) - d%p)
        #print("X_min, X_max", X_min, X_max)
        Y_min = 5000 - (Ly/p)*(p - d//p)
        Y_max = Y_min + (Ly/p)
        #print("Y_min, Y_max", Y_min, Y_max)

        X_list_b.append(r.uniform(X_min, X_max))
        Y_list_b.append(r.uniform(Y_min, Y_max))

        n_b -= 1
        i += 1

    # print("^^^^^^^^^^^^^^ n_b = ", n_b)
    # print("X_list_b = ", X_list_b)
    # print("---"*10)
    # print("Y_list_b = ", Y_list_b)
    return sp.array(X_list_b), sp.array(Y_list_b)

#----------------------------------------------

def PSO_2(pop, USlist, gbest, MaxIter, c1, c2, w, wdamp, Pt_max, K):
    """
    c1: personal learning coefficient
    c2: global learning coefficient
    w: weight -> speed of convergance
    wdamp: 
    """

    X_list_u = USlist[0]
    Y_list_u = USlist[1]    
    
    # flag specifies the utility function
    flag = 1
    check = 0

    sp.random.seed()
    PopSize = len(pop)    
    best_array = []
    best_cost_array = []
    # PSO_2 Loop
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
            pop[i]['cost'] = CostFunctionK(pop[i]['position'], [X_list_u, Y_list_u], flag, check, Pt_max);
            
            if pop[i]['cost'] > pop[i]['best_cost'] and (flag == 1 or CostFunctionK(pop[i]['position'], [X_list_u, Y_list_u], flag, 1, Pt_max)):
                pop[i]['best_position'] = pop[i]['position'].copy();
                pop[i]['best_cost'] = pop[i]['cost'];

                if pop[i]['best_cost'] > gbest['cost']:
                    gbest['position'] = pop[i]['best_position'].copy();
                    gbest['cost'] = pop[i]['best_cost'];

        print('Iteration {} (Utility {}) : Best Cost = {}'.format(it+1, flag, gbest['cost']));


        if flag == 1 and gbest['cost'] - xi >= 0: #when Utility 1 satisfies
            flag = 2
            gbest['cost'] = sp.NINF
            print("Iteration [",it+1,"] Best Cost Utility 1 (Avg_Cov_Prob) = ",gbest['cost'])
            for j in range(0, PopSize):   #refresh velocities
                pop[j]['velocity'] = sp.zeros((3,Nu//K))
        
                pop[j]['cost'] = CostFunctionK(pop[j]['position'], [X_list_u, Y_list_u], flag, check, Pt_max)
                pop[j]['best_position'] = pop[j]['position'].copy()
                pop[j]['best_cost'] = pop[j]['cost']

                print("Cost of",j+1,"th particle (Utility",flag,") =" ,pop[j]['cost'])

                if pop[j]['best_cost'] > gbest['cost'] and CostFunctionK(pop[j]['position'], [X_list_u, Y_list_u], flag, 1, Pt_max):
                    gbest['position'] = pop[j]['best_position'].copy()
                    gbest['cost'] = pop[j]['best_cost']


        if flag == 2 and gbest['cost'] - beta >= 0:
            print("Iteration {",it+1,"} Best Cost Utility 2 (Avg_Rate) = ",gbest['cost'])
            best_array.append(gbest['position'])
            best_cost_array.append(gbest['cost'])
            break

        w *= wdamp;
        
        best_array.append(gbest['position'])
        best_cost_array.append(gbest['cost'])

        if it == MaxIter - 1:
            print("COULDN'T SATISFY CONTRAINTS\n")

    return purge([X_list_u, Y_list_u], best_array[-1], Pt_max), best_array

#--------------------------------------------------------------

def purge(Arr_u, Arr_b, Pt_max): #global_best_BS_locations: should have BS positions

    X_list_b = Arr_b[0]
    Y_list_b = Arr_b[1]
    Z_list_b = Arr_b[2]
    
    #Nb BSs initially taken - uniformly random
    BS_list = [BaseStation(p[0], p[1], p[2]) for p in zip(X_list_b, Y_list_b, Z_list_b)]

    X_list_u = Arr_u[0]
    Y_list_u = Arr_u[1]

    #Nu users uniformly distributed
    US_list = [User(p[0], p[1]) for p in zip(X_list_u, Y_list_u) ]
    

    i = 0
    while i < len(BS_list):
        print("i = ", i)
        BS_locs = (BS_list[:i] + BS_list[i+1:]).copy()
        S = System(US_list, BS_locs)

        
        if S.avg_cov_prob(Pt_max) - xi >= 0 and S.avg_rate(Pt_max) - beta >= 0:
            BS_list = BS_locs
            #i = 0
        else:
            i += 1

        if S.num_bs == 1:
            return BS_list

    return BS_list

#--------------------------------------------------------------
                                
if __name__ == "__main__":

    beta = 1e6
    best_array, best_cost_array, pop, Arr_u = PSO_2(30, 50, 1.3962, 1.3962, 0.3298, 1.0, K)

    print("best array = ", best_array)
    Arr_b = best_array[-1]
    BS_list = purge(Arr_u, Arr_b)
    print(len(BS_list))

    Xs = sp.array([b.x for b in BS_list])
    Ys = sp.array([b.y for b in BS_list])
    Zs = sp.array([b.z for b in BS_list])

    XYZs = sp.array((Xs, Ys, Zs))
    best_array.append(XYZs)


    with open('BS_array.txt', 'w') as thefile:
        thefile.write(str(best_array) + '\n')

    with open("best_cost_array.txt", 'w') as best_cost_file:
        best_cost_file.write(str(best_cost_array))
    
