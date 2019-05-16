from Parameters import *
import scipy as sp
import multiprocessing as mp
from math import ceil
import random as r

class BaseStation(object):
    def __init__(self, x, y, z):
        """
        Pt : transmit power
        x,y,w: position cordinates of the UAV
        """
        self.x = x
        self.y = y
        self.z = z

        self.users_serve = 0
        #self.list_of_slaves = []

#----------------------------------
class User(object):
        def __init__(self, x, y):
            """
            Pr: power received by the user
            """
            self.x = x
            self.y = y

            #crap follows
            self.rec_max_sinr =  0
            self.master = None

        @property
        def got_good_sinr(self):
            return (self.rec_max_sinr >= gamma)
#---------------------------------

class System(object):
        def __init__(self, US_list, BS_list):
                self.BSs = BS_list
                self.USs = US_list

                self.num_bs = len(BS_list)
                self.num_u = len(US_list)

#-------------------------------------------

        def rec_power_users_matrix(self, Pt):

            dist = lambda u, b: sp.sqrt((u.x - b.x)**2 + (u.y - b.y)**2 + (0 - b.z)**2)

            theta = lambda u,b: (180/sp.pi) * sp.arcsin(b.z/dist(u,b))

            Plos = lambda u, b: 1 / (1 + C * sp.exp(- D*(theta(u,b) - C)))

            # L = lambda u, b: Plos(u,b)*(dist(u,b)**s) + eta*(1 - Plos(u,b))*(dist(u,b)**s)

            L = lambda u, b: 20*sp.log10(4*fc*sp.pi*dist(u, b)/c) + Plos(u, b) * etaLOS + (1 - Plos(u, b))*etaNLOS

            rec_power = lambda u,b: 10**((10*sp.log10(Pt) - L(u, b))/10)

            res = sp.zeros((self.num_u, self.num_bs))

            i = 0
            while i < self.num_u:
                sinr_i = 0
                j = 0
                while j < self.num_bs:
                    junk = rec_power(self.USs[i], self.BSs[j])
                    res[i][j] = junk
                    j+= 1

                i += 1

            return res

#-------------------------------------------

        def associate(self, Pt):
            res = self.rec_power_users_matrix(Pt)

            for BS in self.BSs:
                BS.users_serve = 0

            for user in self.USs:
                user.master = None
                user.rec_max_sinr =  0

            i = 0
            while i < self.num_u:
                u_i = self.USs[i]
                interf = sum(res[i]) - res[i]
                sinr_i = res[i]/(N0 + interf)

                max_sinr_i = max(sinr_i)

                def find_master():
                    pocket = []
                    ix = 0
                    while ix < self.num_bs:
                        if sinr_i[ix] >= gamma:
                            pocket.append((self.BSs[ix], sinr_i[ix]))
                        ix += 1
                    res = None
                    try:
                        res = max(pocket, key=lambda x:x[1])
                    except ValueError:
                        pass
                    return res

                pot_master = find_master()
                if pot_master != None:
                    u_i.master = pot_master[0]
                    u_i.rec_max_sinr = pot_master[1]
                    pot_master[0].users_serve += 1
                

                # for pot_master in pot_masters:
                #     if pot_master[0].users_serve < slave_limit:
                #         u_i.master = pot_master[0]
                #         u_i.rec_max_sinr = pot_master[1]
                #         pot_master[0].users_serve += 1
                #         break
                #     else:
                #         pass

                i += 1

#-------------------------------------------

        # def slave_distribution(self):
        #     for bs in self.BSs:
        #         print(bs.users_serve)

        def avg_cov_prob(self, Pt):
            self.associate(Pt)
            ans = 0
            for u_i in self.USs:
                if u_i.got_good_sinr:
                    ans += 1

            return ans/self.num_u

#-------------------------------------------
        def avg_rate(self, Pt):
            self.associate(Pt)
            ans = 0
            avg = []
            for u_i in self.USs:
                if u_i.master:
                    bw = B/u_i.master.users_serve
                    #print("slaves: ", u_i.master.users_serve)
                    avg.append(bw*sp.log2(1 + u_i.rec_max_sinr))
                else:
                    avg.append(0)

            return sum(avg)/self.num_u

#--------------------------------------

        def optimal_power(self, PMAX):
            eps = 0.01

            # metric1 = lambda Pt: self.avg_cov_prob(Pt) - xi
            # PMAX = 20
            Pt_min = 0
            Pt_max = PMAX

            metric = self.avg_cov_prob(Pt_max) - xi
            if metric < -eps  :
                self.op_flag = metric
                print("crappy covprob = ", metric)
                print("---"*10)
                return sp.inf

            metric2 = self.avg_rate(Pt_max) - beta
            if metric2 < -eps:
                self.op_flag = metric2
                print("crappy avg_rate =  ", metric2)
                print("---"*10)
                return sp.inf

            N = 0; Nmax = 25
            while N < Nmax:
                mid = (Pt_min + Pt_max)/2
                metric = self.avg_cov_prob(mid) - xi
                # print("metric = ", metric)
                # print("mid = ", mid)
                if metric > eps:
                    Pt_max = mid
                elif metric < -eps:
                    Pt_min = mid
                else:
                    print("metric = 0")
                    break
                N += 1

            Pt_cp_min = mid
            print("mid1 = ", Pt_cp_min)
            Pt_min = 0
            Pt_max = PMAX

            #print("Pt_min = %f \npt_max = %f"%(Pt_min, Pt_max))
            N = 0;
            while N < Nmax:
                mid2 = (Pt_min + Pt_max)/2
                if mid2 < Pt_cp_min:
                    break
                metric = self.avg_rate(mid2) - beta
                if metric > eps:
                    Pt_max = mid2
                elif metric < -eps:
                    Pt_min = mid2
                else:
                    break
                N += 1
            # print("2cov_prob = %s; rate = %s"%(self.avg_cov_prob(Pt_cp_min), self.avg_rate(mid2)))
            Pt_r_min = mid2

            print("mid2 = ", mid2)
            print("---"*10)
            return max(Pt_cp_min, Pt_r_min)

#----------------------------------------------

def users_per_cell(X_list_u, Y_list_u):
 #   L = [[0,0,0,0] for i in range(4)]
    L = [0 for i in range(16)]

    c = 0
    while c < Nu:
        
        x = X_list_u[c]
        y = Y_list_u[c]

        i = ceil(x/2500) + 1
        j = ceil(y/2500) + 1
        
        try:
            L[4*j + i] += 1
        except IndexError:
            pass
        
        c += 1

    # print("L = ", L)
    return L

def BS_factory(L, K):
    d = 0
    X_list_b = []
    Y_list_b = []
    n_b = Nu//K #no. of basestations to distribute
    while d < len(L) and n_b > 0:
        # print("d = ", d)
        X_min = -5000 + 2500*(d%4)
        X_max = 5000 - 2500*(3 - d%4)
        # print("X_min, X_max", X_min, X_max)
        Y_min = 5000 - 2500*(4 - d//4)
        Y_max = Y_min + 2500
        # print("Y_min, Y_max", Y_min, Y_max)

        
        num_bs_in_cell = L[d]//K
        if num_bs_in_cell == 0:
            num_bs_in_cell = 1
        if num_bs_in_cell > n_b:
            num_bs_in_cell = n_b

        for i in range(num_bs_in_cell):
            X_list_b.append(r.uniform(X_min, X_max))
            Y_list_b.append(r.uniform(Y_min, Y_max))

        n_b = n_b - num_bs_in_cell
        d += 1

    # print("X_list_b = ", X_list_b)
    # print("---"*10)
    # print("Y_list_b = ", Y_list_b)
    return X_list_b, Y_list_b
#----------------------------------------------

if __name__ == "__main__":
    import pylab as p
    sp.random.seed(17081986)
    X_list_u = sp.random.uniform(-Lx/2, Lx/2, Nu)
    Y_list_u = sp.random.uniform(-Ly/2, Ly/2, Nu)
    US_list = [User(p[0], p[1]) for p in zip(X_list_u, Y_list_u) ]
    K = 30
    # X_list_b = sp.random.uniform(-Lx/2, Lx/2, Nu//K)
    # Y_list_b = sp.random.uniform(-Ly/2, Ly/2, Nu//K)
    X_list_b, Y_list_b = BS_factory(users_per_cell(X_list_u, Y_list_u), K)
    
    Z_list_b = sp.random.uniform(Lz_min, Lz, Nu//K)
    BS_list = [BaseStation(p[0], p[1], p[2]) for p in zip(X_list_b, Y_list_b, Z_list_b)]

    S = System(US_list, BS_list)
    print("coverage probability  = ",  S.avg_cov_prob(9e-2))
    print("avg rate  = ",  S.avg_rate(1))

    # S.optimal_power_bisec()

    #print("avg rate = ", S.avg_rate(1))
    #print("optimal power = ", S.optimal_power())
    #pool = mp.Pool(processes = 14)
    #results = [pool.apply(S.avg_cov_prob, args=(0.1*power,)) for power in range(1,2)]
    # results = [pool.apply(S.avg_rate, args=(power,)) for power in range(100,200, 10)]
    #
    # print(results)
    # #print(results2)
    # p.plot(range(100,200, 10), results)
    # p.show()
    # p.clear('all')
