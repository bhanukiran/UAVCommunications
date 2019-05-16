import scipy as sp
from Parameters import *

if __name__ == "__main__":
    X_list_b = sp.random.uniform(-Lx/2, Lx/2, Nu//K)
    Y_list_b = sp.random.uniform(-Ly/2, Ly/2, Nu//K)
    Z_list_b = sp.random.uniform(Lz_min, Lz, Nu//K)
    best_array = []
    pos = sp.array([ X_list_b, Y_list_b, Z_list_b])
    best_array.append(pos)
    with open('afterBS.txt', 'w') as thefile:
        thefile.write("x y label\n")
        for i in range(len(X_list_b)):
            thefile.write(str(X_list_b[i])+' '+str(Y_list_b[i])+' b\n')
