import math as m

#-------------

# slave_limit = 35
K = 38 #Nu/Nb
#-------------
Lx = 10000 #width of the playground
Ly = 10000 #length of the playground
Lz = 600 #maximum height of a UAV
Lz_min = 100 #minimum height of a UAV

Nu = 1000 #number of users

B = 20 * 10**6 #bandwidth MHz
spec_eff = 1.7 #bps/Hz

# Nb = Nu//K #number of basestations -- to be corrected
#----------
#environment variables
C = 9.61
D = 0.16
c = 3e8 #vel of light
fc = 2e9
#----------

etaLOS = 1
etaNLOS = 20
gamma = 0.2
#---------

beta = 10**6 # target (minimum) average rate per user
xi = 0.95 #target (minimum) coverage probability
#--------
#s = 2 # pathloss exponent
N0 = 10**(-12) #noise variance
#https://ieeexplore.ieee.org/abstract/document/7511353/

# PT = 0.5
#--------------------------------------------
