# Imports
from roca import *
from params import *
from sage.all_cmdline import *
from decimal import Decimal
import time

# Define vulnerable keys
keylength = 2048
#p512 = 90100284916570139638439028489735281413386074204403670730318637933879173958391
#q512 = 89316585341239200103164798629168230185983346548558222489515734210664382833579
p = 135879036921529661794648581653002330298301044224526679653380767028908108252308273197382392628515754461497140112085352276569074111872088188367336757057332590938346879044292991775026289443754785127606230777989486075849384095736865778026395017314284500188674246388734465652666728075877428904646726042443084490733
q = 136030166899916836494910593158841550636266310029556929683174827580476574762487106877006810987126725903225945843864212303796002840361299997548544768590518964089753416844749381816714973552330950849352052797513575852750175731227705787558580111648617297716840633123746097117675990517245812564173658065172087693179
#p = int(roca_primegen(param[keylength]['M'], param[keylength]['k_max'], param[keylength]['a_max']))
#q = int(roca_primegen(param[keylength]['M'], param[keylength]['k_max'], param[keylength]['a_max']))

# Define params fo ROCA
M = param[keylength]['M']
M_prime = param[keylength]['M_prime']
a, k = get_a_k(p, M)
b, l = get_a_k(q, M)
a_prime, k_prime = get_a_k(p, M_prime)
b_prime, l_prime = get_a_k(q, M_prime)
N= p*q

print "[+] Found RSA-%i key" % keylength
print "[+] N = %i\n" % N
print "[+] M' = %i\n" % M_prime
print "[+] p = %i\n" % p
print "[+] q = %i\n" % q
print "[+] a = %i, k = %i, b = %i , l = %i\n" % (a, k, b, l)
print "[+] a' = %i, k' = %i, b' = %i , l' = %i\n" % (a_prime, k_prime, b_prime, l_prime)

# Construct polynomial
mm = param[keylength]['m'] 
tt = param[keylength]['t'] 
beta = 0.5 
XX = int((2*pow(N, Decimal(beta))) / M_prime)
m_inv = int(inverse_mod(M_prime, N))
k_tmp = int(pow(65537, a_prime, M_prime))
known_part_pol = int(k_tmp * m_inv)
F = PolynomialRing(Zmod(N), implementation='NTL', names=('x',))
(x,) = F._first_ngens(1)
pol = x + known_part_pol

# Factor
start = time.time() 
roots = coppersmith_howgrave_univariate(pol, N, beta, mm, tt, XX)
print "[+] Took: ", time.time() - start
for root in roots:
  factor1 = k_tmp + abs(root) * M_prime
  if mod(N, factor1) == 0:
    factor2 = N // factor1
    print "\n############### FACTORED ################"
    print "[+] p =" , factor1
    print "[+] q =" , factor2
    print "#########################################"
    assert (p == factor1 and q == factor2) or (q == factor1 and p == factor2)
