# Hardcoded parameters for efficiency
# Found using params.py
param = \
{
  512: {
    "n": 39,
    "a_max": 62,
    "k_max": 37,
    "M": 0x924cba6ae99dfa084537facc54948df0c23da044d8cabe0edd75bc6,
    "M_prime": 0x1b3e6c9433a7735fa5fc479ffe4027e13bea,
    "m": 5,
    "t": 6,
    "c_a": 0x80000
  },
  1024: {
    "n": 71,
    "a_max": 134,
    "k_max": 37,
    "M": 0x7923ba25d1263232812ac930e9683ac0b02180c32bae1d77aa950c4a18a4e660db8cc90384a394940593408f192de1a05e1b61673ac499416088382,
    "M_prime": 0x24683144f41188c2b1d6a217f81f12888e4e6513c43f3f60e72af8bd9728807483425d1e,
    "m": 4,
    "t": 5,
    "c_a": 0x40000000
  },
  2048: {
    "n": 126,
    "a_max": 434,
    "k_max": 53,
    "M": 0x7cda79f57f60a9b65478052f383ad7dadb714b4f4ac069997c7ff23d34d075fca08fdf20f95fbc5f0a981d65c3a3ee7ff74d769da52e948d6b0270dd736ef61fa99a54f80fb22091b055885dc22b9f17562778dfb2aeac87f51de339f71731d207c0af3244d35129feba028a48402247f4ba1d2b6d0755baff6,
    "M_prime": 0x16928dc3e47b44daf289a60e80e1fc6bd7648d7ef60d1890f3e0a9455efe0abdb7a748131413cebd2e36a76a355c1b664be462e115ac330f9c13344f8f3d1034a02c23396e6,
    "m": 7,
    "t": 8,
    "c_a": 0x400000000
  }
}

# https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/coppersmith.sage
def coppersmith_howgrave_univariate(pol, N, beta, mm, tt, XX):
    """
    Coppersmith revisited by Howgrave-Graham
    
    finds a solution if:
    * b|N, b >= N^beta , 0 < beta <= 1
    * |x| < XX
    """
    #
    # init
    #
    dd = pol.degree()
    nn = dd * mm + tt
    
    #
    # checks
    #
    if not 0 < beta <= 1 :
        raise ValueError("beta should belongs in (0, 1]")

    if not pol.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    
    #
    # Coppersmith revisited algo for univariate
    #

    # change ring of pol and x
    polZ = pol.change_ring(ZZ)
    x = polZ.parent().gen()

    # compute polynomials
    gg = []
    for ii in range(mm):
        for jj in range(dd):
            gg.append((x * XX)**jj * N**(mm - ii) * polZ(x * XX)**ii)
    for ii in range(tt):
        gg.append((x * XX)**ii * polZ(x * XX)**mm)
    
    # construct lattice B
    BB = Matrix(ZZ, nn)
    
    for ii in range(nn):
        for jj in range(ii+1):
            BB[ii, jj] = gg[ii][jj]

    # LLL
    BB = BB.LLL(early_red=True, use_siegel=True)

    # transform shortest vector in polynomial    
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[0, ii] / XX**ii

    # factor polynomial
    potential_roots = new_pol.roots()

    return [i[0] for i in potential_roots]

# Top level of the attack, feeds the queue for the workers
def roca(N):
  
  # Key is not always of perfect size, infer from size
  keylength = int(log(N, 2))
  if keylength < 1000 :
    keylength = 512
  elif  keylength < 2000 :
    keylength = 1024 
  elif keylength < 4000 :
    keylength = 2048 
  else:
    keylength = 4096 
  
  # bruteforce
  M_prime = param[keylength]['M_prime']
  c_prime = discrete_log(N, Mod(65537, M_prime))
  ord_prime = Zmod(M_prime)(65537).multiplicative_order()
  top = (c_prime + ord_prime)/2
  beta = 0.5 
  mm = param[keylength]['m']
  tt = param[keylength]['t']

  XX = int((2*pow(N, beta)) / M_prime) 

  # Bruteforce until p, q are found
  a_prime = floor(c_prime/2)
  while a_prime < top:
      
      # Construct polynomial
      m_inv = int(inverse_mod(M_prime, N))
      k_tmp = int(pow(65537, a_prime, M_prime))
      known_part_pol = int(k_tmp * m_inv)
      F = PolynomialRing(Zmod(N), implementation='NTL', names=('x',))
      (x,) = F._first_ngens(1)
      pol = x + known_part_pol
      
      # Get roots of polynomial using coppersmith
      roots = coppersmith_howgrave_univariate(pol, N, beta, mm, tt, XX)
     
      # Check if roots are p, q
      for root in roots:
        factor1 = k_tmp + abs(root) * M_prime
        if mod(N, factor1) == 0:
          factor2 = N // factor1
          return int(factor1), int(factor2)
      a_prime += 1


# Roca
p = 135879036921529661794648581653002330298301044224526679653380767028908108252308273197382392628515754461497140112085352276569074111872088188367336757057332590938346879044292991775026289443754785127606230777989486075849384095736865778026395017314284500188674246388734465652666728075877428904646726042443084490733
q = 136030166899916836494910593158841550636266310029556929683174827580476574762487106877006810987126725903225945843864212303796002840361299997548544768590518964089753416844749381816714973552330950849352052797513575852750175731227705787558580111648617297716840633123746097117675990517245812564173658065172087693179
N = p*q
print ("[+] Factoring %i" % N)

factor1, factor2 = roca(N)

print ("[+] Found factors of N:")
print ("[+] p =" , factor1)
print ("[+] q =" , factor2)

