from sage.all_cmdline import *
import random

# Hardcoded parameters for efficiency
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

def n_first_primes(n):
  primes = []
  p=0
  while len(primes)!=n:
    p+=1 
    primes = list(filter(is_prime, range(1, p+1)))
  return primes

def prime_power_divisor(x):
  b=[]
  for i in list(factor(x)):
    if i[1 ] != 1 :
      for j in range(1 ,i[1 ]+1 ):
        b.append(i[0]**j)
    else:
      b.append(i[0])
  b.sort()
  return b[::-1 ]

# Generates M according to key size
# Slow, use the parameters
def gen_M(n):
  return prod(n_first_primes(n))

# Generates the multiplicative order of 65537 mod M_prime
def gen_order_prime(M_prime, n):
  primes = n_first_primes(n)
  orders = [Zmod(i)(65537).multiplicative_order() for i in primes]
  return lcm(orders)

# Generates M' according to key size
# Slow, use the parameters
def gen_M_prime(M_prime_old, order):
  primes = prime_power_divisor(M_prime_old)
  M_prime = M_prime_old
  for p_i in primes:
    ord_p_i = Zmod(p_i)(65537).multiplicative_order()
    if order % ord_p_i !=0:
      M_prime = M_prime/p_i
  return M_prime

# Generate bad primes according to the ROCA paper
def roca_primegen(M, k_max, a_max):
  k = random.randint(2, 2**k_max)
  a = random.randint(2, 2**a_max)
  p = int(pow(65537,a,M))+(int(M)*int(k))
  while not is_prime(p):
    k = random.randint(2, 2**k_max)
    a = random.randint(2, 2**a_max)
    p = int(pow(65537,a,M))+(int(M)*int(k))
  return sage.rings.integer.Integer(p)

# Recover a and k from prime
def get_a_k(p, M):
  k = (p - p % M) / M
  a = discrete_log(int(p-k*M), Mod(65537, M))
  assert(int(pow(65537,a,M))+(int(M)*int(k))==p)
  return a, k

# Metric for calculating M'
def reward_at_cost_divisor(ord_M_prime_old, M_prime_old, ord_M_prime_new, M_prime_new):
  return RR((log(ord_M_prime_old, 2) - log(ord_M_prime_new, 2)) / (log(M_prime_old, 2) - log(M_prime_new, 2)))

# Get M' from M
def M_prime_keysize(keylength):
  M_prime_new = param[keylength]['M']
  ord_M_prime_new = Zmod(M_prime_new)(65537).multiplicative_order()
  divisors = []
  while int(log(M_prime_new, 2)) > keylength/4:
    M_prime_old = M_prime_new
    ord_M_prime_old = ord_M_prime_new
    possible_prime_power_divisor = prime_power_divisor(ord_M_prime_old)
    costs_per_prime = {}
    M_prime_candidates = {}
    for p_j in possible_prime_power_divisor: 
      M_prime_candidates[p_j] = gen_M_prime(M_prime_old, ord_M_prime_old / p_j)
      costs_per_prime[p_j] = reward_at_cost_divisor(ord_M_prime_old, M_prime_old, Zmod(M_prime_candidates[p_j])(65537).multiplicative_order(), M_prime_candidates[p_j])

    M_prime_new = M_prime_candidates[max(costs_per_prime, key=costs_per_prime.get)]
    ord_M_prime_new = ord_M_prime_new/max(costs_per_prime, key=costs_per_prime.get)
    divisors.append(max(costs_per_prime, key=costs_per_prime.get))
  return int(M_prime_old)

