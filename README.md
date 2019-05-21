# Implementation of the ROCA attack (CVE-2017-15361)
This is the implementation of the paper [Return of the Coppersmith attack](https://roca.crocs.fi.muni.cz/).

The implementation is in python 2.7 and uses the Howgrave-Graham code from [RSA-and-LLL-attacks](https://github.com/mimoo/RSA-and-LLL-attacks).

For the detection of vulnerable keys, the code from the original authors of the paper is used (detect.py) [crocs-muni](https://github.com/crocs-muni/roca)

# Usage
```
$ python roca.py <path to key> -j <number of cores>
```
```
$ python optimization.py <path to key> -j <number of cores>
```
# Test
```
$ cd src
$ cp ../test/test_roca.py .
$ python test_roca.py
```
# Optimization
The optimization is based on an analysis of properties observed from real keys exported from affected card (Infineon JavaCard SLJ52GCA150). The parameters a' and k' turn out to have a lower entropy than stated in the original paper. This is used to change the iteration order and speedup the attack.


# HPC
```
$ sbatch slurm.sh <path to key>
```

# Organization of the Code
data/512.pem                 -> Vulnerable RSA 512-bit key

LICENSE                     -> APACHE 2.0

README.md                   -> This file

requirements.txt            -> Python 2.7 requirements

HPC/optimization_hpc.py -> modified attack implemetation for HPC

HPC/slurm.sh            -> SLURM start script for the attack

HPC/split_iteration.py  -> helper for SCRUM

src/optimization.py         -> Optimized ROCA attack

src/params.py               -> Used to calculate the parameters

src/roca.py                 -> Non-optimized attack

src/roca.sage               -> Pure sage version of the attack

src/detect.py               -> [crocs-muni](https://github.com/crocs-muni/roca)

test/test_roca.py            -> Test the attack
