#!/usr/bin/env python
"""
Algo
in 0
create A, B
in not 0 
Create none

broadcast A
Scatter B 

Pad B : np.pad(B, ((rang, rang - 1- N),(0,0)), 'constant') 
# padding is not working*

C = A @ B

sum (gather C)

"""

import mpi4py.MPI as MPI
import numpy as np

comm = MPI.COMM_WORLD
rang = comm.Get_rank()
nprocs = comm.Get_size()

N = 3

if rang== 0:
	A = np.random.random((N,N))
        B = np.random.random((N,N))
else:
     	A, B = None, None

A = comm.bcast(A, root=0)
B = comm.scatter(B, root=0)

C = np.empty(A.shape)
for i in range(len(A)):
        for j in range(len(B)):
                C[i,j] = np.dot(A[:,i], B)
C = sum(comm.gather(C, root=0))
print(rang, "has", C)
