#!/usr/bin/env python

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
        # dummy	values for transpose to	work if	rank != 0
     	A, B = np.zeros((N,N)), np.zeros((N,N))

A = comm.bcast(A, root=0)
B = comm.scatter(B.transpose(), root=0)
C = A[rang]* B

C = np.array(comm.gather(C, root=0))

print(rang, "has", C)
