#!/usr/bin/env python
"""
for python3
"""

import mpi4py.MPI as MPI
import numpy as np

comm = MPI.COMM_WORLD
rang = comm.Get_rank()

if (comm.Get_size() != 2):
	if (rang == 0):
		print(’Stop.\n Required: 2 MPI processes’)
exit()

n = 10
x = np.arange(n,dtype=np.float64)

data = (1+rang)*x
if rang == 0:
	req = comm.Isend([data,n,MPI.DOUBLE],dest=1,tag=101)
elif rang == 1:
	print(’on task’,rang,’before recv: data = ’,data)
	req = comm.Irecv(data,source=0,tag=101)
	re = False
	while re == False :
		re = MPI.Request.Test(req)
	print(’test result’,re)
	re = MPI.Request.Wait(req)
	print(’wait result’,re)
	print(’on task’,rang,’after recv: data = ’,data)
