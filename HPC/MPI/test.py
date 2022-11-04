#!/usr/bin/env python
"""
for python3
"""
import mpi4py.MPI as MPI
import numpy as np

comm = MPI.COMM_WORLD
rang = comm.Get_rank()
nprocs = comm.Get_size()

data = np.array([1,2,3,4,6],dtype='i')
if rang == master:
	comm.Sendrecv([data,MPI.INT],dest=suiv,recvbuf=recv_data,source=prec)
	print("Process ",rang,"has sent ",data[0],"and received ",recv_data[0])
