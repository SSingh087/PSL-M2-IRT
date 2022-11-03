#!/usr/bin/env python
"""
for python3
"""
import mpi4py.MPI as MPI
import numpy as np

comm = MPI.COMM_WORLD
rang = comm.Get_rank()
nprocs = comm.Get_size()

suiv = (rang+1)%nprocs
prec = (nprocs+rang-1)%nprocs
y = 0
x = 100 + rang

recv_data = np.empty(1,dtype='i')
data = np.array([x],dtype='i')
comm.Sendrecv([data,MPI.INT],dest=suiv,recvbuf=recv_data,source=prec)

print("Process ",rang,"has sent ",data[0],"and received ",recv_data[0])
