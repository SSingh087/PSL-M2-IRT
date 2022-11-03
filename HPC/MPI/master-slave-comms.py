#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np

# master send to slave a values and slave
# operates on it and returns to master
	

rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
master = 0
tagm = 101
tagr = 201

if rank == master:
	n = int(12)
	for i in rank(1, nprocs):
		MPI.COMM_WORLD.send(n,dest=i,tag=tagm)
	for i in rank(1, nprocs):
		n = MPI.COMM_WORLD.recv(source=i,tag=tagr)
		print("Master ",rank," received from slave ",i,": n=",n)
else:
	n = MPI.COMM_WORLD.recv(source=0,tag=tagm)
	print("Slave ",rank," has received n=",n," from ",master)
	n = n*rank
	MPI.COMM_WORLD.send(n,dest=0,tag=tagr)
