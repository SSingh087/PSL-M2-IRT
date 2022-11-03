#!/usr/bin/env python
"""
for python3
"""
import mpi4py.MPI as MPI
rank = MPI.COMM_WORLD.Get_rank()
numtasks = MPI.COMM_WORLD.Get_size()
hostname = MPI.Get_processor_name()
mess = "Hello! I am process %d of %d on %s."
print(mess % (rank, numtasks, hostname))
