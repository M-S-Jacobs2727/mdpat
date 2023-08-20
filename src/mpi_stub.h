#ifndef OMPI_MPI_H
#define OMPI_MPI_H

#define MPI_COMM_WORLD 0

void MPI_Comm_rank(int comm, int* me)
{
    *me = 0;
}

void MPI_Comm_size(int comm, int* nprocs)
{
    *nprocs = 1;
}

#endif