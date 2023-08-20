#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

namespace MDPAT
{
    void bcast(
        std::string st,
        int source,
        MPI_Comm comm);

    void bcast(
        std::vector<std::string> vec,
        int source,
        MPI_Comm comm);

    template <typename T>
    void bcast(
        std::vector<T> vec,
        MPI_Datatype datatype,
        int source,
        MPI_Comm comm);

    template <typename T>
    void bcast(
        std::vector<std::vector<T>> vec,
        MPI_Datatype datatype,
        int source,
        MPI_Comm comm);
}