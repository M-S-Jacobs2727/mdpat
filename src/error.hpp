#pragma once

#include <iostream>

#include <mpi.h>

namespace MDPAT
{
    enum class Error {NONE, IOERROR, SYNTAXERROR, ARGUMENTERROR};

    template<typename... Args>
    void errorAll(const Error error, const char message[], Args... args)
    {
        if (me == 0)
            errorOne(error, message, args...);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    template<typename... Args>
    void errorOne(const Error error, const char message[], Args... args)
    {
        char output[1024];
        snprintf(output, 1023, message, args...);
        std::cerr << output << '\n';
        MPI_Abort(MPI_COMM_WORLD, code);
    }
}
