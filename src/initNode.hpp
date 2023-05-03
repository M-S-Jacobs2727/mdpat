#pragma once

#include <iostream>

#include <mpi.h>
#include <openacc.h>

namespace MDPAT
{
    void initNode(
        int argc,
        char **argv,
        int &me,
        int &nprocs);
}