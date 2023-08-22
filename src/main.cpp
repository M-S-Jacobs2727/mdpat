// get function to perform as an argument. leave room for multiple functions in one call?
#include "readInput.hpp"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MDPAT::InputReader input(argv[1]);
    input.runFile();
    MPI_Finalize();
    return 0;
}
