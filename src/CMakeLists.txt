find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
# find_library(openmpi PATHS $ENV{LD_LIBRARY_PATH}REQUIRED)

find_file(MPI_CXX NAMES mpiCC mpicxx mpic++)
set(CMAKE_CXX_COMPILER ${MPI_CXX})

add_library(mdpatlib
    bcastContainers.cpp
    initNode.cpp
    output.cpp
    permuteDims.cpp
    readTrajectories.cpp
    selectFromArray.cpp
    splitValues.cpp
    msd.cpp
    readInput.cpp
)
target_compile_options(mdpatlib PRIVATE -Wall)
target_link_libraries(mdpatlib PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX)

add_executable(mdpat main.cpp)

target_include_directories(mdpat PRIVATE $ENV{HOME}/.local/include/**)
# target_include_directories(mdpat PRIVATE "/home/jacobs/software/HIP/include/**")

target_link_libraries(mdpat PUBLIC mdpatlib)
