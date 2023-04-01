find_package(MPI REQUIRED)
# find_library(openmpi PATHS $ENV{LD_LIBRARY_PATH}REQUIRED)

add_library(mdpatlib
    bcastContainers.cpp
    initNode.cpp
    permuteDims.cpp
    readInput.cpp
    readTrajectories.cpp
    splitValues.cpp
    selectFromArray.cpp
)
target_link_libraries(mdpatlib PUBLIC MPI::MPI_CXX)

add_executable(mdpat main.cpp)

target_link_libraries(mdpat PUBLIC mdpatlib)
