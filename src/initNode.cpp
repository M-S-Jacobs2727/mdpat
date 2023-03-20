#include "initNode.hpp"

void initNode(int argc, char **argv, int & me, int & nprocs) {
    // MPI variables
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int local_rank = 0, local_size = 1;
    MPI_Comm MPI_COMM_SHARE;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &MPI_COMM_SHARE);
    MPI_Comm_rank(MPI_COMM_SHARE, &local_rank);
    MPI_Comm_size(MPI_COMM_SHARE, &local_size);
    int nnodes = nprocs / local_size;
    
    // Set rank-gpu affinity
    int ngpus = acc_get_num_devices(acc_device_nvidia);
    std::cout << "# Number of GPUs per node: (" << me << ") " << ngpus << "\n";
    
    if (me==0){ 
        std::cout << "# Number of GPUs per node: " << ngpus << "\n";
        std::cout << "# Number of nodes: " << nnodes << "\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&ngpus, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int gpunum = (me / local_size) * ngpus + local_rank % ngpus;
    acc_set_device_num(gpunum, acc_device_nvidia);
    std::cout << "# me: " << me << ", gpunum: " << gpunum << "\n";
}
