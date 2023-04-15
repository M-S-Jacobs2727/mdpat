#include "bcastContainers.hpp"

namespace MDPAT
{
void bcast(std::string st, int source, MPI_Comm comm) {
    int strlen = st.length();
    MPI_Bcast(&strlen, 1, MPI_INT, source, comm);
    st.resize(strlen);
    MPI_Bcast(const_cast<char*>(st.data()), strlen, MPI_CHAR, source, comm);
}

template<typename T>
void bcast(std::vector<T> vec, MPI_Datatype datatype, int source, MPI_Comm comm) {
    int size = vec.size();
    MPI_Bcast(&size, 1, MPI_INT, source, comm);
    vec.resize(size);
    MPI_Bcast(&vec[0], size, datatype, source, comm);
}

template<typename T>
void bcast(std::vector<std::vector<T>> vec, MPI_Datatype datatype, int source, MPI_Comm comm) {
    int dim1 = vec.size();
    MPI_Bcast(&dim1, 1, MPI_INT, source, comm);
    vec.resize(dim1);

    std::vector<int> sizes(dim1);
    for (int i = 0; i < dim1; ++i) sizes[i] = vec[i].size();

    MPI_Bcast(&sizes[0], dim1, MPI_INT, source, comm);

    for (int i = 0; i < dim1; ++i) {
        vec[i].resize(sizes[i]);
        MPI_Bcast(&vec[i][0], sizes[i], datatype, source, comm);
    }
}
}