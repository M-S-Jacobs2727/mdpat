#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace MDPAT
{
template <typename T>
int writeColumns(const std::vector<T> & values,
                 int numColumns,
                 const std::filesystem::path outfile,
                 char delimiter=' ');

template <typename T>
int writeTable(const std::vector<T> & xvalues,
               const std::vector<T> & yvalues,
               const std::vector<T> & values,
               const std::filesystem::path outfile);
}