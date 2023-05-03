#include "output.hpp"

namespace fs = std::filesystem;
using std::vector;

namespace MDPAT
{
    template <typename T>
    int writeColumns(
        const std::vector<T> &values,
        const int numColumns,
        const std::filesystem::path outfile,
        const char delimiter)
    {
        std::ofstream out(outfile);
        if (out.bad())
        {
            std::cerr << "Couldn't open file " << outfile << " for writing.\n";
            return -1;
        }
        if (numColumns <= 0)
        {
            std::cerr << "Invalid number of columns: " << numColumns << '\n';
            return -2;
        }

        int numRows = values.size() / numColumns;
        if (numRows == 0 || numColumns * numRows != values.size())
        {
            std::cerr << "Number of columns (" << numColumns
                      << ") not compatible with size of data ("
                      << values.size() << ").\n";
            return -3;
        }

        for (int i = 0; i < numRows; ++i)
        {
            for (int j = 0; j < numColumns - 1; ++j)
                out << values[j * numRows + i] << delimiter;

            out << values[(numColumns - 1) * numRows + i] << '\n';
        }

        return 0;
    }

    template <typename T>
    int writeTable(
        const std::vector<T> &xvalues,
        const std::vector<T> &yvalues,
        const std::vector<T> &values,
        const std::filesystem::path outfile)
    {
        return 0;
    }

} // namespace MDPAT
