#include "readTrajectories.hpp"

using std::cerr;
using std::string;
using std::vector;

namespace MDPAT
{
    Frame<double> readDumpBinary(string filepath)
    {
        int i, n;
        int nchunk, triclinic;
        long int ntimestep;
        double xy, xz, yz;
        int boundary[3][2];

        int maxbuf = 0;
        vector<double> buf;

        Frame<double> dump;
        dump.nAtoms = 0L;
        dump.nCols = 0;
        for (i = 0; i < 6; ++i)
            dump.box[i] = 0.0;

        // loop over files

        fflush(stdout);
        std::fstream dumpFile(filepath, std::ios::binary);
        if (!dumpFile)
        {
            cerr << "ERROR: Could not open " << filepath << "\n";
            return dump;
        }

        // detect newer format
        char *magic_string = nullptr;
        char *columns = nullptr;
        char *unit_style = nullptr;

        // loop over snapshots in file

        while (true)
        {
            int endian = 0x0001;
            int revision = 0x0001;

            dumpFile.read(reinterpret_cast<char *>(&ntimestep), sizeof(long int));

            // detect end-of-file

            if (dumpFile.eof())
            {
                dumpFile.close();
                break;
            }

            // detect newer format
            if (ntimestep < 0)
            {
                // first bigint encodes negative format name length
                long int magic_string_len = -ntimestep;

                delete[] magic_string;
                magic_string = new char[magic_string_len + 1];
                dumpFile.read(magic_string, sizeof(char) * magic_string_len);
                magic_string[magic_string_len] = '\0';

                // read endian flag
                dumpFile.read(reinterpret_cast<char *>(&endian), sizeof(int));

                // read revision number
                dumpFile.read(reinterpret_cast<char *>(&revision), sizeof(int));

                // read the real ntimestep
                dumpFile.read(reinterpret_cast<char *>(&ntimestep), sizeof(long int));
            }

            dumpFile.read(reinterpret_cast<char *>(&(dump.nAtoms)), sizeof(long int));
            dumpFile.read(reinterpret_cast<char *>(&triclinic), sizeof(int));
            dumpFile.read(reinterpret_cast<char *>(boundary), sizeof(int) * 6);
            dumpFile.read(reinterpret_cast<char *>(dump.box), sizeof(double) * 6);
            if (triclinic)
            {
                dumpFile.read(reinterpret_cast<char *>(&xy), sizeof(double));
                dumpFile.read(reinterpret_cast<char *>(&xz), sizeof(double));
                dumpFile.read(reinterpret_cast<char *>(&yz), sizeof(double));
            }
            dumpFile.read(reinterpret_cast<char *>(&(dump.nCols)), sizeof(int));

            if (magic_string && revision > 0x0001)
            {
                // newer format includes units string, columns string
                // and time
                int len = 0;
                dumpFile.read(reinterpret_cast<char *>(&len), sizeof(int));

                if (len > 0)
                {
                    // has units
                    delete[] unit_style;
                    unit_style = new char[len + 1];
                    dumpFile.read(unit_style, sizeof(char) * len);
                    unit_style[len] = '\0';
                }

                char flag = 0;
                dumpFile.read(&flag, sizeof(char));

                if (flag)
                {
                    double time;
                    dumpFile.read(reinterpret_cast<char *>(&time), sizeof(double));
                }

                dumpFile.read(reinterpret_cast<char *>(&len), sizeof(int));
                delete[] columns;
                columns = new char[len + 1];
                dumpFile.read(columns, sizeof(char) * len);
                columns[len] = '\0';
            }

            dumpFile.read(reinterpret_cast<char *>(&nchunk), sizeof(int));

            // loop over processor chunks in file
            dump.atoms.reserve(dump.atoms.size() + dump.nAtoms * dump.nCols);
            for (i = 0; i < nchunk; i++)
            {
                dumpFile.read(reinterpret_cast<char *>(&n), sizeof(int));

                // extend buffer to fit chunk size

                if (n > maxbuf)
                {
                    buf.resize(n);
                    maxbuf = n;
                }

                // read chunk and write as nCols values per line
                dumpFile.read(reinterpret_cast<char *>(buf.data()), sizeof(double) * n);
                dump.atoms.insert(dump.atoms.begin(), buf.begin(), buf.end());
            }
        }
        delete[] columns;
        delete[] magic_string;
        delete[] unit_style;
        vector<double>().swap(buf);
        return dump;
    }

    string getFilepath(const string &directory, const int step)
    {
        std::stringstream filepath(directory);
        filepath << "dump.";
        filepath << std::setfill('0') << std::setw(9) << step;
        filepath << ".txt";
        return filepath.str();
    }

    string getFilename(const uint64_t step)
    {
        std::stringstream filename;
        filename << "dump.";
        filename << std::setfill('0') << std::setw(9) << step;
        filename << ".txt";
        return filename.str();
    }

    int32_t findCoordColumn(
        const Frame<float> & frame,
        const char coord,
        const bool wrapped)
    {
        auto colbegin = frame.columnLabels.begin();
        auto colend = frame.columnLabels.end();
        auto it = colbegin;
        while (true)
        {
            it = std::find_if(it+1, colend, [coord](string l){return l[0] == coord;});
            if (it == colend)
            {
                // std::cerr << "Couldn't find " << coord << " coordinate!\n";
                return 0;
            }
            std::string label = *it;
            auto col = (int32_t)(it-colbegin);
            if (wrapped && label[label.size()-1] != 'u')
                return col * (2*(label[1] != 's') - 1);
            else if (!wrapped && label[label.size()-1] == 'u')
                return col * (2*(label[1] != 's') - 1);
        }
        return 0;  // never touched
    }

    int readAtomsSection(
        std::istream &dumpFile,
        vector<int> columnFlag,
        uint32_t nValidCols,
        uint64_t nAtoms,
        vector<float> &atoms)
    {
        auto oldSize = atoms.size();
        atoms.resize(oldSize + nAtoms * nValidCols);

        string tmp;
        float tmpval;

        for (uint64_t i = 0; i < nAtoms; ++i)
        {
            uint64_t id, offset;
            dumpFile >> id;
            offset = (id - 1) * nValidCols;
            uint64_t j = 0;
            for (auto flag : columnFlag)
            {
                if (flag)
                {
                    dumpFile >> tmpval;
                    atoms[oldSize + offset + j] = tmpval;
                    ++j;
                }
                else
                    dumpFile >> tmp;
            }
        }
        return 0;
    }

    int readAtomsSectionByType(
        std::istream &dumpFile,
        vector<int> columnFlag,
        uint32_t nValidCols,
        vector<float> &atoms,
        vector<uint32_t> typeflag,
        uint64_t nValidAtoms)
    {
        auto oldSize = atoms.size();
        atoms.resize(oldSize + nValidAtoms * nValidCols);

        string tmp;
        uint64_t id, offset;
        uint32_t j;

        for (uint64_t i = 0; i < typeflag.size(); ++i)
        {
            dumpFile >> id;

            if (!typeflag[id - 1])
            {
                // for (j = 0; j < columnFlag.size(); ++j)
                //     dumpFile >> tmp;
                dumpFile.ignore(10000, '\n');
                continue;
            }

            offset = (id - 1) * nValidCols;
            j = 0;
            for (auto flag : columnFlag)
            {
                if (flag)
                    dumpFile >> atoms[oldSize + offset + j];
                else
                    dumpFile >> tmp;
                ++j;
            }
        }
        return 0;
    }

    Frame<float> readDumpTextHeader(std::istream &dumpFile)
    {
        string word;
        int i;

        Frame<float> frame;

        dumpFile >> word;
        if (word != "ITEM:")
        {
            cerr << "Syntax error while reading dump file!\n";
            throw -1;
        }

        while (true)
        {
            dumpFile >> word;
            if (word == "TIMESTEP")
            {
                dumpFile >> frame.timestep;
            }
            else if (word == "BOX")
            {
                dumpFile >> word >> word >> word >> word; // ITEM: BOX BOUNDS pp pp pp
                for (i = 0; i < 6; ++i)
                    dumpFile >> frame.box[i];
            }
            else if (word == "NUMBER")
            {
                dumpFile >> word >> word; // ITEM: NUMBER OF ATOMS
                dumpFile >> frame.nAtoms;
            }
            else if (word == "ATOMS")
            {
                std::getline(dumpFile, word); // ITEM: ATOMS id type x y z ...
                std::istringstream labelstream(word);
                while (labelstream)
                {
                    labelstream >> word;
                    frame.columnLabels.push_back(word);
                }
                frame.columnLabels.pop_back();
                frame.nCols = frame.columnLabels.size();
                return frame;
            }
            else
            {
                while (word != "ITEM:" && dumpFile.good())
                    dumpFile >> word;
                if (dumpFile.bad())
                    throw -3;
                continue;
            }
            dumpFile >> word;
        }

        throw -2;
    }

    void skipDumpTextHeader(std::istream &dumpFile)
    {
        string word;

        dumpFile >> word;
        if (word != "ITEM:")
        {
            cerr << "Syntax error while reading dump file!\n";
            throw -1;
        }

        while (true)
        {
            dumpFile >> word;
            if (word == "ATOMS")
            {
                dumpFile.ignore(10000, '\n');
                return;
            }
            while (word != "ITEM:" && dumpFile.good())
                dumpFile >> word;
            if (dumpFile.bad())
                throw -2;
        }
    }

    /*
     * This function assumes that
     * 1. Each frame is a separate file named "dump.{timestep}.txt", where
     *    timestep is a 9-digit integer with leading 0's
     * 2. The frames are located in a single directory
     * 3. The timesteps are evenly spaced: initStep, initStep+dumpStep, initStep+2*dumpStep, ..., endStep
     * 4. The atoms are contiguously numbered 1-nAtoms
     * 5. The number of atoms does not change
     * 6. The atom types do not change
     * It ignores any changes in box size or atom type. The types are defined by the first frame.
     */
    vector<float> getTrajectories(
        std::filesystem::path directory,
        const StepRange stepRange,
        const bool wrapped,
        const uint32_t atomType,
        const uint32_t dim)
    {
        if (dim == 0)
            ;
        if (dim > 3)
            ;

        auto initStep = stepRange.initStep;
        auto endStep = stepRange.endStep;
        auto dumpStep = stepRange.dumpStep;
        auto nSteps = stepRange.nSteps;

        // Read header of first frame
        string filename = getFilename(initStep);
        std::ifstream dumpstream((directory / filename).string());
        auto frame = readDumpTextHeader(dumpstream);
        auto colbegin = frame.columnLabels.begin();
        auto colend = frame.columnLabels.end();

        // Get columns of coordinates
        vector<int> columnFlags(frame.nCols-1, 0);

        auto ind = findCoordColumn(frame, 'x', wrapped);
        columnFlags[abs(ind)-1] = 1;

        if (dim > 1)
        {
            ind = findCoordColumn(frame, 'y', wrapped);
            columnFlags[abs(ind)-1] = 1;
        }

        if (dim > 2)
        {
            ind = findCoordColumn(frame, 'z', wrapped);
            columnFlags[abs(ind)-1] = 1;
        }

        bool shrink = ind < 0;

        vector<float> atoms;

        // Read the rest
        if (atomType == 0)
        {
            atoms.reserve(frame.nAtoms * dim * nSteps);
            dumpstream.close();
            for (uint64_t step = initStep; step <= endStep; step += dumpStep)
            {
                filename = getFilename(step);
                dumpstream.open((directory / filename).string());

                skipDumpTextHeader(dumpstream);
                readAtomsSection(dumpstream, columnFlags, dim, frame.nAtoms, atoms);

                dumpstream.close();
            }
        }
        else
        {
            // Get column of atom types
            auto it = std::find(colbegin, colend, "type");
            if (it == colend)
            {
                std::cerr << "Couldn't find atom types!\n";
                throw -2;
            }
            uint32_t typeCol = it - colbegin;

            uint32_t typeColRelative = 0;
            for (uint32_t i = 0; i < frame.nCols; ++i)
            {
                if (i == typeCol)
                    break;
                if (columnFlags[i])
                    ++typeColRelative;
            }

            // Finish reading first frame, including atom types
            atoms.reserve(frame.nAtoms * (dim + 1));
            columnFlags[typeCol] = 1;
            readAtomsSection(dumpstream, columnFlags, dim + 1, frame.nAtoms, atoms);
            columnFlags[typeCol] = 0;

            dumpstream.close();

            // Flag the atoms by ID-1 if the type is atomType
            vector<uint32_t> typeflag(frame.nAtoms, 0);
            uint64_t nValidAtoms = 0UL;
            for (uint64_t i = 0; i < frame.nAtoms; ++i)
            {
                if ((uint32_t)(atoms[i * frame.nCols + typeColRelative] + 0.5) == atomType)
                {
                    typeflag[i] = 1;
                    ++nValidAtoms;
                }
            }
            atoms.clear();
            atoms.reserve(nValidAtoms * dim * nSteps);

            for (uint64_t step = initStep; step <= endStep; step += dumpStep)
            {
                filename = getFilename(step);
                dumpstream.open((directory / filename).string());

                skipDumpTextHeader(dumpstream);
                readAtomsSectionByType(dumpstream, columnFlags, dim, atoms, typeflag, nValidAtoms);

                dumpstream.close();
            }
        }

        if (shrink)
        {
            for (uint64_t i = 0; i < atoms.size() / dim; ++i)
                for (uint32_t j = 0; j < dim; ++j)
                    atoms[i*dim + j] = atoms[i*dim + j] * (frame.box[j*2+1] - frame.box[j*2]) + frame.box[j*2];
        }

        return atoms;
    }
}
