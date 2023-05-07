#include "readTrajectories.hpp"

using std::cerr;
using std::string;
using std::vector;

namespace MDPAT
{
    Frame<double> readDumpBinary(string filepath)
    {
        int i, j, k, m, n;
        int nchunk, triclinic, oldsize;
        long int ntimestep;
        double xy, xz, yz;
        int boundary[3][2];
        char boundstr[9];

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

        for (uint64_t i = 0; i < nAtoms; ++i)
        {
            uint64_t id, offset;
            dumpFile >> id;
            offset = (id - 1) * nValidCols;
            uint64_t j = 0;
            for (auto flag : columnFlag)
            {
                if (flag)
                    dumpFile >> atoms[oldSize + offset + (j++)];
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
        long id, offset;
        int j;

        for (long i = 0; i < typeflag.size(); ++i)
        {
            dumpFile >> id;

            if (!typeflag[id - 1])
            {
                for (j = 0; j < columnFlag.size(); ++j)
                    dumpFile >> tmp;
                continue;
            }

            offset = (id - 1) * nValidCols;
            j = -1;
            for (auto flag : columnFlag)
            {
                if (flag)
                    dumpFile >> atoms[oldSize + offset + (++j)];
                else
                    dumpFile >> tmp;
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

    int readDumpTextFrame(std::istream & dumpFile, vector<int>columnFlag, int nCols, Frame<float> * frame) {
        string word;
        int i;

        while (true) {
            dumpFile >> word;
            if (word != "ITEM:") {
                cerr << "Syntax error while reading dump file!\n";
                return 1;
            }

            dumpFile >> word;
            if (word == "TIMESTEP") {
                dumpFile >> frame->timestep;
            } else if (word == "BOX") {
                dumpFile >> word >> word >> word >> word;  // Boundary conditions
                for (i = 0; i < 6; ++i)
                    dumpFile >> frame->box[i];
            } else if (word == "NUMBER") {
                dumpFile >> word >> word;  // ITEM: NUMBER OF ATOMS
                dumpFile >> frame->nAtoms;
            } else if (word == "ATOMS") {
                std::getline(dumpFile, word);
                std::istringstream labelstream(word);
                while (labelstream)
                {
                    labelstream >> word;
                    frame->columnLabels.push_back(word);
                }
                frame->nCols = frame->columnLabels.size();
                // dumpFile.ignore(1000, '\n');
                break;  // This will be the last section of the frame
            } else {
                while (word != "ITEM:" && dumpFile.good())
                    dumpFile >> word;
            }
        }
        auto err = readAtomsSection(dumpFile, columnFlag, nCols, frame->nAtoms, frame->atoms);
        if (err) {
            cerr << "Error reading atoms section of frame!\n";
            return err;
        }

        return 0;
    }

    int readDumpTextFile(std::istream & dumpFile, vector<int> columnFlag, vector<Frame<float>*> & dumpFrames) {
        int nCols = std::count(columnFlag.begin(), columnFlag.end(), 1);

        for (Frame<float> * frame : dumpFrames) {
            cerr << frame << '\n';
            auto err = readDumpTextFrame(dumpFile, columnFlag, nCols, frame);
            if (err) {
                cerr << "Error reading frame!\n";
                return err;
            }
        }
        cerr << "Successfully read file!\n";

        return 0;
    }

    vector<Frame<float>*> readDumpFiles(
        int firstStep, int nSteps, int dumpStep, int nAtoms, int totalCols,
        vector<int> columns, const string & directory
    ) {
        vector<Frame<float>*> dumpFrames;
        dumpFrames.resize(nSteps);

        vector<int> columnFlag(totalCols, 0);
        for (auto c : columns)
            columnFlag[c] = 1;

        string filename;
        for (int step = firstStep, i = 0; i < nSteps; step += dumpStep, ++i) {
            filename = getFilepath(directory, step);
            std::ifstream dumpFile(filename);
            if (dumpFile.bad()) break;

            auto err = readDumpTextFile(dumpFile, columnFlag, dumpFrames);

            if (err) break;
        }
        return dumpFrames;
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
    vector<float> getTrajectories(bool wrapped,
                                  StepRange stepRange,
                                  std::filesystem::path directory,
                                  uint32_t atomType,
                                  uint32_t dim)
    {
        if (dim == 0)
            ;
        if (dim > 3)
            ;

        auto initStep = stepRange.initStep;
        auto endStep = stepRange.endStep;
        auto dumpStep = stepRange.dumpStep;

        uint64_t nSteps = (stepRange.endStep - stepRange.initStep) / stepRange.dumpStep + 1;

        // Read first frame
        string filename;
        filename = getFilename(stepRange.initStep);
        std::ifstream dumpstream((directory / filename).string());
        auto frame = readDumpTextHeader(dumpstream);

        // Get column of atom types
        auto it = std::find(frame.columnLabels.begin(), frame.columnLabels.end(), "type");
        if (it == frame.columnLabels.end())
        {
            std::cerr << "Couldn't find atom types!\n";
            throw -2;
        }
        uint32_t typeCol = it - frame.columnLabels.begin();

        // Get columns of coordinates
        vector<int> columnFlags(frame.nCols, 0);

        std::string label = wrapped ? "x" : "xu";
        it = std::find(frame.columnLabels.begin(), frame.columnLabels.end(), label);
        if (it == frame.columnLabels.end())
        {
            std::cerr << "Couldn't find unwrapped coordinates!\n";
            throw -2;
        }
        columnFlags[(int)(it - frame.columnLabels.begin())] = 1;

        if (dim > 1)
        {
            label[0] = 'y';
            it = std::find(frame.columnLabels.begin(), frame.columnLabels.end(), label);
            if (it == frame.columnLabels.end())
            {
                std::cerr << "Couldn't find unwrapped coordinates!\n";
                throw -2;
            }
            columnFlags[(int)(it - frame.columnLabels.begin())] = 1;
        }

        if (dim > 2)
        {
            label[0] = 'z';
            it = std::find(frame.columnLabels.begin(), frame.columnLabels.end(), label);
            if (it == frame.columnLabels.end())
            {
                std::cerr << "Couldn't find unwrapped coordinates!\n";
                throw -2;
            }
            columnFlags[(int)(it - frame.columnLabels.begin())] = 1;
        }

        uint32_t typeColRelative = 0;
        for (uint32_t i = 0; i < frame.nCols; ++i)
        {
            if (i == typeCol)
                break;
            if (columnFlags[i])
                ++typeColRelative;
        }

        // Finish reading first frame, including atom types
        vector<float> atoms;
        atoms.reserve(frame.nAtoms * (dim + 1) * nSteps);

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
        dumpstream.seekg(0);

        // Read the rest
        for (uint64_t step = initStep; step <= endStep; step += dumpStep)
        {
            filename = getFilename(step);
            dumpstream.open((directory / filename).string());

            skipDumpTextHeader(dumpstream);
            readAtomsSectionByType(dumpstream, columnFlags, dim, atoms, typeflag, nValidAtoms);

            dumpstream.close();
        }

        return atoms;
    }
}
