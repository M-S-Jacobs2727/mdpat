#include "readTrajectories.hpp"

using std::vector;
using std::string;
using std::cerr;

namespace MDPAT
{
Frame<double> readDumpBinary(string filepath) {
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
    if (!dumpFile) {
        cerr << "ERROR: Could not open " << filepath << "\n";
        return dump;
    }

    // detect newer format
    char *magic_string = nullptr;
    char *columns = nullptr;
    char *unit_style = nullptr;

    // loop over snapshots in file

    while (true) {
        int endian = 0x0001;
        int revision = 0x0001;

        dumpFile.read(reinterpret_cast<char*>(&ntimestep), sizeof(long int));

        // detect end-of-file

        if (dumpFile.eof()) {
            dumpFile.close();
            break;
        }

        // detect newer format
        if (ntimestep < 0) {
            // first bigint encodes negative format name length
            long int magic_string_len = -ntimestep;

            delete[] magic_string;
            magic_string = new char[magic_string_len + 1];
            dumpFile.read(magic_string, sizeof(char)*magic_string_len);
            magic_string[magic_string_len] = '\0';

            // read endian flag
            dumpFile.read(reinterpret_cast<char*>(&endian), sizeof(int));

            // read revision number
            dumpFile.read(reinterpret_cast<char*>(&revision), sizeof(int));

            // read the real ntimestep
            dumpFile.read(reinterpret_cast<char*>(&ntimestep), sizeof(long int));
        }

        dumpFile.read(reinterpret_cast<char*>(&(dump.nAtoms)), sizeof(long int));
        dumpFile.read(reinterpret_cast<char*>(&triclinic), sizeof(int));
        dumpFile.read(reinterpret_cast<char*>(boundary), sizeof(int)*6);
        dumpFile.read(reinterpret_cast<char*>(dump.box), sizeof(double)*6);
        if (triclinic) {
            dumpFile.read(reinterpret_cast<char*>(&xy), sizeof(double));
            dumpFile.read(reinterpret_cast<char*>(&xz), sizeof(double));
            dumpFile.read(reinterpret_cast<char*>(&yz), sizeof(double));
        }
        dumpFile.read(reinterpret_cast<char*>(&(dump.nCols)), sizeof(int));

        if (magic_string && revision > 0x0001) {
            // newer format includes units string, columns string
            // and time
            int len = 0;
            dumpFile.read(reinterpret_cast<char*>(&len), sizeof(int));

            if (len > 0) {
                // has units
                delete[] unit_style;
                unit_style = new char[len + 1];
                dumpFile.read(unit_style, sizeof(char)*len);
                unit_style[len] = '\0';
            }

            char flag = 0;
            dumpFile.read(&flag, sizeof(char));

            if (flag) {
                double time;
                dumpFile.read(reinterpret_cast<char*>(&time), sizeof(double));
            }

            dumpFile.read(reinterpret_cast<char*>(&len), sizeof(int));
            delete[] columns;
            columns = new char[len + 1];
            dumpFile.read(columns, sizeof(char)*len);
            columns[len] = '\0';
        }

        dumpFile.read(reinterpret_cast<char*>(&nchunk), sizeof(int));

        // loop over processor chunks in file
        dump.atoms.reserve(dump.atoms.size() + dump.nAtoms*dump.nCols);
        for (i = 0; i < nchunk; i++) {
            dumpFile.read(reinterpret_cast<char*>(&n), sizeof(int));

            // extend buffer to fit chunk size

            if (n > maxbuf) {
                buf.resize(n);
                maxbuf = n;
            }

            // read chunk and write as nCols values per line
            dumpFile.read(reinterpret_cast<char*>(buf.data()), sizeof(double)*n);
            dump.atoms.insert(dump.atoms.begin(), buf.begin(), buf.end());
        }
    }
    delete[] columns;
    delete[] magic_string;
    delete[] unit_style;
    vector<double>().swap(buf);
    return dump;
}

string getFilename(string directory, int step) {
    std::stringstream filepath(directory);
    filepath << "dump.";
    filepath << std::setfill('0') << std::setw(9) << step;
    filepath << ".txt";
    return filepath.str();
}

template<typename T>
int readAtomsSection(std::istream & dumpFile, vector<int> columnFlag, int nCols, long nAtoms, vector<T> & atoms) {
    atoms.resize(nAtoms*nCols);

    string tmp;
    long id, offset;
    int j;

    for (long i = 0; i < nAtoms; ++i) {
        dumpFile >> id;
        offset = (id-1) * nCols;
        j = -1;
        for (auto flag : columnFlag) {
            if (flag)
                dumpFile >> atoms[offset + (++j)];
            else
                dumpFile >> tmp;
        }
    }
    return 0;
}

template<typename T>
int readDumpTextFrame(std::istream & dumpFile, vector<int>columnFlag, int nCols, Frame<T> * frame) {
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
            dumpFile.ignore(1000, '\n');  // Skip the column labels
            break;  // This will be the last section of the frame
        } else {
            while (word != "ITEM:" && dumpFile.good())
                dumpFile >> word;
        }
    }
    auto err = readAtomsSection(dumpFile, columnFlag, nCols, frame->nAtoms, frame->atoms);
    if (err) {
        cerr << "Error reading file!\n";
        return err;
    }

    return 0;
}

template<typename T>
int readDumpTextFile(std::istream & dumpFile, vector<int> columnFlag, vector<Frame<T>*> & dumpFrames) {
    int nCols = std::count(columnFlag.begin(), columnFlag.end(), 1);

    for (Frame<T> * frame : dumpFrames) {
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

template<typename T>
vector<Frame<T>*> readDumpFiles(
    int firstStep, int nSteps, int dumpStep, int nAtoms, int totalCols,
    vector<int> columns, string directory
) {
    vector<Frame<T>*> dumpFrames;
    dumpFrames.reserve(nSteps);

    vector<int> columnFlag(totalCols, 0);
    for (auto c : columns)
        columnFlag[c] = 1;

    string filename;        
    for (int step = firstStep, i = 0; i < nSteps; step += dumpStep, ++i) {
        filename = getFilename(directory, step);
        std::ifstream dumpFile(filename);
        if (dumpFile.bad()) break;

        auto err = readDumpTextFile<T>(dumpFile, columnFlag, dumpFrames);

        if (err) break;
    }
    return dumpFrames;
}

template int readAtomsSection(std::istream & dumpFile, vector<int> columnFlag, int nCols, long nAtoms, std::vector<float> & atoms);
template int readAtomsSection(std::istream & dumpFile, vector<int> columnFlag, int nCols, long nAtoms, std::vector<double> & atoms);

template int readDumpTextFrame(std::istream & dumpFile, vector<int>columnFlag, int nCols, Frame<float> * dump);
template int readDumpTextFrame(std::istream & dumpFile, vector<int>columnFlag, int nCols, Frame<double> * dump);

template int readDumpTextFile(std::istream & dumpFile, vector<int> columnFlag, vector<Frame<float>*> & dumpFrames);
template int readDumpTextFile(std::istream & dumpFile, vector<int> columnFlag, vector<Frame<double>*> & dumpFrames);

}