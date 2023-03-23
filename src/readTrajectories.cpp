#include "readTrajectories.hpp"

using std::vector;
using std::string;
using std::cerr;

Dump<double> readDumpBinary(string filepath) {
    int i, j, k, m, n;
    int nchunk, triclinic, oldsize;
    long int ntimestep;
    double xy, xz, yz;
    int boundary[3][2];
    char boundstr[9];

    int maxbuf = 0;
    vector<double> buf;

    Dump<double> dump;
    dump.nAtoms = 0L;
    dump.nCols = 0;
    for (i = 0; i < 6; ++i)
        dump.box[i] = 0.0;

    // loop over files

    fflush(stdout);
    std::fstream dumpfile(filepath, std::ios::binary);
    if (!dumpfile) {
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

        dumpfile.read(reinterpret_cast<char*>(&ntimestep), sizeof(long int));

        // detect end-of-file

        if (dumpfile.eof()) {
            dumpfile.close();
            break;
        }

        // detect newer format
        if (ntimestep < 0) {
            // first bigint encodes negative format name length
            long int magic_string_len = -ntimestep;

            delete[] magic_string;
            magic_string = new char[magic_string_len + 1];
            dumpfile.read(magic_string, sizeof(char)*magic_string_len);
            magic_string[magic_string_len] = '\0';

            // read endian flag
            dumpfile.read(reinterpret_cast<char*>(&endian), sizeof(int));

            // read revision number
            dumpfile.read(reinterpret_cast<char*>(&revision), sizeof(int));

            // read the real ntimestep
            dumpfile.read(reinterpret_cast<char*>(&ntimestep), sizeof(long int));
        }

        dumpfile.read(reinterpret_cast<char*>(&(dump.nAtoms)), sizeof(long int));
        dumpfile.read(reinterpret_cast<char*>(&triclinic), sizeof(int));
        dumpfile.read(reinterpret_cast<char*>(boundary), sizeof(int)*6);
        dumpfile.read(reinterpret_cast<char*>(dump.box), sizeof(double)*6);
        if (triclinic) {
            dumpfile.read(reinterpret_cast<char*>(&xy), sizeof(double));
            dumpfile.read(reinterpret_cast<char*>(&xz), sizeof(double));
            dumpfile.read(reinterpret_cast<char*>(&yz), sizeof(double));
        }
        dumpfile.read(reinterpret_cast<char*>(&(dump.nCols)), sizeof(int));

        if (magic_string && revision > 0x0001) {
            // newer format includes units string, columns string
            // and time
            int len = 0;
            dumpfile.read(reinterpret_cast<char*>(&len), sizeof(int));

            if (len > 0) {
                // has units
                delete[] unit_style;
                unit_style = new char[len + 1];
                dumpfile.read(unit_style, sizeof(char)*len);
                unit_style[len] = '\0';
            }

            char flag = 0;
            dumpfile.read(&flag, sizeof(char));

            if (flag) {
                double time;
                dumpfile.read(reinterpret_cast<char*>(&time), sizeof(double));
            }

            dumpfile.read(reinterpret_cast<char*>(&len), sizeof(int));
            delete[] columns;
            columns = new char[len + 1];
            dumpfile.read(columns, sizeof(char)*len);
            columns[len] = '\0';
        }

        dumpfile.read(reinterpret_cast<char*>(&nchunk), sizeof(int));

        // loop over processor chunks in file
        dump.atoms.reserve(dump.atoms.size() + dump.nAtoms*dump.nCols);
        for (i = 0; i < nchunk; i++) {
            dumpfile.read(reinterpret_cast<char*>(&n), sizeof(int));

            // extend buffer to fit chunk size

            if (n > maxbuf) {
                buf.resize(n);
                maxbuf = n;
            }

            // read chunk and write as nCols values per line
            dumpfile.read(reinterpret_cast<char*>(buf.data()), sizeof(double)*n);
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
void readAtomsSection(std::ifstream &dumpfile, int nAtoms, vector<int> columnFlag, Dump<T> & dump) {
    auto nCols = std::count(columnFlag.begin(), columnFlag.end(), 1);
    dump.atoms.resize(dump.atoms.size() + nAtoms*nCols);
    string tmp;
    T tmpdata;
    int id;
    for (int i = 0; i < nAtoms; ++i) {
        dumpfile >> id;
        
        for (auto c : columnFlag) {
            if (c) dumpfile >> dump.atoms[id];
            else dumpfile >> tmp;
        }
    }
}

template<typename T>
int readDumpText(const string & filepath, int nAtoms, vector<int> columnFlag, Dump<T> & dump) {
    std::ifstream dumpfile(filepath);
    if (!dumpfile) {
        cerr << "Couldn't open dump file '" << filepath << "'.\n";
        dumpfile.close();
        return 1;
    }
    // cerr << "Opened dump file '" << filepath << "'.\n";
    
    string word;
    int id;
    int i, j, ii, jj;

    dumpfile >> word;
    while (dumpfile.good()) {
        if (word != "ITEM:")
            continue;

        dumpfile >> word;
        if (word == "TIMESTEP") {
            dumpfile >> word;
        } else if (word == "BOX" && dump.box[0] == 0.0 && dump.box[1] == 0.0) {
            for (i = 0; i < 4; ++i) dumpfile >> word;
            for (i = 0; i < 6; ++i) dumpfile >> dump.box[i];
        } else if (word == "NUMBER" && dump.nAtoms == 0) {
            for (i = 0; i < 2; ++i) dumpfile >> word;
            dumpfile >> dump.nAtoms;
        } else if (word == "ATOMS") {
            dumpfile.ignore(1000, '\n');
            readAtomsSection(dumpfile, nAtoms, columnFlag, dump);
        }
        dumpfile >> word;
    }

    dumpfile.close();
    return 0;
}

template<typename T>
Dump<T> readDumpFiles(
    int firstStep, int nSteps, int dumpStep, int nAtoms, int totalCols,
    vector<int> columns, string directory
) {
    Dump<T> dump;

    vector<int> columnFlag(totalCols, 0);
    for (auto c : columns)
        columnFlag[c] = 1;

    string filename;        
    for (int step = firstStep, i = 0; i < nSteps; step += dumpStep, ++i) {
        filename = getFilename(directory, step);
        auto err = readDumpText<T>(filename, nAtoms, columnFlag, dump);
        if (err) break;
    }
    return dump;
}
