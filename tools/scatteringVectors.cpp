#include "scatteringVectors.hpp"

using std::vector;

int main(int nargs, char *args[]) {
    int maxGrid = 0;
    std::string candidateFilename = "";
    int dim = 0;
    int targetLength = 0;
    std::string QFilename = "";

    int i = 0, j;

    while (i < nargs) {
        if (args[i] == "-g" || args[i] == "--max-grid") {
            maxGrid = std::stoi(args[i+1]);
            i += 2;
        } else if (args[i] == "-d" || args[i] == "--dim") {
            dim = std::stoi(args[i+1]);
            i += 2;
        } else if (args[i] == "-t" || args[i] == "--target") {
            targetLength = std::stoi(args[i+1]);
            i += 2;
        } else if (args[i] == "-c" || args[i] == "--candidate-file") {
            candidateFilename = args[i+1];
            i += 2;
        } else if (args[i] == "-o" || args[i] == "--output-file") {
            QFilename = args[i+1];
            i += 2;
        } else if (args[i] == "-h" || args[i] == "--help" || args[i] == "-?") {
            showhelp();
            return 0;
        } else {
            std::cerr << "Unrecognized option: " << args[i] << std::endl;
            showhelp();
            return 1;
        }
    }

    if (dim <= 0 || targetLength <= 0 || QFilename.size() == 0 || candidateFilename.size() == 0) {
        std::cerr << "Unspecified values! Must specify dimension," << 
            " target number of vectors, and both filenames (for candidate and final vectors)." << 
            std::endl;
        return 2;
    }

    vector<vector<int>> candidateQs;
    if (maxGrid == 0) {
        std::ifstream fpr(candidateFilename);

        if (!fpr) {
            std::cerr << "Couldn't open file: " << candidateFilename << std::endl;
            return 3;
        }

        std::string line;
        vector<int> tmpq(dim, 0);
        int ii = -1;
        while (true) {
            candidateQs.push_back(tmpq);
            fpr >> candidateQs[++ii][0];
            if (!fpr) {
                candidateQs.pop_back();
                break;
            }
            for (int i = 1; i < dim; ++i) fpr >> candidateQs[ii][i];
        }
        fpr.close();
        vector<int>().swap(tmpq);
    } else {
        candidateQs = generateCandidates(maxGrid, dim);
        writeAllVecs(candidateFilename, candidateQs);
    }

    auto selectedQs = selectVectors(candidateQs, targetLength);
    vector<vector<int>>().swap(candidateQs);

    auto finalQs = includePermutations(selectedQs);
    writeAllVecs(QFilename, finalQs);

    vector<vector<int>>().swap(selectedQs);
    vector<vector<int>>().swap(finalQs);

    return 0;
}

double norm(vector<int> vec) {
    int sum = 0;
    for (auto v : vec)
        sum += v*v;
    return sqrt((double)sum);
}

vector<double> computeNorms(vector<vector<int>> vecs) {
    vector<double> norms;
    norms.resize(vecs.size());
    std::transform(
        vecs.begin(),
        vecs.end(),
        norms.begin(),
        norm
    );
    return norms;
}

vector<vector<int>> generateCandidates(int gridMax, int dim) {
    int i, j, n, q;
    int nCandidates;

    for (i = 0; i < dim; ++i) {
        nCandidates *= gridMax + 1;
        nCandidates /= i + 1;
    }

    vector<vector<int>> candidateQs;
    vector<int> tempQ;
    tempQ.resize(dim);
    candidateQs.resize(nCandidates);

    for (i = 0; i < nCandidates; ++i) {
        candidateQs[i] = tempQ;
        
        ++tempQ[dim-1];
        for (j = dim-1; j > 0; --j) {
            if (tempQ[j] <= tempQ[j-1]) break;
            ++tempQ[j-1];
            tempQ[j] = 0;
        }
    }

    vector<std::pair<vector<int>, double>> pairs;
    pairs.resize(nCandidates);
    for (i = 0; i < nCandidates; ++i) {
        pairs[i].first = candidateQs[i];
        pairs[i].second = norm(candidateQs[i]);
    }
    std::sort(
        pairs.begin(),
        pairs.end(),
        [](auto & a, auto & b) {
            return a.second < b.second
        }    
    );

    for (i = 0; i < nCandidates; ++i)
        candidateQs[i] = pairs[i].first;

    return candidateQs;
}

vector<vector<int>> selectVectors(vector<vector<int>> candidateQs, int targetLength) {
    vector<double> norms = computeNorms(candidateQs);

    vector<vector<int>> selectedQs;
    selectedQs.reserve(targetLength);
    selectedQs.push_back(candidateQs[0]);

    double r = pow(10.0, log10(norms.back()-1.0) / (targetLength - 1));
    double binlo, binhi;

    int i, j;
    i = j = 1;
    while (j < targetLength) {
        while (i < candidateQs.size()) {
            binlo = pow(r, (double)(j-1));
            binhi = binlo * r;
            if (binlo <= norms[i] && norms[i] < binhi) break;
            if (norms[i] < binlo) ++i;
            else ++j;
        }
        selectedQs.push_back(candidateQs[i]);
        ++i;
        ++j;
    }
    return selectedQs;
}

int writeAllVecs(std::string filename, const vector<vector<int>> & vectors) {
    std::ofstream qfile(filename);
    if (!qfile) {
        std::cerr << "Couldn't open file for writing: " << filename << std::endl;
        return 1;
    }

    for (auto vec : vectors) {
        for (auto v: vec) qfile << v;
        qfile << "\n";
    }
    return 0;
}

vector<vector<int>> generatePermMatrices(int dim) {
    int i, j, k;
    int twotodim = 2, dimfactorial = 1;
    for (i = 1; i < dim; ++i) {
        twotodim *= 2;
        dimfactorial *= (i+1);
    }
    int nmatrices = twotodim * dimfactorial;

    vector<vector<int>> matrices;
    matrices.resize(nmatrices);

    vector<int> order_indices(dim, 0);
    for (i = 0; i < dim; ++i) order_indices[i] = i;
    vector<int> index_values(dim, 1);

    int ii = -1;
    for (i = 0; i < twotodim; ++i) {
        for (j = 0; j < dimfactorial; ++j) {
            matrices[++ii].resize(dim*dim);
            for (k = 0; k < dim; ++k) {
                matrices[ii][k*dim + order_indices[k]] = index_values[k];
            }
            std::next_permutation(order_indices.begin(), order_indices.end());
        }
        index_values[dim-1] *= -1;
        for (j = dim-1; j > 0; ++j) {
            if (index_values[j] == -1) break;
            index_values[j-1] *= -1;
        }
    }

    return matrices;
}

vector<int> matMul(vector<int> mat, vector<int> vec) {
    int dim = vec.size();
    vector<int> out(dim, 0);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            out[i] += mat[i*dim+j] * vec[j];
    return out;
}

vector<vector<int>> includePermutations(vector<vector<int>> selectedQs) {
    int dim = selectedQs[0].size();
    auto perm_matrices = generatePermMatrices(dim);
    vector<vector<int>> final_vectors;
    vector<int> vec(dim, 0);

    for (auto q : selectedQs) {
        auto old_end = final_vectors.end();
        final_vectors.reserve(final_vectors.size() + perm_matrices.size());
        for (auto mat : perm_matrices) {
            vec = matMul(mat, q);
            final_vectors.push_back(vec);
        }
        std::sort(
            old_end,
            final_vectors.end(),
            [](auto & a, auto & b){
                for (int i = 0; i < *a.size(); ++i) {
                    if (a[i] < b[i]) return true;
                    else if (a[i] > b[i]) return false;
                }
                return true;
            }
        );
        std::unique(old_end, final_vectors.end());
    }
    return final_vectors;
}

void showhelp() {
    std::cout << "A tool for generating scattering vectors at a logarithmic scale" <<
        " and saving the results for later use.\n\n";
    std::cout << "-d <D>\n";
    std::cout << "--dim <D>                     " << 
        "Number of columns in the file (number of dimensions in the vectors)\n";
    std::cout << "-g <N>\n";
    std::cout << "--max-grid <N>                " << 
        "Maximum value in each dimension of the vectors\n";
    std::cout << "-t <N>\n";
    std::cout << "--target <N>                  " << 
        "Target number of vectors in final output file\n";
    std::cout << "-c <filename>\n";
    std::cout << "--candidate-file <filename>   " << 
        "Filename for the intermediate/candidate vectors (all vectors" <<
        " with components less than max value).\n";
    std::cout << "-o <filename>\n";
    std::cout << "--output-file <filename>      " <<
        "Filename for the final vectors\n";
}

