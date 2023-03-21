#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

double norm(std::vector<int> vec);
std::vector<double> computeNorms(std::vector<std::vector<int>> vecs);
std::vector<int> matMul(std::vector<int> mat, std::vector<int> vec);
int writeAllVecs(std::string filename, const std::vector<std::vector<int>> & vectors);
std::vector<std::vector<int>> generateCandidates(int gridMax, int dim);
std::vector<std::vector<int>> selectVectors(std::vector<std::vector<int>> candidateQs, int targetLength);
std::vector<std::vector<int>> generatePermMatrices(int dim);
std::vector<std::vector<int>> includePermutations(std::vector<std::vector<int>> selectedQs);
void showhelp();
