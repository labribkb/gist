#include <tuple>
#include <string>
#include <vector>
#include <cassert>

// #include "SGD2.hpp"
#include "randomkit.h"
#include "mystructs.hpp"




void loadParams(Parametrizer &params, char *argv[]);

// Loading and saving

// Main functions, terms = delta ; adaptation of SGD2
// std::vector<term> layoutToTerms(Layout &g, Layout &init_g, Parametrizer &params);
// void OPTIMIZATION_PASS(Layout &g, Layout &init_g, std::vector<term> &terms, const std::vector<double> &etas, Parametrizer &params);
double circleDIJ(Coord const &ci, Coord const &cj, const Size &si, const Size &sj);
void fisheryates_shuffle(std::vector<term> &terms, rk_state &rstate);
std::vector<double> schedule(const std::vector<term> &terms, int t_max, double eps);

void saveImage(std::string savePath, std::vector<std::vector<int>> &v);
double do_gist(Layout & g,Parametrizer & params);
