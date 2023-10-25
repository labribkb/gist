#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <numeric>   // std::iota
#include <algorithm> // std::sort, std::stable_sort
#include <float.h>

#include "GIST.hpp"

using namespace std;


double SIDE2SIDE_delta(Coord const &ci, Coord const &cj, const Size &si, const Size &sj, double intersec_width, double intersec_height)
{
    double d_eucl = eucl(ci, cj);
    return d_eucl + min(abs(intersec_width), abs(intersec_height));
}

double PRISM_delta(Coord const &ci, Coord const &cj, const Size &si, const Size &sj, double intersec_width, double intersec_height)
{
    double t_ij = max(
        1.,
        min((si.width() / 2 + sj.width() / 2) / abs(ci.x() - cj.x()),
            (si.height() / 2 + sj.height() / 2) / abs(ci.y() - cj.y())));

    return t_ij * eucl(ci, cj);
}

double circleDIJ(Coord const &ci, Coord const &cj, const Size &si, const Size &sj)
{
    double r1 = si.width()/2;
    double r2 = sj.width()/2;
    return r1+r2;
}
void loadParams(Parametrizer &params, char *argv[])
{
    params.ALPHA = atof(argv[2]);
    params.K = atof(argv[3]);
    params.MINIMUM_MOVEMENT = atof(argv[4]);
    params.MAX_ITER = atoi(argv[5]);
    params.MAX_PASSES = atoi(argv[6]);
    params.SCALE_STEP = atof(argv[7]);
    params.PRIME = atoi(argv[8]);
    params.tolerance = atof(argv[9]);
    params.maxScale = atof(argv[10]);
    params.R = atoi(argv[11]);

}

// S_GD2 function, taken from https://github.com/jxz12/s_gd2/blob/master/cpp/s_gd2/layout.cpp
vector<double> schedule(const vector<term> &terms, int t_max, double eps)
{
    double w_min = terms[0].w, w_max = terms[0].w;
    for (unsigned i = 1; i < terms.size(); i++)
    {
        const double &w = terms[i].w;
        if (w < w_min)
            w_min = w;
        if (w > w_max)
            w_max = w;
    }
    double eta_max = 1.0 / w_min;
    double eta_min = eps / w_max;
    double lambda = log(eta_max / eta_min) / ((double)t_max - 1);

    // initialize step sizes
    vector<double> etas;
    etas.reserve(t_max);
    for (int t = 0; t < t_max; t++)
    {
        etas.push_back(eta_max * exp(-lambda * t));
        // cout << eta_max * exp(-lambda * t) << endl;
    }
    return etas;
}

void fisheryates_shuffle(vector<term> &terms, rk_state &rstate)
{
    int n = terms.size();
    for (unsigned i = n - 1; i >= 1; i--)
    {
        unsigned j = rk_interval(i, &rstate);
        term temp = terms[i];
        terms[i] = terms[j];
        terms[j] = temp;
    }
}

void saveImage(string savePath, vector<vector<int>> &v)
{
    ofstream f(savePath);
    int height = v.size();
    int width = v[0].size();

    f << width << "," << height << endl;

    for (int h = 0; h < v.size(); ++h)
    {
        for (int r = 0; r < v[h].size(); ++r)
        {
            f << v[h][r];
            if(r == v[h].size()-1)
                f << endl;
            else
                f << ",";
        }
    }
    // f << n->coord.x() << sep << n->coord.y() << sep << n->size.width() << sep << n->size.height() << endl;
    f.close();
}

tuple<double, double, double, double> getBB(vector<INode *> &superNodes)
{
    double minX, maxX, minY, maxY;
    maxY = maxX = -DBL_MAX;
    minY = minX = DBL_MAX;
    for (int i = 0; i < superNodes.size(); ++i)
    {
        auto [left, right, bot, top] = superNodes[i]->getBB();
        minX = min(minX, left);
        maxX = max(maxX, right);
        minY = min(minY, bot);
        maxY = max(maxY, top);
    }
    return make_tuple(minX, maxX, minY, maxY);
}
