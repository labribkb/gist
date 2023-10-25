#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <sstream>
#include <omp.h>
#include "GIST.hpp"

using namespace std;
vector<term> layoutToTerms(Layout &g, Layout &init_g, Parametrizer &params)
{
    vector<term> terms;
    bool overlap;
    double d_ij, w_ij;
    for (int i = 0; i < g.N(); ++i)
    {
        for (int j = i + 1; j < g.N(); ++j)
        {
            overlap = overlapCheck(g.nodes[i], g.nodes[j], params.tolerance);
            if (overlap)
            {
                d_ij = circleDIJ(g.nodes[i]->coord, g.nodes[j]->coord, g.nodes[i]->size, g.nodes[j]->size);
                w_ij = pow(d_ij, params.ALPHA * params.K);
            }
            else
            {
                d_ij = eucl(g.nodes[i]->coord, g.nodes[j]->coord);
                w_ij = pow(d_ij, params.ALPHA);
            }
            terms.push_back(term(i, j, d_ij, w_ij, overlap));
        }
    }
    return terms;
}

vector<term> overlapTerms(Layout &g, Layout &init_g, Parametrizer &params)
{
    vector<tuple<int, int>> overlaps;
    getAllOverlaps(g, overlaps, params.tolerance);

    vector<term> terms;
    double d_ij, w_ij;
    vector<vector<term>*> termsMultiThread;
    for(int i =0;i<omp_get_max_threads();++i)
        termsMultiThread.push_back(new vector<term>());
#pragma omp parallel for private(d_ij,w_ij)
    for (int cpt = 0; cpt < overlaps.size(); ++cpt)
    {
        auto [i, j] = overlaps[cpt];
        d_ij = circleDIJ(g.nodes[i]->coord, g.nodes[j]->coord, g.nodes[i]->size, g.nodes[j]->size);
        w_ij = pow(d_ij, params.ALPHA * params.K);
        termsMultiThread[omp_get_thread_num()]->push_back(term(i, j, d_ij, w_ij, true));
        // terms.push_back(term(i, j, d_ij, w_ij, true));
    }

    int size=0;
    for(unsigned int i =0;i<termsMultiThread.size();++i)
        size+=termsMultiThread[i]->size();
    terms.reserve(size);
    for(unsigned int i =0;i<termsMultiThread.size();++i){
        terms.insert(terms.end(),termsMultiThread[i]->begin(),termsMultiThread[i]->end());
        delete termsMultiThread[i];
    }
    return terms;
}

// S_GD2 optim algorithm, adapted from https://github.com/jxz12/s_gd2/blob/master/cpp/s_gd2/layout.cpp
void OPTIMIZATION_PASS(Layout &g, Layout &init_g, vector<term> &terms, const vector<double> &etas, Parametrizer &params)
{
    rk_state rstate;
    rk_seed(0, &rstate);
    double mvt_sum;
    unsigned int i_eta;
    for (i_eta = 0; i_eta < etas.size(); i_eta++)
    {
        const double eta = etas[i_eta];

        unsigned n_terms = terms.size();
        if (n_terms == 0)
            return;
        fisheryates_shuffle(terms, rstate);

        mvt_sum = 0;
        for (unsigned i_term = 0; i_term < n_terms; i_term++)
        {
            const term &t = terms[i_term];
            const int &i = t.i, &j = t.j;
            const double &w_ij = t.w;
            const double &d_ij = t.d;
            if (t.o)
            {
                // cap step size
                double mu = eta * w_ij;
                if (mu > 1)
                    mu = 1;

                double dx = g.nodes[i]->coord.x() - g.nodes[j]->coord.x();
                double dy = g.nodes[i]->coord.y() - g.nodes[j]->coord.y();
                if (dx == 0 && dy == 0)
                {
                    // jitter
                    dx = fRand(1e-10, 1e-9);
                    dy = fRand(1e-10, 1e-9);
                }
                double mag = sqrt(dx * dx + dy * dy);

                double r = 0;
                if (mag != 0)
                    r = (mu * (mag - d_ij)) / (2 * mag);
                double r_x = r * dx;
                double r_y = r * dy;
                mvt_sum += abs(r_x) + abs(r_y);

                g.nodes[i]->coord.set_x(g.nodes[i]->coord.x() - r_x);
                g.nodes[i]->coord.set_y(g.nodes[i]->coord.y() - r_y);

                g.nodes[j]->coord.set_x(g.nodes[j]->coord.x() + r_x);
                g.nodes[j]->coord.set_y(g.nodes[j]->coord.y() + r_y);
            }
        }
        if (mvt_sum < params.MINIMUM_MOVEMENT)
        {
            double o;
            return;
        }
        // terms = layoutToTerms(g, init_g, params);
        
        terms = overlapTerms(g, init_g, params);
    }

    return;
}

void passInOptim(Layout &g, Layout &init_g, Parametrizer &params)
{
    // vector<term> orig_terms = layoutToTerms(g, init_g, params);
    vector<term> orig_terms = overlapTerms(g, init_g, params);
    if (orig_terms.size() == 0){
        return;
    }
    vector<double> etas = schedule(orig_terms, params.MAX_ITER, params.eps);
    OPTIMIZATION_PASS(g, init_g, orig_terms, etas, params);
}

double getTolerance(Layout &g, double pix_tolerance, double R)
{
    auto [l, r, b, t] = g.getBB();
    double bb_w = r - l;
    double bb_h = t - b;
    double wD = max(bb_h, bb_w);

    double visual2geom = ((1 / R) * wD);
    double geom2visual = ((1 / wD) * R);
    double nodes_radius = g.nodes[0]->size.width() / 2;

    double nodes_radius_visual = nodes_radius * geom2visual;
    // cout << "nodesRadiusVisual : " << nodes_radius_visual << endl;

    double tolerance_geom = pix_tolerance * visual2geom;
    // cout << "tolerance visuel : " << tolerance_geom * geom2visual << endl;
    return tolerance_geom;

    // if (nodes_radius_visual - pix_tolerance > 1)
    //     return tolerance_geom;
    // else
    // {
    //     if (nodes_radius_visual < 1)
    //         tolerance_geom = 0;
    //     else
    //     {
    //         pix_tolerance = nodes_radius_visual - 1;
    //         tolerance_geom = pix_tolerance * visual2geom;
    //     }
    // }
    // cout << "pixTolerance " << pix_tolerance << endl;
    // return tolerance_geom;
}

double nodeRadiusVisual(Layout &g, Parametrizer &params)
{
    auto [l, r, b, t] = g.getBB();
    double bb_w = r - l;
    double bb_h = t - b;
    double wD = max(bb_h, bb_w);

    double geom2visual = ((1 / wD) * params.R);
    double nodes_radius = g.nodes[0]->size.width() / 2;

    double nodes_radius_visual = nodes_radius * geom2visual;
    return nodes_radius_visual;
}

tuple<double, double> screen_space_ratios(Layout &g, int R)
{
    auto [l, r, b, t] = g.getBB();
    double bb_w = r - l;
    double bb_h = t - b;
    double wD, mD;
    if (bb_h > bb_w)
    {
        wD = bb_h;
        mD = b;
    }
    else
    {
        wD = bb_w;
        mD = l;
    }
    // cout << " R : " << R << " ; wD : " << wD << endl;
    double visual2geom = ((1. / R) * wD);
    double geom2visual = ((1. / wD) * R);
    // cout << " visual2geom : " << visual2geom << " ; geom2visual : " << geom2visual << endl;
    return make_tuple(visual2geom, geom2visual);
}

void setDiam(Layout &g, double curDiam)
{
    for (int i = 0; i < g.nodes.size(); ++i)
    {
        g.nodes[i]->size.a = curDiam;
        g.nodes[i]->size.b = curDiam;
    }
}

bool full_pass_gist(Layout &g, Layout &init_g, Parametrizer &params, double pix_tolerance, Layout &best,double & best_diam,double & best_tol){
    auto [visual2geom, geom2visual] = screen_space_ratios(g, params.R);
    double nodes_diameter = g.nodes[0]->size.width();
    double min_diameter = 1;
    double max_diameter = nodes_diameter * geom2visual+1;
    const double init_pix_tolerance = pix_tolerance;

    while(max_diameter < min_diameter){
        params.R *=1.5;
        auto [v2g0, g2v0] = screen_space_ratios(g, params.R);
        max_diameter = nodes_diameter * g2v0;
    }

    // cout << "Upper diameter bound : " << max_diameter << endl;
    // cout << "lower diameter bound : " << min_diameter << endl;
    // cout << "R = " << params.R << endl;
    double upperDiam = max_diameter;
    double lowerDiam = min_diameter;
    double oldDiam = min_diameter;
    double curDiam = min_diameter;
    bool foundAtLeastOneLayout = false;
    bool lastWasWithoutOverlaps = false;
    int n_passes = 0;
    double scaleFactor;

    // double wanted_pix_tol = max(min(pix_tolerance, max_diameter-1), 0.);
    pix_tolerance = init_pix_tolerance;
    double wanted_pix_tol = (max_diameter-2*pix_tolerance) >= 1 ? pix_tolerance : 0;
    while(pix_tolerance > 1 && wanted_pix_tol == 0){
        pix_tolerance-=0.1;
        wanted_pix_tol = (max_diameter-2*pix_tolerance) >= 1 ? pix_tolerance : 0;
    }
    params.tolerance = wanted_pix_tol * visual2geom;
    cout << "tolerance : " << params.tolerance * geom2visual << endl;
    bool stop = !scanLineOverlapCheck(g, params.tolerance); // do not enter if there is no overlap
    
    if (!stop && g.isCurrentScaleSolvable(params.tolerance))
    {
        passInOptim(g, init_g, params);
        stop = !scanLineOverlapCheck(g, params.tolerance);
    }
    foundAtLeastOneLayout = stop;
    lastWasWithoutOverlaps = stop;
    // vector<tuple<int, int>> overlaps;
    best = Layout(g);
    best_diam = nodes_diameter *geom2visual;
    best_tol = wanted_pix_tol;

    while (!stop)
    {
        curDiam = (upperDiam + lowerDiam) / 2;
        cout << "==== pass " << n_passes << " ; diam : " << curDiam << " ====" << endl;
        scaleFactor = curDiam / oldDiam;
        oldDiam = curDiam;
        auto [v2g, g2v] = screen_space_ratios(g, params.R);
        /*loann-fred*/
        g = Layout(best);
        setDiam(g, v2g * curDiam);
        // wanted_pix_tol = min(pix_tolerance, curDiam-1);
        // params.tolerance = pix_tolerance < curDiam / 3 ? pix_tolerance * v2g : 0;
        // wanted_pix_tol = max(min(pix_tolerance, curDiam-1), 0.);
        pix_tolerance = init_pix_tolerance;
        wanted_pix_tol = (curDiam-2*pix_tolerance) >= 1 ? pix_tolerance : 0;
        while(pix_tolerance > 1 && wanted_pix_tol == 0){
            pix_tolerance-=0.1;
            wanted_pix_tol = (curDiam-2*pix_tolerance) >= 1 ? pix_tolerance : 0;
        }
        params.tolerance = wanted_pix_tol * visual2geom;
        // wanted_pix_tol = (curDiam-2*pix_tolerance) >= 1 ? pix_tolerance : 0;
        // params.tolerance = wanted_pix_tol * visual2geom;
        // cout << "vis tolerance : " << params.tolerance * g2v << endl;
        // cout << "geom tol : " << params.tolerance << endl;
        // cout << "geom diam : " << v2g*curDiam << endl;

        passInOptim(g, init_g, params);

        if (scanLineOverlapCheck(g, params.tolerance))
        {
            lastWasWithoutOverlaps = false;
            upperDiam = curDiam;
            cout << "KO" << endl;
        }
        else // no overlap
        {
            cout << "OK" << endl;

            lastWasWithoutOverlaps = true;
            lowerDiam = curDiam;
            best = Layout(g);
            best_diam = curDiam;
            best_tol=wanted_pix_tol;
            foundAtLeastOneLayout = true;
        }
        if (upperDiam - lowerDiam < params.SCALE_STEP)
        {
            stop = true;
        }
        ++n_passes;
        if (n_passes >= params.MAX_PASSES)
            break;

    }
    return foundAtLeastOneLayout;
}
double do_gist(Layout & input,Parametrizer & params){
    Layout g(input);
    int n_passes = 0;
    double scaleFactor;
    auto start = std::chrono::steady_clock::now();
    Layout init_g(g);
    Layout best(g);
    double init_pix_tolerance = params.tolerance;
    const double rescaleFactor = 2.;
    int res_iter=0;
    double best_diam,best_tol;
    while(!full_pass_gist(g, init_g, params, init_pix_tolerance, best,best_diam,best_tol)){
        params.R *= rescaleFactor;
        cout << "R changed to " << params.R;
        g = Layout(init_g);
        // g.scale(rescaleFactor);
    }
    for(int i =0;i<input.nodes.size();++i){
	    input.nodes[i]->coord.set_x(best.nodes[i]->coord.x());
	    input.nodes[i]->coord.set_y(best.nodes[i]->coord.y());
	    input.nodes[i]->size.set(best.nodes[i]->size.a,best.nodes[i]->size.b);
    }
    auto end = std::chrono::steady_clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count();;


}
int main(int argc, char *argv[])
{
    if (argc != 12)
    {
        cerr << "Wrong parameters" <<argc << endl;
        return EXIT_FAILURE;
    }
    string input_layout_path = argv[1];
    Layout g;
    Parametrizer params;
    loadParams(params, argv);
    g.loadLayout(input_layout_path);
    double elapsed = do_gist(g,params);
    string savepath = g.save(input_layout_path, "gist", elapsed, params.PRIME);

    cout << "DONE in " << elapsed << " ms; saved to " << savepath << endl;

    return EXIT_SUCCESS;
}