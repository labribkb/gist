#include <tuple>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <unordered_set>


#include "mystructs.hpp"

using namespace std;

double maxOfLeft(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
{
    double left1 = p1.x() - s1.width() / 2;
    double left2 = p2.x() - s2.width() / 2;
    return left1 > left2 ? left1 : left2;
}

double minOfRight(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
{
    double right1 = p1.x() + s1.width() / 2;
    double right2 = p2.x() + s2.width() / 2;
    return right1 < right2 ? right1 : right2;
}

double minOfTop(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
{
    double top1 = p1.y() + s1.height() / 2;
    double top2 = p2.y() + s2.height() / 2;
    return top1 < top2 ? top1 : top2;
}

double maxOfBot(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
{
    double bot1 = p1.y() - s1.height() / 2;
    double bot2 = p2.y() - s2.height() / 2;
    return bot1 > bot2 ? bot1 : bot2;
}

tuple<double, double> nodeRectanglesIntersection(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
{
    double min_of_right = minOfRight(p1, s1, p2, s2);
    double max_of_left = maxOfLeft(p1, s1, p2, s2);
    double min_of_top = minOfTop(p1, s1, p2, s2);
    double max_of_bot = maxOfBot(p1, s1, p2, s2);

    double intersection_width = min_of_right - max_of_left;
    double intersection_height = min_of_top - max_of_bot;
    return make_tuple(intersection_width, intersection_height);
}

bool contains(Coord &p, Coord &p2, Size &s2)
{
    if (p.x() >= p2.x() - s2.width() / 2 && p.x() <= p2.x() + s2.width() / 2 && p.y() >= p2.y() - s2.height() / 2 && p.y() >= p2.y() + s2.height() / 2)
        return true;
    return false;
}


vector<size_t> sort_indexes(const vector<double> &v)
{

    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    stable_sort(
            idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2)
                { return v[i1] < v[i2]; });

    return idx;
}

double vecNorm2D(double vec_x, double vec_y)
{
    return sqrt(vec_x * vec_x + vec_y * vec_y);
}

double eucl(double x1, double y1, double x2, double y2)
{
    return sqrt(pow(x1 - x2, 2.) + pow(y1 - y2, 2.));
}

double eucl(Coord const &p1, Coord const &p2)
{
    return sqrt(pow(p1.x() - p2.x(), 2.) + pow(p1.y() - p2.y(), 2.));
}

vector<size_t> sortNodesByX(vector<INode *> &nodes_group)
{
    int N = nodes_group.size();
    vector<INode *> tmp(nodes_group);

    vector<double> xs(N);
    for (int i = 0; i < N; ++i){
        xs[i] = nodes_group[i]->coord.x() - nodes_group[i]->size.width() /2;
    }

    vector<size_t> idx = sort_indexes(xs);
    
    for (int i = 0; i < N; ++i)
        nodes_group[i] = tmp[idx[i]];
    return idx;
}

std::vector<size_t> sort_indexes(const double *v, int size)
{

    return sort_indexes(std::vector<double>(v,v+size));
}

bool overlapCheck(INode *ni, INode *nj, double tolerance)
{
    double o;
    return overlapCheck(ni,nj,tolerance,o);
}
bool overlapCheck(INode * ni, INode * nj, double tolerance,double & overlap)
{
    double dist = eucl(ni->coord, nj->coord);
    double r1 = ni->size.width()/2;
    double r2 = nj->size.width()/2;
    if(dist > (r1+r2)-max(0.0,tolerance)){
        overlap=max(0.0,r1+r2-dist);
        return false;
    }
    return true;
}

bool scanLineOverlapCheck(NodesComposite &nodes_group, double tolerance){
    double o;
    return scanLineOverlapCheck(nodes_group,tolerance,o);
}
bool scanLineOverlapCheck(NodesComposite &nodes_group, double tolerance,double &maxOverlap)
{
    maxOverlap=0.0;
    double o;
    NodesComposite tmp_nodes_group(nodes_group);
    int N = tmp_nodes_group.nodes.size();
    sortNodesByX(tmp_nodes_group.nodes);

    INode* ni;
    INode* nj;
    double left_j, right_i;
    for (int i = 0; i < N; ++i)
    {
        ni = tmp_nodes_group.nodes[i];
        right_i = ni->coord.x() + ni->size.width() / 2;
        for (int j = i + 1; j < N; ++j)
        {
            nj = tmp_nodes_group.nodes[j];
            left_j = nj->coord.x() - nj->size.width() / 2;
            if (overlapCheck(ni, nj, tolerance,o)){
                return true;
}
            else {
                maxOverlap=max(maxOverlap,o);
                if (left_j > right_i)
                    break;
            }
        }
    }
    return false;
}

bool scanLineOverlapCheck(vector<INode *> &nodes_group, double tolerance)
{
    vector<INode *> tmp_nodes_group(nodes_group);
    int N = tmp_nodes_group.size();
    sortNodesByX(tmp_nodes_group);

    INode* ni;
    INode* nj;
    double left_j, right_i;
    for (int i = 0; i < N; ++i)
    {
        ni = tmp_nodes_group[i];
        right_i = ni->coord.x() + ni->size.width() / 2;
        for (int j = i + 1; j < N; ++j)
        {
            nj = tmp_nodes_group[j];
            left_j = nj->coord.x() - nj->size.width() / 2;
            if (overlapCheck(ni, nj, tolerance))
                return true;
            else if (left_j > right_i)
                break;
        }
    }
    return false;
}

// want a copy because we are going to sort it and we don't want to remember the sorted order
void getAllOverlaps(NodesComposite &nodes_group, vector<tuple<int, int>> &overlaps, double tolerance) 
{
    INode* ni;
    INode* nj;
    double left_j, right_i;
    int N = nodes_group.nodes.size();
    NodesComposite tmp_nodes_group(nodes_group);
    vector<size_t> mapping = sortNodesByX(tmp_nodes_group.nodes);

    vector<vector<tuple<int,int>>*> overlapsMT;
    for(unsigned int i =0;i<omp_get_max_threads();++i)
        overlapsMT.push_back(new vector<tuple<int,int>>());
#pragma omp parallel for  private(ni,nj,left_j,right_i)
    for (int i = 0; i < N; ++i)
    {
        ni = tmp_nodes_group.nodes[i];
        right_i = ni->coord.x() + ni->size.width() / 2;
        for (int j = i + 1; j < N; ++j)
        {
            nj = tmp_nodes_group.nodes[j];
            left_j = nj->coord.x() - nj->size.width() / 2;
            if (overlapCheck(ni, nj, tolerance)){
                overlapsMT[omp_get_thread_num()]->emplace_back(make_tuple(mapping[i], mapping[j]));
                // overlaps.emplace_back(make_tuple(mapping[i], mapping[j]));
            }
            else if (left_j > right_i)
                break;
        }
    }
    
    int size=0;
    for(unsigned int i =0;i<overlapsMT.size();++i){
        size+=overlapsMT[i]->size();
    }
    overlaps.reserve(size);
    for(unsigned int i =0;i<overlapsMT.size();++i){
        overlaps.insert(overlaps.end(),overlapsMT[i]->begin(),overlapsMT[i]->end());
        delete overlapsMT[i];
    }

}

tuple<float, float, float, float> getSuperShapeBB(INode *n1, INode *n2)
{
    auto [l1, r1, b1, t1] = n1->getBB();
    auto [l2, r2, b2, t2] = n2->getBB();

    return make_tuple(
        min(l1, l2),
        max(r1, r2),
        min(b1, b2),
        max(t1, t2));
}
