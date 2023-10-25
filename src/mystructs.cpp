#include <iostream>
#include <fstream>
#include <float.h>
#include <sstream>
#include <cmath>

#include "mystructs.hpp"

using namespace std;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

Node::Node(const Coord &c, const Size &s, int id)
{
    this->coord = Coord(c);
    this->size = Size(s);
    this->id = id;
}

NodesComposite::NodesComposite(vector<INode *> &nodes)
{
    this->nodes = nodes;
}

NodesComposite::NodesComposite(const Coord &c, const Size &s)
{
    this->coord = Coord(c);
    this->size = Size(s);
}

void Node::move(double mvt_x, double mvt_y)
{
    this->coord.set_x(this->coord.x() + mvt_x);
    this->coord.set_y(this->coord.y() + mvt_y);
}

void NodesComposite::move(double mvt_x, double mvt_y)
{
    this->coord.set_x(this->coord.x() + mvt_x);
    this->coord.set_y(this->coord.y() + mvt_y);
    for (int i = 0; i < this->nodes.size(); ++i)
        this->nodes[i]->move(mvt_x, mvt_y);
}

NodesComposite::NodesComposite(const NodesComposite &other)
{
    this->coord = Coord(other.coord);
    this->size = Size(other.size);
    this->nodes = vector<INode *>(other.nodes);
}

Layout::~Layout()
{
    // cout << "LAYOUT DESTROYER" << endl;
    for (int i = 0; i < this->N(); ++i)
        delete this->nodes[i];
}

INode *NodesComposite::remove(int i)
{
    INode *ret;
    if (i > 0 && i < this->nodes.size())
    {
        ret = this->nodes.at(i);
        this->nodes.erase(this->nodes.begin() + i);
    }
    return ret;
}

Layout::Layout(vector<Coord> &pos, vector<Size> &sizes)
{
    for (int i = 0; i < pos.size(); ++i)
        this->nodes.push_back(new Node(pos[i], sizes[i], i));
}

Layout::Layout(const Layout &other)
{
    this->deepcopy(other);
}
void Layout::deepcopy(const Layout &other)
{
    this->nodes.clear();
    this->nodes.reserve(other.nodes.size());
    for (int i = 0; i < other.nodes.size(); ++i)
    {
        this->nodes.push_back(new Node(Coord(other.nodes[i]->coord.x(), other.nodes[i]->coord.y()), Size(other.nodes[i]->size.width(), other.nodes[i]->size.height()), i));
    }
}
Layout &Layout::operator=(const Layout &other)
{
    this->deepcopy(other);
    return *this;
}
void NodesComposite::scale(double scaleFactor)
{
    for (int i = 0; i < this->nodes.size(); ++i)
    {
        this->nodes[i]->coord.set(this->nodes[i]->coord.x() * scaleFactor, this->nodes[i]->coord.y() * scaleFactor);
    }
}

void Node::print()
{
    cout << "ID = " << id << "\tX: " << this->coord.x() << " ; Y: " << this->coord.y() << "\t\t"
         << "W: " << this->size.width() << " ; H: " << this->size.height() << endl;
}

void Node::print(string prefix)
{
    cout << prefix;
    this->print();
}

void NodesComposite::print()
{
    cout << "X: " << this->coord.x() << " ; Y: " << this->coord.y() << "\t\t "
         << "W: " << this->size.width() << " ; H: " << this->size.height() << "\t\t " << endl;
    for (int i = 0; i < this->N(); i++)
    {
        cout << "\t i: " << i << "\t\t";
        this->nodes[i]->print();
    }
    cout << "\n"
         << endl;
}

void NodesComposite::print(string prefix)
{
    cout << prefix << "X: " << this->coord.x() << " ; Y: " << this->coord.y() << "\t\t "
         << "W: " << this->size.width() << " ; H: " << this->size.height() << "\t\t " << endl;
    for (int i = 0; i < this->N(); i++)
    {
        cout << prefix << "i: " << i << "\t";
        this->nodes[i]->print(prefix + "\t\t");
    }
}

void NodesComposite::print(int n)
{
    cout << "X: " << this->coord.x() << " ; Y: " << this->coord.y() << "\t\t "
         << "W: " << this->size.width() << " ; H: " << this->size.height() << "\t\t " << endl;

    int bound = min(n, this->N());
    for (int i = 0; i < bound; i++)
    {
        cout << "\t i: " << i << "\t\t";
        this->nodes[i]->print();
    }
}

// left, right, bot, top
tuple<double, double, double, double> Node::getBB()const
{
    return make_tuple(
        this->coord.x() - this->size.width() / 2,
        this->coord.x() + this->size.width() / 2,
        this->coord.y() - this->size.height() / 2,
        this->coord.y() + this->size.height() / 2);
}

// left, right, bot, top
tuple<double, double, double, double> NodesComposite::getBB()const
{
    double leftmost = DBL_MAX;
    double rightmost = -DBL_MAX;
    double topmost = -DBL_MAX;
    double botmost = DBL_MAX;

    for (int i = 0; i < this->N(); ++i)
    {
        auto [left, right, bot, top] = this->nodes[i]->getBB();
        leftmost = min(leftmost, left);
        rightmost = max(rightmost, right);
        botmost = min(botmost, bot);
        topmost = max(topmost, top);
    }
    return make_tuple(leftmost, rightmost, botmost, topmost);
}

int NodesComposite::N(bool deep)
{
    if (!deep)
        return this->N();

    int n = 0;
    for (int i = 0; i < this->nodes.size(); ++i)
    {
        n += this->nodes[i]->N(true);
    }
    return n;
}

double Layout::maxScaleRatio(double tolerance)
{
    double padding = 1e-4;
    double maxRatio = 1.;
    double optimalDist, actualDist, ratio, unoverlapRatio;
    double actualX, actualY, desiredWidth, desiredHeight, widthRatio, heightRatio;
    vector<tuple<int, int>> overlaps;
    getAllOverlaps(*this, overlaps, tolerance);
    INode *nu, *nv;
    for (unsigned int i = 0; i < overlaps.size(); ++i)
    {
        auto [u, v] = overlaps[i];
        nu = this->nodes[u];
        nv = this->nodes[v];
        actualDist = eucl(nu->coord.x(), nu->coord.y(), nv->coord.x(), nv->coord.y());

        actualX = abs(nu->coord.x() - nv->coord.x());
        actualY = abs(nu->coord.y() - nv->coord.y());
        if(actualX == 0 && actualY == 0){
            // jitter
            double r1 = nu->size.width()/2;
            actualX = fRand(r1/100, r1/10); //fRand(1e-10, 1e-9);
            actualY = fRand(r1/100, r1/10);//fRand(1e-10, 1e-9);
            nu->coord.a += actualX;
            nu->coord.b += actualY;
            actualDist = sqrt(actualX*actualX + actualY*actualY);
        }
        actualDist += tolerance;
        desiredWidth = (nu->size.width() + nv->size.width()) / 2 + padding;//-tolerance;
        desiredHeight = (nu->size.height() + nv->size.height()) / 2 + padding;//-tolerance;

        widthRatio = desiredWidth / actualX;
        heightRatio = desiredHeight / actualY;

        unoverlapRatio = min(widthRatio, heightRatio);
        actualX *= unoverlapRatio;
        actualY *= unoverlapRatio;

        optimalDist = vecNorm2D(actualX, actualY);
        ratio = optimalDist / actualDist;
        maxRatio = max(maxRatio, ratio);
    }
    return maxRatio;
}

bool Layout::isCurrentScaleSolvable(double tolerance)
{
    double areas_sum = 0;
    double min_x = DBL_MAX;
    double min_y = DBL_MAX;
    double max_x = -DBL_MAX;
    double max_y = -DBL_MAX;

    double left, right, top, bot;
    INode *n;
    double minRadius=DBL_MAX;
    double maxRadius=-DBL_MAX;
    for (int i = 0; i < this->nodes.size(); ++i)
    {
        n = this->nodes[i];
        //areas_sum += n->size.width() * n->size.height();
        areas_sum += n->size.width() * n->size.height()*M_PI/4;
        if(minRadius<n->size.width()/2)
            minRadius=n->size.width()/2;
        if(maxRadius>n->size.width()/2)
            maxRadius=n->size.width()/2;
        left = n->coord.x() - n->size.width() / 2;
        right = n->coord.x() + n->size.width() / 2;
        top = n->coord.y() + n->size.height() / 2;
        bot = n->coord.y() - n->size.height() / 2;
        if (left < min_x)
            min_x = left;
        if (right > max_x)
            max_x = right;
        if (bot < min_y)
            min_y = bot;
        if (top > max_y)
            max_y = top;
    }
    double bb_area = (max_x - min_x) * (max_y - min_y);

    //if(tolerance > 0){
        //double packing_diameter = 2*minRadius-tolerance;
        //if(packing_diameter>0){
            //double packing_area_sum = this->nodes.size()*M_PI*packing_diameter/4;
            //return packing_area_sum>=bb_area;
        //}
        return bb_area >= areas_sum;
    //}else{
        //return bb_area >= areas_sum;
    //}
}
void Layout::read(string s){
    stringstream sstream;
    sstream<<s;
    nodes.clear();
    loadLayout(sstream,false);
}
string Layout::write(const std::string & sep,const std::string & linesep ){
    stringstream sstream;
    //string sep = " ";
    sstream<<this->nodes.size()<<linesep;
    sstream.precision(20);
    INode *n;
    for (unsigned int i = 0; i < this->nodes.size(); ++i)
    {
        n = this->nodes[i];
        sstream << n->coord.x() << sep << n->coord.y() << sep << n->size.width() << sep << n->size.height() << linesep;
    }
    return sstream.str();
}
string Layout::save(string path, string extension, double elapsed, bool isPrime)
{
    string sep = " ";
    if (isPrime)
        extension += "p";
    ofstream f(path + "." + extension);
    f << elapsed << endl;
    INode *n;
    for (unsigned int i = 0; i < this->nodes.size(); ++i)
    {
        n = this->nodes[i];
        f << n->coord.x() << sep << n->coord.y() << sep << n->size.width() << sep << n->size.height() << endl;
    }
    f.close();
    cout << "save : " << path + "." + extension << endl;
    return path + "." + extension;
}

void Layout::loadLayout(string input_path){
    ifstream f;
    f.open(input_path);
    loadLayout(f);
    f.close();
}
std::vector<int> Layout::loadLayout(string input_path,bool hasLabels){
    ifstream f;
    f.open(input_path);
    std::vector<int> labels = loadLayout(f,hasLabels);
    f.close();
    return labels;
}
std::vector<int> Layout::loadLayout(istream & f,bool hasLabels){
    int N;
    f >> N;

    vector<Coord> pos;
    pos.reserve(N);
    vector<Size> sizes;
    sizes.reserve(N);
    vector<int> labels;
    labels.reserve(N);

    double x, y, w, h;
    double l;
    for (int i = 0; i < N; i++)
    {
        if(hasLabels){
            f >> x >> y >> w >> h >> l;
            labels.push_back(l);
        }
        else
            f >> x >> y >> w >> h;
        pos.emplace_back(Coord(x, y));
        sizes.emplace_back(Size(w, h));
    }

    this->nodes.clear();
    this->nodes.resize(pos.size());
    for (int i = 0; i < pos.size(); ++i)
        this->nodes[i] = new Node(pos[i], sizes[i], i);
    return labels;
}
void Layout::loadLayout(istream & f)
{
    int N;
    f >> N;

    vector<Coord> pos;
    pos.reserve(N);
    vector<Size> sizes;
    sizes.reserve(N);

    double x, y, w, h;
    string firstline;
    f>> x >> y >> w >> h;
    std::getline(f,firstline);
    bool label = firstline.size()>0;
    
    pos.emplace_back(Coord(x, y));
    sizes.emplace_back(Size(w, h));
    double l;
    for (int i = 1; i < N; i++)
    {
        if(label)
            f >> x >> y >> w >> h >> l;
        else
            f >> x >> y >> w >> h;
        pos.emplace_back(Coord(x, y));
        sizes.emplace_back(Size(w, h));
    }

    this->nodes.clear();
    this->nodes.resize(pos.size());
    for (int i = 0; i < pos.size(); ++i)
        this->nodes[i] = new Node(pos[i], sizes[i], i);
}


void NodesComposite::add(INode *n)
{
    this->nodes.push_back(n);
}
double Layout::minScale()const{
    auto [l,r,b,t] = getBB();
    double area = 0.0;
    for(int i=0;i<nodes.size();++i){
        area+=M_PI*nodes[i]->size.width()*nodes[i]->size.width()/4;
    }
    double bbArea = (r-l)*(t-b);
    if(area<bbArea)
        return 1.0;
    else
        return sqrt(area/bbArea);
}
