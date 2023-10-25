#ifndef LAYOUT_H
#define LAYOUT_H

#include <vector>
#include <string>
#include <tuple>

struct Parametrizer
{
    double K = -2;
    double ALPHA = 1;
    double MINIMUM_MOVEMENT = 1e-4;
    int MAX_ITER = 200;
    int MAX_PASSES = 1000;
    double SCALE_STEP = 0.1;
    double eps = 0.01;
    double tolerance=1;
    double maxScale=-1;
    bool PRIME = false;
    int R=2000;
};

struct Vec2d
{
    double a;
    double b;
    Vec2d(double x, double y) : a(x), b(y) {}
    Vec2d() = default;

    inline void set(double a, double b)
    {
        this->a = a;
        this->b = b;
    }
};

struct Coord : Vec2d
{
    Coord(const Coord &c) : Vec2d(c.x(), c.y()) {}
    Coord(double x, double y) : Vec2d(x, y) {}
    Coord() = default;

    inline double x() const { return a; }
    inline double y() const { return b; }
    inline void set_x(double x) { a = x; }
    inline void set_y(double y) { b = y; }
};

struct Size : Vec2d
{
    Size(const Size &s) : Vec2d(s.width(), s.height()) {}
    Size(double width, double height) : Vec2d(width, height) {}
    Size() = default;

    inline double width() const { return this->a; }
    inline double height() const { return this->b; }
    inline void set_width(double x) { a = x; }
    inline void set_height(double y) { b = y; }
};

struct SizeCircle : Size
{
    SizeCircle(double radius) : Size(radius, -666){}
    SizeCircle() = default;

    inline double radius() const {return this->a;}
    inline double diameter() const {return this->a*2;}
    inline double width() const { return diameter(); }
    inline double height() const { return diameter(); }
    
    inline void set_radius(double r){this->a = r;}
};

struct term
{
    int i, j;
    double d, w;
    bool o;
    bool f;
    term(int i, int j, double d, double w, bool o, bool f) : i(i), j(j), d(d), w(w), o(o), f(f) {}
    term(int i, int j, double d, double w, bool o) : i(i), j(j), d(d), w(w), o(o) {}
    term(int i, int j, double d, double w) : i(i), j(j), d(d), w(w) {}
};

class INode
{
public:
    Coord coord;
    Size size;
    INode *parent;
    virtual ~INode() {}
    virtual void move(double mvt_x, double mvt_y) = 0;
    virtual void print() = 0;
    virtual void print(std::string prefix) = 0;
    virtual int N() { return -1; }
    virtual int depth() =0;
    virtual int N(bool deep) {return 1;}
    virtual std::tuple<double, double, double, double> getBB()const = 0;
};

class Node : public INode
{
private:
public:
    int id;
    Node(const Coord &c, const Size &s, int id);
    Node() = default;
    ~Node() = default;
    void print();
    void print(std::string prefix);
    void move(double mvt_x, double mvt_y);
    int depth() {return 1;}
    std::tuple<double, double, double, double> getBB()const;
};

class NodesComposite : public INode
{
public:
    std::vector<INode *> nodes;
    NodesComposite(const Coord &c, const Size &s);
    NodesComposite() = default;
    NodesComposite(std::vector<INode *> &nodes);
    NodesComposite(const NodesComposite &other);
    void move(double mvt_x, double mvt_y);
    virtual void add(INode *n);
    void scale(double scaleFactor);
    int N() const{ return nodes.size(); }
    int N(bool deep);
    void print();
    void print(int n);
    int depth(){return this->N() > 0 ? 1+this->nodes[0]->depth() : 1;};
    std::tuple<double, double, double, double> getBB()const;
    void print(std::string prefix);
    INode *remove(int i);
};

class Layout : public NodesComposite
{
public:
    Layout(std::vector<Coord> &pos, std::vector<Size> &sizes);
    Layout(const Layout &other);
    Layout() = default;
    ~Layout();
    void deepcopy(const Layout &other);
    Layout & operator=(const Layout &other);
    void loadLayout(std::string input_path);
    std::vector<int> loadLayout(std::string input_path,bool hasLabels);
    void loadLayout(std::istream & f);
    std::vector<int> loadLayout(std::istream & f,bool hasLabels);
    double maxScaleRatio(double tolerance);
    void move(double mvt_x, double mvt_y){}; // interdire
    bool isCurrentScaleSolvable(double tolerance);
//lower bound of solvable scale (might not be solvable)
    double minScale()const;
    std::string save(std::string path, std::string extension, double elapsed, bool isPrime);
    void read(std::string s);
    std::string write(const std::string & sep=" ",const std::string & linesep=" ");
};

double maxOfLeft(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);
double minOfRight(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);
double minOfTop(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);
double maxOfBot(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);
std::tuple<double, double> nodeRectanglesIntersection(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2);


bool contains(Coord &p, Coord &p2, Size &s2);
double vecNorm2D(double vec_x, double vec_y);

double eucl(double x1, double y1, double x2, double y2);
double eucl(Coord const &p1, Coord const &p2);

// Getting overlaps
std::vector<size_t> sort_indexes(const double *v, int size);
bool overlapCheck(INode *ni, INode *nj, double tolerance);
bool overlapCheck(INode *ni, INode *nj, double tolerance,double &m);
bool scanLineOverlapCheck(NodesComposite &nodes_group, double tolerance,double & m);
bool scanLineOverlapCheck(NodesComposite &nodes_group, double tolerance);
bool scanLineOverlapCheck(std::vector<INode *> &nodes_group, double tolerance);
void getAllOverlaps(NodesComposite &nodes_group, std::vector<std::tuple<int, int>> &overlaps, double tolerance);
std::tuple<double, double, double, double> getBB(std::vector<INode *> &superNodes);

std::tuple<float, float, float, float> getSuperShapeBB(INode *n1, INode *n2);
double fRand(double fMin, double fMax);
#endif
