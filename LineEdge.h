#ifndef LINEEDGE_H
#define LINEEDGE_H
#include <vector>
#include <cmath>

#define UNIFORM -1
#define EXPREFINE0 0 
#define EXPREFINE1 1
#define SINREFINE0 2
#define SINREFINE1 3
#define SINREFINE2 4
#define QUDREFINE0 5
#define QUDREFINE1 6

class LineEdge {
public: 
    LineEdge(double* p0, double* p1, int N, int refineType, double h0, double h1);
    std::vector<double> Evaluate(double s);
    int m_N;
private:
    double* m_p0;
    double* m_p1;
    int m_refineType;
    double m_g0;
    double m_g1;
    double m_h0;
    double m_h1;
};

#endif
