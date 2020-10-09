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
#define BOUNDARYLAYER0 7
#define BOUNDARYLAYER1 8
#define BOUNDARYLAYER2 9

class LineEdge {
public: 
    LineEdge(double* p0, double* p1, int N, int refineType, double h0, double h1);
    LineEdge(double* p0, double* p1, int N, int refineType,
             double h0, double q0, int NBlayers0,
             double h1, double q1, int NBlayers1);
    std::vector<double> Evaluate(double s);
    double DiscreteStretch(double s, double h0, double h1);
    double BuildDiscretes(double ds, double ds1);
    int m_N;
private:
    double* m_p0;
    double* m_p1;
    int m_refineType;
    double m_g0;
    double m_g1;
    double m_h0;
    double m_h1;
    double m_q0;
    double m_q1;
    double m_nBLayers0;
    double m_nBLayers1;
    std::vector<double> m_discretes;
};

#endif
