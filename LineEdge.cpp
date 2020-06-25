#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

#include "LineEdge.h"

static double alphaExpStretch(double h) {
	double alpha = 1.2;
	for(int i=0; i<50; ++i) {
		alpha = 0.5*log(1.+2./h*alpha);
	}
	return alpha;
}

static double alphaSinStretch(double h) {
	double alpha = 1.4;
	for(int i=0; i<50; ++i) {
		alpha = 0.5*acos(h*sin(2.*alpha)/(2.*alpha));
	}
	return alpha;
}

static double alphaSin2Stretch(double h0, double h1, double &alpha, double &beta) {
	alpha = 1.4;
	beta = 0.;
	double tatb = (h0 - h1)/(h0 + h1);
	for(int i=0; i<50; ++i) {
	    alpha = beta - acos(h0*(sin(alpha+beta) - sin(-alpha+beta))/(2.*alpha));
	    beta = atan(tatb/tan(alpha));
	}
	return alpha;
}

double ExpStretch0(double s, double h0) {
	double alpha = alphaExpStretch(h0);
	return -1. + (exp(s*alpha)-exp(-alpha))/(exp(alpha)-exp(-alpha))*2.;
}

double ExpStretch1(double s, double h1) {
	double alpha = alphaExpStretch(h1);
	return -1. + (-exp(-s*alpha)+exp(alpha))/(-exp(-alpha)+exp(alpha))*2.;
}

double QuadStretch0(double s, double h0) {
	double alpha = 0.5*(1. - h0);
	return alpha*(s*s-1.) + s;
}

double QuadStretch1(double s, double h1) {
	double alpha = -0.5*(1. - h1);
	return alpha*(s*s-1.) + s;
}

double SinStretch0(double s, double h0) {
    double alpha = alphaSinStretch(h0);
    double beta = -alpha;
    return -1. + 2.*(sin(alpha*s+beta) - sin(-alpha+beta))/(sin(alpha+beta) - sin(-alpha+beta));
}

double SinStretch1(double s, double h1) {
    double alpha = alphaSinStretch(h1);
    double beta = alpha;
    return -1. + 2.*(sin(alpha*s+beta) - sin(-alpha+beta))/(sin(alpha+beta) - sin(-alpha+beta));
}

double SinStretch2(double s, double h0, double h1) {
    double alpha, beta;
    alphaSin2Stretch(h0, h1, alpha, beta);
    return -1. + 2.*(sin(alpha*s+beta) - sin(-alpha+beta))/(sin(alpha+beta) - sin(-alpha+beta));
}

LineEdge::LineEdge(double* p0, double* p1, int N, int refineType, double h0, double h1) {
    m_p0 = p0;
    m_p1 = p1;
    m_N = N;
    m_refineType = refineType;
    m_h0 = h0*N;
    m_h1 = h1*N;
}

std::vector<double> LineEdge::Evaluate(double s) {
    double Len = sqrt((m_p0[0]-m_p1[0])*(m_p0[0]-m_p1[0]) + (m_p0[1]-m_p1[1])*(m_p0[1]-m_p1[1]));
    m_g0 = m_h0/Len;
    m_g1 = m_h1/Len;
    std::vector<double> res(2, 0.);
    if(m_refineType == EXPREFINE0) {
        s = ExpStretch0(s, m_g0);
    } 
    if(m_refineType == EXPREFINE1) {
        s = ExpStretch1(s, m_g1);
    } 
    if(m_refineType == SINREFINE0) {
        s = SinStretch0(s, m_g0);
    } 
    if(m_refineType == SINREFINE1) {
        s = SinStretch1(s, m_g1);
    } 
    if(m_refineType == SINREFINE2) {
        s = SinStretch2(s, m_g0, m_g1);
    } 
    if(m_refineType == QUDREFINE0) {
        s = QuadStretch0(s, m_g0);
    } 
    if(m_refineType == QUDREFINE1) {
        s = QuadStretch1(s, m_g1);
    } 
    for(int i=0; i<2; ++i) {
        res[i] = 0.5*(1.-s)*m_p0[i] + 0.5*(s+1.)*m_p1[i];
    }
    return res;
}
