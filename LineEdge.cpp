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

double LineEdge::BuildDiscretes(double ds0, double ds1) {
    double ds;
    if(m_refineType == BOUNDARYLAYER0) {
        ds = ds0;
        m_discretes = std::vector<double>(m_N+1, -1.);
        for(int i=1; i<=m_nBLayers0; ++i) {
            m_discretes[i] = m_discretes[i-1] + ds;
            ds *= m_q0;
            if(ds > (1. - m_discretes[i])/(m_N - i)) {
                m_nBLayers0 = i;
                break;
            }
        }
        ds = (1. - m_discretes[m_nBLayers0])/(m_N - m_nBLayers0);
        for(int i=m_nBLayers0 + 1;i<=m_N;++i) {
            m_discretes[i] = m_discretes[i-1] + ds;
        }
        m_discretes[m_N] = 1.;
    } else if(m_refineType == BOUNDARYLAYER1) {
        ds = ds1;
        m_discretes = std::vector<double>(m_N+1, 1.);
        for(int i=1; i<=m_nBLayers1; ++i) {
            m_discretes[m_N - i] = m_discretes[m_N - i + 1] - ds;
            ds *= m_q1;
            if(ds > (1. + m_discretes[m_N - i])/(m_N - i)) {
                m_nBLayers1 = i;
                break;
            }
        }
        ds = (1. + m_discretes[m_N - m_nBLayers1])/(m_N - m_nBLayers1);
        for(int i=m_nBLayers1 + 1;i<=m_N;++i) {
            m_discretes[m_N - i] = m_discretes[m_N - i + 1] - ds;
        }
        m_discretes[0] = -1.;
    } else if(m_refineType == BOUNDARYLAYER2) {
        m_discretes = std::vector<double>(m_N+1, -1.);
        m_discretes[m_N] = 1.;
        ds = ds0;
        for(int i=1; i<=m_nBLayers0; ++i) {
            m_discretes[i] = m_discretes[i-1] + ds;
            ds *= m_q0;
            if(ds > (1. - m_discretes[i])/(m_N - i)) {
                m_nBLayers0 = i;
                break;
            }
        }
        ds = ds1;
        for(int i=1; i<=m_nBLayers1; ++i) {
            m_discretes[m_N - i] = m_discretes[m_N - i + 1] - ds;
            ds *= m_q1;
            if(ds > (1. + m_discretes[m_N - i])/(m_N - i)) {
                m_nBLayers1 = i;
                break;
            }
        }
        ds = (m_discretes[m_N - m_nBLayers1] - m_discretes[m_nBLayers0])
            /(m_N - m_nBLayers0 - m_nBLayers1);
        if(ds <= 0) {
            std::cout << "error: unsupported BL refine type for edge number [r] " << m_N << std::endl;
        }
        for(int i=m_nBLayers0 + 1;i<m_N - m_nBLayers1;++i) {
            m_discretes[i] = m_discretes[i-1] + ds;
        }
    }
}

double LineEdge::DiscreteStretch(double s, double h0, double h1) {
    double x = (s+1.)/2.*m_N;
    if(x<0.) x = 0.;
    int ilower = floor(x);
    if(ilower==m_N) return 1.;
    BuildDiscretes(h0, h1);
    return (x-ilower)*(m_discretes[ilower+1] - m_discretes[ilower])
           + m_discretes[ilower];
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
    if(m_refineType == BOUNDARYLAYER0
    || m_refineType == BOUNDARYLAYER1
    || m_refineType == BOUNDARYLAYER2) {
        s = DiscreteStretch(s, m_g0*2., m_g1*2.);
    }
    for(int i=0; i<2; ++i) {
        res[i] = 0.5*(1.-s)*m_p0[i] + 0.5*(s+1.)*m_p1[i];
    }
    return res;
}

LineEdge::LineEdge(double* p0, double* p1, int N, int refineType,
                   double h0, double q0, int NBlayers0,
                   double h1, double q1, int NBlayers1) {
    m_p0 = p0;
    m_p1 = p1;
    m_N = N;
    m_refineType = refineType;
    if((BOUNDARYLAYER0 == refineType && (q0 < 1. || N <= NBlayers0))
    || (BOUNDARYLAYER1 == refineType && (q1 < 1. || N <= NBlayers1))
    || (BOUNDARYLAYER2 == refineType && (q0 < 1. || q1 < 1. || N <= NBlayers0 + NBlayers1))) {
        std::cout << "error: unsupported BL refine type for edge number " << N << std::endl;
    }
    m_h0 = h0;
    m_h1 = h1;
    m_q0 = q0;
    m_q1 = q1;
    m_nBLayers0 = NBlayers0;
    m_nBLayers1 = NBlayers1;
}