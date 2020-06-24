#include"MeshRegion.h"
#include<cmath>
#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<set>
MeshRegion::MeshRegion(std::string name, double tolerance)
{
    m_name = name;
    m_tolerance = tolerance;
    m_ptsIndexMax = -1;
    m_edgeIndexMax = -1;
    m_faceIndexMax = -1;
    m_cellIndexMax = -1;
}

void MeshRegion::rebuildEdgesIndex(){
    m_edgesIndex.clear();
    for(auto it=m_edges.begin(); it!=m_edges.end(); ++it) {
        std::set<int> p;
        p.insert((it->second)[0]);
        p.insert((it->second)[1]);
        m_edgesIndex[p] = it->first;
    }
}

void MeshRegion::rebuildFacesIndex(){
    m_facesIndex.clear();
    for(auto it=m_faces.begin(); it!=m_faces.end(); ++it) {
        std::set<int> p;
        for(int j=0; j<(it->second).size(); ++j) {
            p.insert((it->second)[j]);
        }
        m_facesIndex[p] = it->first;
    }
}

int MeshRegion::pointIsExist(std::vector<double> p, int &pId) {
    for(auto it = m_bndPts.begin(); it!=m_bndPts.end(); ++it) {
        pId = *it;
        if((fabs(m_pts[pId][0] - p[0]) + fabs(m_pts[pId][1] - p[1]) + fabs(m_pts[pId][2] - p[2]))/3.<m_tolerance) {
            return 1;
        }
    }
    return 0;
}
