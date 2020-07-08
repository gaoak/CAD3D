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

void MeshRegion::extractBndPts() {
    m_bndPts.clear();
    for(auto it=m_bndComposite.begin(); it!=m_bndComposite.end(); ++it) {
        std::set<int> comp = it->second;
        for(auto jt=comp.begin(); jt!=comp.end(); ++jt) {
            std::vector<int> f = m_faces[*jt];
            for(int k=0; k<f.size(); ++k) {
                m_bndPts.insert(m_edges[f[k]][0]);
                m_bndPts.insert(m_edges[f[k]][1]);
            }
        }
    }
}

void MeshRegion::outCompAsGeo(std::string name, std::vector<std::vector<double> > center, std::vector<double> radius) {
    std::set<int> usefaces;
    std::set<int> usepts;
    std::set<int> useedges;
    for(auto it=m_bndComposite.begin(); it!=m_bndComposite.end(); ++it) {
        std::set<int> pts;
        std::cout << "test 0" << std::endl;
        for(auto jt=(it->second).begin(); jt!=(it->second).end(); ++jt) {
            std::vector<int> edge  = m_faces[*jt];
            for(int k=0; k<edge.size(); ++k) {
                pts.insert(m_edges[edge[k]][0]);
                pts.insert(m_edges[edge[k]][1]);
            }
        }
        int select = 0;
        std::cout << "test 1" << std::endl;
        for(auto jt=pts.begin(); jt!=pts.end() && select == 0; ++jt) {
            for(int n=0; n<center.size() && select == 0; ++n) {
                if(fabs(m_pts[*jt][0]-center[n][0]) + fabs(m_pts[*jt][1]-center[n][1]) + fabs(m_pts[*jt][2]-center[n][2]) < 3.*radius[n]) {
                    select = 1;
                    break;
                }
            }
        }
        if(select) {
            for(auto kt=pts.begin(); kt!=pts.end(); ++kt) usepts.insert(*kt);
            for(auto jt=(it->second).begin(); jt!=(it->second).end(); ++jt) {
                std::vector<int> edge  = m_faces[*jt];
                for(int k=0; k<edge.size(); ++k) {
                    useedges.insert(edge[k]);
                }
            }
        }
    }
    //output
    char buffer[1000];
    std::ofstream geo(name.c_str());
    for(auto it=usepts.begin(); it!=usepts.end(); ++it) {
        sprintf(buffer, "Point(%d) = {%20.12e, %20.12e, %20.12e};\n", (*it), m_pts[*it][0], m_pts[*it][1], m_pts[*it][2]);
        geo << buffer ;
    }
    for(auto it=useedges.begin(); it!=useedges.end(); ++it) {
        sprintf(buffer, "Line(%d) = {%d, %d};\n", (*it), m_edges[*it][0], m_edges[*it][1]);
        geo << buffer ;
    }
    for(auto it=usefaces.begin(); it!=usefaces.end(); ++it) {
        std::vector<int> edges = m_faces[*it];
        std::vector<int> sign(edges.size());
        for(int i=0; i<edges.size(); ++i) {
            if(m_edges[edges[i]][1] == m_edges[edges[(i+1)%edges.size()]][0] || m_edges[edges[i]][1] == m_edges[edges[(i+1)%edges.size()]][1]) {
                sign[i] = 1;
            } else {
                sign[i] =-1;
            }
        }
        sprintf(buffer, "Line Loop(%d) = {", (*it));
        geo << buffer ;
        for(int i=0; i<edges.size(); ++i) {
            if(sign[i]<0) geo << "-";
            geo << edges[i];
            if(i<edges.size()-1) geo << ",";
        }
        geo << "};\n";
        sprintf(buffer, "Plane Surface(%d) = {%d};\n", (*it), (*it));
        geo << buffer;
    }
    auto it=usefaces.begin();
    geo << "Surface Loop(0) = {" << (*it);
    for(++it; it!=usefaces.end(); ++it) {
        geo << "," << (*it);
    }
    geo << "};\n";
    geo << "Volume(0)={0};\n";
}
