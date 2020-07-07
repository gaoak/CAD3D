#ifndef MESHREGION_H
#define MESHREGION_H
#include<vector>
#include<set>
#include<map>
#include<string>
#define NFACEMAX 6
//only topology structures
//no physics

class MeshRegion {
public:
    MeshRegion(std::string name, double tolerance = 1.E-6);
    double m_tolerance;
    std::string m_name;
    std::map<int, std::vector<double> > m_pts;
    std::map<int, std::vector<int   > > m_edges;
    std::map<int, std::vector<int   > > m_faces;
    std::map<int, std::vector<int   > > m_cells;
    std::map<int, char> m_cellsType;
    int m_ptsIndexMax;
    int m_edgeIndexMax;
    int m_faceIndexMax;
    int m_cellIndexMax;
    std::map<int, std::set<int> > m_bndComposite;
    std::map<int, std::set<int> > m_domain;
    std::map<int, char > m_domainType;

    std::map<std::set<int>, int> m_edgesIndex;
    std::map<std::set<int>, int> m_facesIndex;
    std::set<int> m_bndPts;

    int pointIsExist(std::vector<double> p, int &pId);
    void rebuildEdgesIndex();
    void rebuildFacesIndex();
    void extractBndPts();
};
#endif // MESHREGION_H
