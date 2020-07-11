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
    int PointIsExist(std::vector<double> p, int &pId);
    void RebuildEdgesIndex();
    void RebuildFacesIndex();
    void ExtractBndPts();
    void OutCompAsGeo(std::string name, std::vector<std::vector<double> > center, std::vector<double> radius);
    void AddMeshRegion(MeshRegion &doc);
    void GetCellCenter(int i, std::vector<double> & c);
    void GetFacePts(int index, std::vector<int> & pts);
    void GetCellPts(int index, std::vector<int> & pts, char &type);
    void ReorgDomain(std::vector<void*> condition);
    
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

protected:
    void MergePts(MeshRegion &doc, std::map<int, int> &ptsMap);
    void MergeEdge(MeshRegion &doc, std::map<int, int> &ptsMap, std::map<int, int> &edgeMap);
    void MergeFace(MeshRegion &doc, std::map<int, int> &edgeMap, std::map<int, int> &faceMap);
    void AddCell(MeshRegion &doc, std::map<int, int> &faceMap, std::map<int, int> &cellMap);
    void MergeComposite(MeshRegion &doc, std::map<int, int> &faceMap, std::map<int, int> &cellMap);
    
    std::map<std::set<int>, int> m_edgesIndex;
    std::map<std::set<int>, int> m_facesIndex;
    std::set<int> m_bndPts;
};
#endif // MESHREGION_H
