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
    void OutCompAsGeo3D(std::string name, std::vector<std::vector<double> > center, std::vector<double> radius);
    void AddMeshRegion(MeshRegion &doc);
    void GetCellCenter(int i, std::vector<double> & c);
    void GetPts(int index, std::vector<int> & pts, char &type);
    void GetEdgePts(int index, std::vector<int> & pts, char &type);
    void GetFacePts(int index, std::vector<int> & pts, char &type);
    void GetCellPts(int index, std::vector<int> & pts, char &type);
    void GetElementPts(int index, std::vector<int> & pts, char &type);
    void GetBndElementPts(int index, std::vector<int> & pts, char &type);
    void ReorgDomain(std::vector<void*> condition);
    void OutPutSU2(std::string name);
    void GetElements(std::map<int, std::vector<int   > > & elements);
    
    double m_tolerance;
    std::string m_name;
    int m_dim;
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
    
    enum ElementType {
        eHexahedron,
        ePrism,
        eTetrahedron,
        ePyramid,
        eQuadrilateral,
        eTriangle,
        eSegment,
        ePoint
    };
    std::vector<char> ElementTag = {'H', 'R', 'A', 'P', 'Q', 'T', 'S', 'V'};
    std::vector<int>  ElementSU2 = { 12,  13,  10,  14,   9,   5,   3,   0};
    std::map<char, int> ElementTypeMap;

protected:
    void MergePts(MeshRegion &doc, std::map<int, int> &ptsMap);
    void MergeEdge(MeshRegion &doc, std::map<int, int> &ptsMap, std::map<int, int> &edgeMap);
    void MergeFace(MeshRegion &doc, std::map<int, int> &edgeMap, std::map<int, int> &faceMap);
    void AddCell(MeshRegion &doc, std::map<int, int> &faceMap, std::map<int, int> &cellMap);
    void MergeComposite(MeshRegion &doc, std::map<int, int> &faceMap, std::map<int, int> &cellMap);
    int PointHash(std::vector<double> p);
    void InsertBndPts(int index);
    std::vector<double> m_minRange;
    std::vector<double> m_maxRange;
    std::map<std::set<int>, int> m_edgesIndex;
    std::map<std::set<int>, int> m_facesIndex;
    std::map<int, std::set<int> > m_bndPts;
};
#endif // MESHREGION_H
