#ifndef NEKTARPP_H
#define NEKTARPP_H
#include"tinyxml2.h"
#include"MeshRegion.h"
#include<string>

class NektarppXml : public MeshRegion {
public:
    NektarppXml(std::string filename, std::string regionname = "xmlRegion_", double tolerance = 1e-6);
    void LoadXml(int nlayers, std::vector<double> targz, double offset0 = 0., double offset1 = 0., bool mapping = false,  int exclude = -1);
    void UpdateXml();
    void OutXml(std::string name);
    void UpdateXmlComposite();
    void UpdateXmlDomainExpansion();
    void LoadWallmapping(std::string filename);
    void DeformPts(void *func);
    tinyxml2::XMLDocument m_doc;
protected:
    void LoadDim();
    void LoadModifyPts(int nlayers, std::vector<double> targz, double offset0 = 0., double offset1 = 0., bool mapping = false,  int exclude = -1);
    void LoadEdge();
    void LoadFace();
    void LoadCell();
    void LoadComposite();
    void LoadModifyCurved(int nlayers, std::vector<double> targz);
    void LoadModifyCurved(void* func);
    void UpdateXmlPts();
    void UpdateXmlEdge();
    void UpdateXmlFace();
    void UpdateXmlCell();
    bool IsWallpoint(std::vector<double> &p, int &id);
    std::string m_edgeTag;
    std::string m_faceTag;
    std::string m_bndTye;
    std::vector<std::vector<double> > m_wallpoints;
    std::map<int, int> m_wallmapping;
};

#endif
