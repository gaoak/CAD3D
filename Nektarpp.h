#ifndef NEKTARPP_H
#define NEKTARPP_H
#include"tinyxml2.h"
#include"MeshRegion.h"
#include<string>

class NektarppXml : public MeshRegion {
public:
    NektarppXml(std::string filename, std::string regionname = "xmlRegion_", double tolerance = 1e-6);
    void LoadXml(int nlayers, std::vector<double> targz);
    void UpdateXml();
    void OutXml(std::string name);
    void UpdateXmlComposite();
    void UpdateXmlDomainExpansion();
    
    tinyxml2::XMLDocument m_doc;
protected:
    void LoadDim();
    void LoadModifyPts(int nlayers, std::vector<double> targz);
    void LoadEdge();
    void LoadFace();
    void LoadCell();
    void LoadComposite();
    void LoadModifyCurved(int nlayers, std::vector<double> targz);
    void UpdateXmlPts();
    void UpdateXmlEdge();
    void UpdateXmlFace();
    void UpdateXmlCell();

    std::string m_edgeTag;
    std::string m_faceTag;
    std::string m_bndTye;
};

#endif
