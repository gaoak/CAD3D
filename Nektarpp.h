#ifndef NEKTARPP_H
#define NEKTARPP_H
#include <tinyxml2.h>
#include "MeshRegion.h"
#include<string>

class NektarppXml : public MeshRegion {
public:
    NektarppXml(std::string name, double tolerance);
    void LoadXml(std::string name, int nlayers, std::vector<double> targz);
    std::map<int, std::set<int> > m_bndComposite;
    int m_bndCompIndexMax;
    void extractBndPts();
    void AddXml(NektarppXml &doc);
    tinyxml2::XMLDocument m_doc;
    void OutXml(std::string name);
};

#endif
