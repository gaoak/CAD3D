#ifndef NEKTARPP_H
#define NEKTARPP_H
#include"tinyxml2.h"
#include"MeshRegion.h"
#include<string>

class NektarppXml : public MeshRegion {
public:
    NektarppXml(std::string name, double tolerance);
    void LoadXml(std::string name, int nlayers, std::vector<double> targz);
    std::map<int, std::set<int> > m_bndComposite;
    int m_CompIndexMax;
    void extractBndPts();
    void AddXml(NektarppXml &doc);
    tinyxml2::XMLDocument m_doc;
    void OutXml(std::string name);
private:
    double transformz(double z, int nlayers, std::vector<double> &targz);
    void loadModifyPts(int nlayers, std::vector<double> targz);
    void loadEdge();
    void loadFace();
    void loadCell();
    void loadComposite();
    void loadModifyCurved(int nlayers, std::vector<double> targz);
    void parserInt(const char * cstr, std::vector<int> & value);
    void parserDouble(const char * cstr, std::vector<double> & value);
    template<typename T>
    void printVector(char * str, const char * format, std::vector<T> values) {
        int pos = 0;
        for(int i=0; i<values.size(); ++i) {
            int len = sprintf(str+pos, format, values[i]);
            pos += len;
        }
    }
    void mergePts(NektarppXml &doc, std::map<int, int> &ptsMap);
    void mergeEdge(NektarppXml &doc, std::map<int, int> &ptsMap, std::map<int, int> &edgeMap);
    void mergeFace(NektarppXml &doc, std::map<int, int> &edgeMap, std::map<int, int> &faceMap);
    void addCell(NektarppXml &doc, std::map<int, int> &faceMap, std::vector<int> &Rcell, std::vector<int> &Hcell);
};

#endif