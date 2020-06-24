#include <tinyxml2.h>
#include "Nektarpp.h"
#include "MeshRegion.h"
#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace tinyxml2;

NektarppXml::NektarppXml(std::string name, double tolerance):MeshRegion(name, tolerance) {
    m_bndCompIndexMax = -1;
}

void NektarppXml::LoadXml(std::string name, int nlayers, std::vector<double> targz) {
    m_doc.LoadFile(name.c_str());
    //load vertex
    const char* tagV = "V";
    XMLElement* ptsEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("VERTEX")->FirstChildElement(tagV);
    while(ptsEle != nullptr) {
        int id = ptsEle->IntAttribute("ID");
        const char * pstr = ptsEle->GetText();
        std::vector<double> p(3, 0.);
        sscanf(pstr, "%lf%lf%lf", &p[0], &p[1], &p[2]);
        p[2] = targz[round(p[2]*nlayers)];
        m_pts[id] = p;
        if(m_ptsIndexMax<id) m_ptsIndexMax = id;
        ptsEle = ptsEle->NextSiblingElement(tagV);
        if(ptsEle==nullptr)
        std::cout << "read pts " << m_ptsIndexMax <<  ": " << id << "[" << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
    }
    //load edge
    const char* tagE = "E";
    XMLElement* edgeEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("EDGE")->FirstChildElement(tagE);
    while(edgeEle != nullptr) {
        int id = edgeEle->IntAttribute("ID");
        const char * estr = edgeEle->GetText();
        std::vector<int> e(2, 0.);
        sscanf(estr, "%d%d", &e[0], &e[1]);
        m_edges[id] = e;
        if(m_edgeIndexMax<id) m_edgeIndexMax = id;
        edgeEle = edgeEle->NextSiblingElement(tagE);
        if(edgeEle==nullptr)
        std::cout << "read edge " << m_edgeIndexMax <<  ": " << id << "[" << e[0] << ", " << e[1] << std::endl;
    }
    //load face
    const char* tagT = "T";
    const char* tagQ = "Q";
    XMLElement* faceEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("FACE")->FirstChildElement();
    while(faceEle != nullptr) {
        int id = faceEle->IntAttribute("ID");
        const char * fstr = faceEle->GetText();
        if(strcmp(tagT, faceEle->Name())==0) {
            std::vector<int> f(3);
            sscanf(fstr, "%d%d%d", &f[0], &f[1], &f[2]);
            m_faces[id] = f;
        } else if(strcmp(tagQ, faceEle->Name())==0) {
            std::vector<int> f(4);
            sscanf(fstr, "%d%d%d%d", &f[0], &f[1], &f[2], &f[3]);
            m_faces[id] = f;
        }
        if(m_faceIndexMax<id) m_faceIndexMax = id;
        faceEle = faceEle->NextSiblingElement();
        if(faceEle==nullptr)
        std::cout << "read face " << m_faceIndexMax <<  ": " << id << "[" << fstr << std::endl;
    }
    //load cells
    const char* tagR = "R";
    const char* tagH = "H";
    XMLElement* cellEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("ELEMENT")->FirstChildElement();
    while(cellEle != nullptr) {
        int id = cellEle->IntAttribute("ID");
        const char * cstr = cellEle->GetText();
        if(strcmp(tagR, cellEle->Name())==0) {
            std::vector<int> f(5);
            sscanf(cstr, "%d%d%d%d%d", &f[0], &f[1], &f[2], &f[3], &f[4]);
            m_cells[id] = f;
        } else if(strcmp(tagH, cellEle->Name())==0) {
            std::vector<int> f(6);
            sscanf(cstr, "%d%d%d%d%d%d", &f[0], &f[1], &f[2], &f[3], &f[4], &f[5]);
            m_cells[id] = f;
        }
        if(m_cellIndexMax<id) m_cellIndexMax = id;
        cellEle = cellEle->NextSiblingElement();
        if(cellEle==nullptr)
        std::cout << "read cell " << m_cellIndexMax <<  ": " << id << "[" << cstr << std::endl;
    }
    //load face composite
    const char* tagF = "F";
    XMLElement* compEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("COMPOSITE")->FirstChildElement();
    while(compEle != nullptr) {
        int id = compEle->IntAttribute("ID");
        const char * cstr = compEle->GetText();
        std::set<int> s;
        m_bndComposite[id] = s;
        if(tagF[0] == cstr[1]) {
            if(m_bndCompIndexMax<id) m_bndCompIndexMax = id;
            std::vector<int> div;
            for(int i=2; i<strlen(cstr); ++i) {
                if(cstr[i]=='[' || cstr[i]==']' || cstr[i]==',') div.push_back(i);
            }
            int k;
            for(int i=0; i<div.size()-1; ++i) {
                std::string cuts(cstr+div[i]+1, div[i+1]-div[i]-1);
                sscanf(cuts.c_str(), "%d", &k);
                m_bndComposite[id].insert(k);
            }
        }
        compEle = compEle->NextSiblingElement();
    }
    rebuildEdgesIndex();
    rebuildFacesIndex();
    extractBndPts();
}

void NektarppXml::extractBndPts() {
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

void NektarppXml::AddXml(NektarppXml &doc) {
    char buffer[10000];
    //merge points
    XMLElement* ptsEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("VERTEX");
    std::map<int, int> ptsMap;
    for(auto it=doc.m_pts.begin(); it!=doc.m_pts.end(); ++it) {
        std::vector<double> p = it->second;
        int pId;
        if(pointIsExist(p, pId)) {
            ptsMap[it->first] = pId;
        } else {
            ptsMap[it->first] = ++m_ptsIndexMax;
            XMLElement* vertex = m_doc.NewElement("V");
            sprintf(buffer, "%14.8lf %14.8lf %14.8lf", p[0], p[1], p[2]);
            XMLText* text = m_doc.NewText(buffer);
            vertex->InsertEndChild(text);
            sprintf(buffer, "%d", m_ptsIndexMax);
            vertex->SetAttribute("ID", buffer);
            ptsEle->InsertEndChild(vertex);
        }
    }
    //merge edges
    XMLElement* edgeEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("EDGE");
    std::map<int, int> edgeMap;
    for(auto it=doc.m_edges.begin(); it!=doc.m_edges.end(); ++it) {
        std::vector<int> e = it->second;
        std::set<int> es;
        for(int i=0; i<e.size(); ++i) {
            e[i] = ptsMap[e[i]];
            es.insert(e[i]);
        }
        if(m_edgesIndex.find(es) != m_edgesIndex.end()) {
            edgeMap[it->first] = m_edgesIndex[es];
        } else {
            edgeMap[it->first] = ++m_edgeIndexMax;
            XMLElement* edge = m_doc.NewElement("E");
            sprintf(buffer, "%d %d", e[0], e[1]);
            XMLText* text = m_doc.NewText(buffer);
            edge->InsertEndChild(text);
            sprintf(buffer, "%d", m_edgeIndexMax);
            edge->SetAttribute("ID", buffer);
            edgeEle->InsertEndChild(edge);
        }
    }
    //merge face
    XMLElement* faceEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("FACE");
    std::map<int, int> faceMap;
    for(auto it=doc.m_faces.begin(); it!=doc.m_faces.end(); ++it) {
        std::vector<int> f = it->second;
        std::set<int> fs;
        for(int i=0; i<f.size(); ++i) {
            f[i] = edgeMap[f[i]];
            fs.insert(f[i]);
        }
        if(m_facesIndex.find(fs)!=m_facesIndex.end()) {
            faceMap[it->first] = m_facesIndex[fs];
            for(auto jt=m_bndComposite.begin(); jt!=m_bndComposite.end(); ++jt) {
                std::set<int> comp = jt->second;
                if(comp.find(m_facesIndex[fs])!=comp.end()) {
                    (jt->second).erase(m_facesIndex[fs]);
                }
            }
            for(auto jt=doc.m_bndComposite.begin(); jt!=doc.m_bndComposite.end(); ++jt) {
                std::set<int> comp = jt->second;
                if(comp.find(it->first)!=comp.end()) {
                    (jt->second).erase(it->first);
                    //std::cout << "remove comp face " << (it->first) << std::endl;
                }
            }
        } else {
            faceMap[it->first] = ++m_faceIndexMax;
            char tag[2];
            tag[1] = 0;
            if(f.size()==3) {
                tag[0] = 'T';
                sprintf(buffer, "%d %d %d", f[0], f[1], f[2]);
            }
            if(f.size()==4) {
                tag[0] = 'Q';
                sprintf(buffer, "%d %d %d %d", f[0], f[1], f[2], f[3]);
            }
            XMLElement* face = m_doc.NewElement(tag);
            XMLText* text = m_doc.NewText(buffer);
            face->InsertEndChild(text);
            sprintf(buffer, "%d", m_faceIndexMax);
            face->SetAttribute("ID", buffer);
            faceEle->InsertEndChild(face);
        }
    }
    //add cells
    XMLElement* cellEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("ELEMENT");
    for(auto it=doc.m_cells.begin(); it!=doc.m_cells.end(); ++it) {
        std::vector<int> c = it->second;
        for(int i=0; i<c.size(); ++i) {
            c[i] = faceMap[c[i]];
        }
        char tag[2];
        tag[1] = 0;
        if(c.size()==5) {
            tag[0] = 'R';
            sprintf(buffer, "%d %d %d %d %d", c[0], c[1], c[2], c[3], c[4]);
        }
        if(c.size()==6) {
            tag[0] = 'H';
            sprintf(buffer, "%d %d %d %d %d %d", c[0], c[1], c[2], c[3], c[4], c[5]);
        }
        XMLElement* cell = m_doc.NewElement(tag);
        XMLText* text = m_doc.NewText(buffer);
        cell->InsertEndChild(text);
        sprintf(buffer, "%d", ++m_cellIndexMax);
        cell->SetAttribute("ID", buffer);
        cellEle->InsertEndChild(cell);
    }
    //modify face composite
    const char* tagF = "F";
    XMLElement* compEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("COMPOSITE")->FirstChildElement();
    while(compEle != nullptr) {
        int id = compEle->IntAttribute("ID");
        const char * cstr = compEle->GetText();
        if(tagF[0] == cstr[1]) {
            std::string text = " F[";
            for(auto jt=m_bndComposite[id].begin(); jt!=m_bndComposite[id].end(); ++jt){
                text += std::to_string(*jt)+",";
            }
            text[text.size()-1] = ']';
            compEle->SetText(text.c_str());
            //std::cout << text << std::endl;
        }
        compEle = compEle->NextSiblingElement();
    }
}

void NektarppXml::OutXml(std::string name) {
    m_doc.SaveFile(name.c_str());
}
