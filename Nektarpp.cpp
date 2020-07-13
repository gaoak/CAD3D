#include"tinyxml2.h"
#include"Nektarpp.h"
#include"MeshRegion.h"
#include"Util.h"
#include<string>
#include<cstring>
#include<iostream>
#include<algorithm>
#include<cmath>

using namespace tinyxml2;

static char buffer[10000];

NektarppXml::NektarppXml(std::string filename, std::string regionname, double tolerance):MeshRegion(regionname, tolerance) {
    m_doc.LoadFile(filename.c_str());
    m_dim = 3;
    m_edgeTag = "EDGE";
    m_faceTag = "FACE";
    m_bndTye  = "F";
}

void NektarppXml::LoadXml(int nlayers, std::vector<double> targz) {
    LoadDim();
    LoadModifyPts(nlayers, targz);
    LoadEdge();
    LoadFace();
    LoadCell();
    LoadComposite();
    LoadModifyCurved(nlayers, targz);
    ExtractBndPts();
}

void NektarppXml::UpdateXml() {
    UpdateXmlPts();
    UpdateXmlEdge();
    UpdateXmlFace();
    UpdateXmlCell();
    UpdateXmlComposite();
    UpdateXmlDomainExpansion();
}

void NektarppXml::LoadDim() {
    m_dim = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->IntAttribute("DIM");
    if(m_dim == 1) {
        m_bndTye = "V";
        m_edgeTag = "ELEMENT";
    } else if(m_dim == 2) {
        m_bndTye = "E";
        m_faceTag = "ELEMENT";
    } else if(m_dim == 3) {
        m_bndTye = "F";
    }
}

void NektarppXml::LoadModifyPts(int nlayers, std::vector<double> targz) {
    XMLElement* ptsEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("VERTEX")->FirstChildElement();
    while(ptsEle != nullptr) {
        int id = ptsEle->IntAttribute("ID");
        const char * pstr = ptsEle->GetText();
        std::vector<double> p;
        parserDouble(pstr, p);
        p[2] = transformz(p[2], nlayers, targz);
        m_pts[id] = p;
        if(m_ptsIndexMax<id) m_ptsIndexMax = id;
        printVector<double>(buffer, "%20.12e ", p);
        ptsEle->SetText(buffer);
        ptsEle = ptsEle->NextSiblingElement();
        if(ptsEle==nullptr) std::cout << "read pts " << m_ptsIndexMax <<  ": " << id << "[" << p[0] << ", " << p[1] << ", " << p[2] << "]" << std::endl;
    }
}

void NektarppXml::LoadModifyCurved(int nlayers, std::vector<double> targz) {
    XMLElement* curvEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("CURVED")->FirstChildElement();
    while(curvEle != nullptr) {
        int numbers = curvEle->IntAttribute("NUMPOINTS");
        numbers *= 3;
        const char * cstr = curvEle->GetText();
        std::vector<double> values;
        parserDouble(cstr, values);
        for(int i=0; i<numbers; ++i) {
            if((i+1)%3==0) values[i] = transformz(values[i], nlayers, targz);;
        }
        printVector<double>(buffer, "%20.12e ", values);
        curvEle->SetText(buffer);
        curvEle = curvEle->NextSiblingElement();
    }
}

void NektarppXml::LoadEdge() {
    XMLElement* edgeEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement(m_edgeTag.c_str())->FirstChildElement();
    while(edgeEle != nullptr) {
        int id = edgeEle->IntAttribute("ID");
        const char * estr = edgeEle->GetText();
        std::vector<int> e;
        parserUInt(estr, e);
        m_edges[id] = e;
        if(m_edgeIndexMax<id) m_edgeIndexMax = id;
        edgeEle = edgeEle->NextSiblingElement();
        if(edgeEle==nullptr) std::cout << "read edge " << m_edgeIndexMax <<  ": " << id << "[" << e[0] << ", " << e[1] << "]" << std::endl;
    }
    RebuildEdgesIndex();
}

void NektarppXml::LoadFace() {
    if(m_dim<2) return;
    XMLElement* faceEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement(m_faceTag.c_str())->FirstChildElement();
    while(faceEle != nullptr) {
        int id = faceEle->IntAttribute("ID");
        const char * fstr = faceEle->GetText();
        std::vector<int> f;
        parserUInt(fstr, f);
        m_faces[id] = f;
        if(m_faceIndexMax<id) m_faceIndexMax = id;
        faceEle = faceEle->NextSiblingElement();
        if(faceEle==nullptr) std::cout << "read face " << m_faceIndexMax <<  ": " << id << "[" << fstr << "]" << std::endl;
    }
    RebuildFacesIndex();
}

void NektarppXml::LoadCell() {
    if(m_dim<3) return;
    XMLElement* cellEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("ELEMENT")->FirstChildElement();
    while(cellEle != nullptr) {
        int id = cellEle->IntAttribute("ID");
        const char * cstr = cellEle->GetText();
        std::vector<int> f;
        parserUInt(cstr, f);
        m_cells[id] = f;
        m_cellsType[id] = (cellEle->Name())[0];
        if(m_cellIndexMax<id) m_cellIndexMax = id;
        cellEle = cellEle->NextSiblingElement();
        if(cellEle==nullptr)
        std::cout << "read cell " << m_cellIndexMax <<  ": " << id << "[" << cstr << "]" << std::endl;
    }
}

void NektarppXml::LoadComposite() {
    XMLElement* compEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("COMPOSITE")->FirstChildElement();
    while(compEle != nullptr) {
        int id = compEle->IntAttribute("ID");
        const char * cstr = compEle->GetText();
        std::vector<int> face;
        parserUInt(cstr, face);
        std::set<int> s;
        for(int i=0; i<face.size(); ++i) {
            s.insert(face[i]);
        }
        if(m_bndTye[0] == cstr[1]) {
            m_bndComposite[id] = s;
        } else {
            m_domain[id] = s;
            m_domainType[id] = cstr[1];
        }
        compEle = compEle->NextSiblingElement();
    }
}

void NektarppXml::UpdateXmlPts() {
    XMLElement* ptsRoot = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("VERTEX");
    ptsRoot->DeleteChildren();
    for(auto it=m_pts.begin(); it!=m_pts.end(); ++it) {
        std::vector<double> p = it->second;
        XMLElement* vertex = m_doc.NewElement("V");
        printVector<double>(buffer, "%20.12e ", p);
        XMLText* text = m_doc.NewText(buffer);
        vertex->InsertEndChild(text);
        sprintf(buffer, "%d", it->first);
        vertex->SetAttribute("ID", buffer);
        ptsRoot->InsertEndChild(vertex);
    }
}

void NektarppXml::UpdateXmlEdge() {
    XMLElement* edgeRoot = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement(m_edgeTag.c_str());
    edgeRoot->DeleteChildren();
    for(auto it=m_edges.begin(); it!=m_edges.end(); ++it) {
        std::vector<int> e = it->second;
        XMLElement* edge = m_doc.NewElement("E");
        printVector<int>(buffer, "%d ", e);
        XMLText* text = m_doc.NewText(buffer);
        edge->InsertEndChild(text);
        sprintf(buffer, "%d", it->first);
        edge->SetAttribute("ID", buffer);
        edgeRoot->InsertEndChild(edge);
    }
}

void NektarppXml::UpdateXmlFace() {
    if(m_dim<2) return;
    XMLElement* faceRoot = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement(m_faceTag.c_str());
    faceRoot->DeleteChildren();
    for(auto it=m_faces.begin(); it!=m_faces.end(); ++it) {
        std::vector<int> f = it->second;
        char tag[2];
        tag[1] = 0;
        if(f.size()==3) {
            tag[0] = 'T';
        }
        if(f.size()==4) {
            tag[0] = 'Q';
        }
        XMLElement* face = m_doc.NewElement(tag);
        printVector<int>(buffer, "%d ", f);
        XMLText* text = m_doc.NewText(buffer);
        face->InsertEndChild(text);
        sprintf(buffer, "%d", it->first);
        face->SetAttribute("ID", buffer);
        faceRoot->InsertEndChild(face);
    }
}

void NektarppXml::UpdateXmlCell() {
    if(m_dim<3) return;
    XMLElement* cellRoot = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("ELEMENT");
    cellRoot->DeleteChildren();
    for(auto it=m_cells.begin(); it!=m_cells.end(); ++it) {
        std::vector<int> c = it->second;
        char tag[2];
        tag[0] = m_cellsType[it->first];
        tag[1] = 0;
        XMLElement* cell = m_doc.NewElement(tag);
        printVector<int>(buffer, "%d ", c);
        XMLText* text = m_doc.NewText(buffer);
        cell->InsertEndChild(text);
        sprintf(buffer, "%d", it->first);
        cell->SetAttribute("ID", buffer);
        cellRoot->InsertEndChild(cell);
    }
}

void NektarppXml::UpdateXmlComposite() {
    XMLElement* compRoot = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("COMPOSITE");
    compRoot->DeleteChildren();
    for(int i=0; i<m_bndComposite.size(); ++i) {
        if(m_bndComposite[i].size()==0) continue;
        std::string list = " F[";
        list[1] = m_bndTye[0];
        for(auto jt=m_bndComposite[i].begin(); jt!=m_bndComposite[i].end(); ++jt){
            list += std::to_string(*jt)+",";
        }
        list[list.size()-1] = ']';
        list += ' ';
        XMLElement* comp = m_doc.NewElement("C");
        XMLText* text = m_doc.NewText(list.c_str());
        comp->InsertEndChild(text);
        comp->SetAttribute("ID", std::to_string(i).c_str());
        compRoot->InsertEndChild(comp);
    }
    for(auto it=m_domain.begin(); it!=m_domain.end(); ++it) {
        std::string list = " H" + printComposite(it->second) + " ";;
        list[1] = m_domainType[it->first];
        XMLElement* comp = m_doc.NewElement("C");
        XMLText* text = m_doc.NewText(list.c_str());
        comp->InsertEndChild(text);
        comp->SetAttribute("ID", std::to_string(it->first).c_str());
        compRoot->InsertEndChild(comp);
    }
}

void NektarppXml::UpdateXmlDomainExpansion() {
    std::set<int> domain;
    for(auto it=m_domain.begin(); it!=m_domain.end(); ++it) domain.insert(it->first);
    {
        XMLElement* domainEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("DOMAIN");
        std::string list=" C" + printComposite(domain) + " ";
        domainEle->SetText(list.c_str());
    }
    //modify expansion
    {
        const char * tag= "E";
        XMLElement* expansionEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("EXPANSIONS");
        expansionEle->DeleteChildren();
        for(auto it=domain.begin(); it!=domain.end(); ++it) {
            XMLElement* exp = m_doc.NewElement(tag);
            sprintf(buffer, "C[%d]", *it);
            exp->SetAttribute("COMPOSITE", buffer);
            exp->SetAttribute("NUMMODES", "2");
            exp->SetAttribute("TYPE", "MODIFIED");
            exp->SetAttribute("FIELDS", "u");
            expansionEle->InsertEndChild(exp);
        }
    }
}

void NektarppXml::OutXml(std::string name) {
    m_doc.SaveFile(name.c_str());
    std::cout << "write xml file " << name << std::endl;
}
