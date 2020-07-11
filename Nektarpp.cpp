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
}

void NektarppXml::LoadXml(int nlayers, std::vector<double> targz) {
    LoadModifyPts(nlayers, targz);
    LoadEdge();
    RebuildEdgesIndex();
    LoadFace();
    RebuildFacesIndex();
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


void NektarppXml::LoadModifyPts(int nlayers, std::vector<double> targz) {
    const char* tagV = "V";
    XMLElement* ptsEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("VERTEX")->FirstChildElement(tagV);
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
        ptsEle = ptsEle->NextSiblingElement(tagV);
        if(ptsEle==nullptr)
        std::cout << "read pts " << m_ptsIndexMax <<  ": " << id << "[" << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
    }
}

void NektarppXml::LoadModifyCurved(int nlayers, std::vector<double> targz) {
    const char* tagE = "E";
    XMLElement* curvEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("CURVED")->FirstChildElement(tagE);
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
        curvEle = curvEle->NextSiblingElement(tagE);
    }
}

void NektarppXml::LoadEdge() {
    const char* tagE = "E";
    XMLElement* edgeEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("EDGE")->FirstChildElement(tagE);
    while(edgeEle != nullptr) {
        int id = edgeEle->IntAttribute("ID");
        const char * estr = edgeEle->GetText();
        std::vector<int> e;
        parserUInt(estr, e);
        m_edges[id] = e;
        if(m_edgeIndexMax<id) m_edgeIndexMax = id;
        edgeEle = edgeEle->NextSiblingElement(tagE);
        if(edgeEle==nullptr)
        std::cout << "read edge " << m_edgeIndexMax <<  ": " << id << "[" << e[0] << ", " << e[1] << std::endl;
    }
}

void NektarppXml::LoadFace() {
    XMLElement* faceEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("FACE")->FirstChildElement();
    while(faceEle != nullptr) {
        int id = faceEle->IntAttribute("ID");
        const char * fstr = faceEle->GetText();
        std::vector<int> f;
        parserUInt(fstr, f);
        m_faces[id] = f;
        if(m_faceIndexMax<id) m_faceIndexMax = id;
        faceEle = faceEle->NextSiblingElement();
        if(faceEle==nullptr)
        std::cout << "read face " << m_faceIndexMax <<  ": " << id << "[" << fstr << std::endl;
    }
}

void NektarppXml::LoadCell() {
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
        std::cout << "read cell " << m_cellIndexMax <<  ": " << id << "[" << cstr << std::endl;
    }
}

void NektarppXml::LoadComposite() {
    const char* tagF = "F";
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
        if(tagF[0] == cstr[1]) {
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
        sprintf(buffer, "%d", m_ptsIndexMax);
        vertex->SetAttribute("ID", buffer);
        ptsRoot->InsertEndChild(vertex);
    }
}

void NektarppXml::UpdateXmlEdge() {
    XMLElement* edgeRoot = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("EDGE");
    edgeRoot->DeleteChildren();
    for(auto it=m_edges.begin(); it!=m_edges.end(); ++it) {
        std::vector<int> e = it->second;
        XMLElement* edge = m_doc.NewElement("E");
        printVector<int>(buffer, "%d ", e);
        XMLText* text = m_doc.NewText(buffer);
        edge->InsertEndChild(text);
        sprintf(buffer, "%d", m_edgeIndexMax);
        edge->SetAttribute("ID", buffer);
        edgeRoot->InsertEndChild(edge);
    }
}

void NektarppXml::UpdateXmlFace() {
    XMLElement* faceRoot = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("FACE");
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
        sprintf(buffer, "%d", m_faceIndexMax);
        face->SetAttribute("ID", buffer);
        faceRoot->InsertEndChild(face);
    }
}

void NektarppXml::UpdateXmlCell() {
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
        sprintf(buffer, "%d", m_cellIndexMax);
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
}
