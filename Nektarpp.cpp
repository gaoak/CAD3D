#include"tinyxml2.h"
#include"Nektarpp.h"
#include"MeshRegion.h"
#include<string>
#include<cstring>
#include<iostream>
#include<algorithm>
#include<cmath>

using namespace tinyxml2;

static char buffer[10000];

NektarppXml::NektarppXml(std::string name, double tolerance):MeshRegion(name, tolerance) {
    m_CompIndexMax = -1;
}

double NektarppXml::transformz(double z, int nlayers, std::vector<double> &targz){
    int i = round(z*nlayers);
    return targz[i];
}

void NektarppXml::loadModifyPts(int nlayers, std::vector<double> targz) {
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

void NektarppXml::loadModifyCurved(int nlayers, std::vector<double> targz) {
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

void NektarppXml::loadEdge() {
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
    rebuildEdgesIndex();
}

void NektarppXml::loadFace() {
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

void NektarppXml::loadCell() {
    XMLElement* cellEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("ELEMENT")->FirstChildElement();
    while(cellEle != nullptr) {
        int id = cellEle->IntAttribute("ID");
        const char * cstr = cellEle->GetText();
        std::vector<int> f;
        parserUInt(cstr, f);
        m_cells[id] = f;
        if(m_cellIndexMax<id) m_cellIndexMax = id;
        cellEle = cellEle->NextSiblingElement();
        if(cellEle==nullptr)
        std::cout << "read cell " << m_cellIndexMax <<  ": " << id << "[" << cstr << std::endl;
    }
    rebuildFacesIndex();
}

void NektarppXml::parserUInt(const char * cstr, std::vector<int> & value) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if(cstr[i]>='0' && cstr[i]<='9') {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    int k;
    for(int i=0; i<digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);// data  in [s e-1]
        if(sscanf(cuts.c_str(), "%d", &k)<1) {
            std::cout << "error: parser int " << cuts << std::endl;
        }
        if(i>0 && (digs[i] - dige[i-1])==1 && cstr[digs[i]-1]=='-') {
            for(int j=value[value.size()-1]+1; j<k; ++j) {
                value.push_back(j);
            }
        }
        value.push_back(k);
    }
}

void NektarppXml::parserDouble(const char * cstr, std::vector<double> & value) {
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i=0;
    int flag = 0; //digit chunk
    while(1) {
        if((cstr[i]>='0' && cstr[i]<='9') ||
            cstr[i]=='.' ||
            cstr[i]=='e' || cstr[i]=='E' ||
            cstr[i]=='+' || cstr[i]=='-') {
            if(flag==0) {
                digs.push_back(i);
            }
            flag = 1;
        } else {
            if(flag==1) {
                dige.push_back(i);
            }
            flag =  0;
        }
        if(cstr[i]==0) break;
        ++i;
    }
    double k;
    for(int i=0; i<digs.size(); ++i) {
        std::string cuts(cstr+digs[i], dige[i]-digs[i]);
        if(sscanf(cuts.c_str(), "%lf", &k)<1) {
            std::cout << "error: parser double " << cuts << std::endl;
        }
        value.push_back(k);
    }
}

void NektarppXml::loadComposite() {
    const char* tagF = "F";
    XMLElement* compEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("COMPOSITE")->FirstChildElement();
    while(compEle != nullptr) {
        int id = compEle->IntAttribute("ID");
        if(m_CompIndexMax<id) m_CompIndexMax = id;
        const char * cstr = compEle->GetText();
        if(tagF[0] == cstr[1]) {
            std::set<int> s;
            std::vector<int> face;
            parserUInt(cstr, face);
            for(int i=0; i<face.size(); ++i) {
                s.insert(face[i]);
            }
            m_bndComposite[id] = s;
        }
        compEle = compEle->NextSiblingElement();
    }
}

void NektarppXml::LoadXml(std::string name, int nlayers, std::vector<double> targz) {
    m_doc.LoadFile(name.c_str());
    //loading and modify vertex
    loadModifyPts(nlayers, targz);
    //load edge
    loadEdge();
    //load face
    loadFace();
    //load cells
    loadCell();
    //load face composite
    loadComposite();
    //modify curved edges
    loadModifyCurved(nlayers, targz);
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


void NektarppXml::mergePts(NektarppXml &doc, std::map<int, int> &ptsMap) {
     XMLElement* ptsEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("VERTEX");
    for(auto it=doc.m_pts.begin(); it!=doc.m_pts.end(); ++it) {
        std::vector<double> p = it->second;
        int pId;
        if(pointIsExist(p, pId)) {
            ptsMap[it->first] = pId;
        } else {
            ptsMap[it->first] = ++m_ptsIndexMax;
            XMLElement* vertex = m_doc.NewElement("V");
            printVector<double>(buffer, "%20.12e ", p);
            XMLText* text = m_doc.NewText(buffer);
            vertex->InsertEndChild(text);
            sprintf(buffer, "%d", m_ptsIndexMax);
            vertex->SetAttribute("ID", buffer);
            ptsEle->InsertEndChild(vertex);
        }
    }
}

void NektarppXml::mergeEdge(NektarppXml &doc, std::map<int, int> &ptsMap, std::map<int, int> &edgeMap) {
    XMLElement* edgeEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("EDGE");
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
            printVector<int>(buffer, "%d ", e);
            XMLText* text = m_doc.NewText(buffer);
            edge->InsertEndChild(text);
            sprintf(buffer, "%d", m_edgeIndexMax);
            edge->SetAttribute("ID", buffer);
            edgeEle->InsertEndChild(edge);
        }
    }
}

void NektarppXml::mergeFace(NektarppXml &doc, std::map<int, int> &edgeMap, std::map<int, int> &faceMap) {
    XMLElement* faceEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("FACE");
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
                }
            }
        } else {
            faceMap[it->first] = ++m_faceIndexMax;
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
            faceEle->InsertEndChild(face);
        }
    }
}

void NektarppXml::addCell(NektarppXml &doc, std::map<int, int> &faceMap, std::vector<int> &Rcell, std::vector<int> &Hcell) {
    XMLElement* cellEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("ELEMENT");
    for(auto it=doc.m_cells.begin(); it!=doc.m_cells.end(); ++it) {
        std::vector<int> c = it->second;
        for(int i=0; i<c.size(); ++i) {
            c[i] = faceMap[c[i]];
        }
        ++m_cellIndexMax;
        char tag[2];
        tag[1] = 0;
        if(c.size()==5) {
            tag[0] = 'R';
            Rcell.push_back(m_cellIndexMax);
        }
        if(c.size()==6) {
            tag[0] = 'H';
            Hcell.push_back(m_cellIndexMax);
        }
        XMLElement* cell = m_doc.NewElement(tag);
        printVector<int>(buffer, "%d ", c);
        XMLText* text = m_doc.NewText(buffer);
        cell->InsertEndChild(text);
        sprintf(buffer, "%d", m_cellIndexMax);
        cell->SetAttribute("ID", buffer);
        cellEle->InsertEndChild(cell);
    }
}

void NektarppXml::mergeComposite(NektarppXml &doc, std::map<int, int> &faceMap, std::vector<int> &Rcell, std::vector<int> &Hcell, std::vector<int> &domain) {
    std::vector<std::string> facecomposite;
    std::vector<std::string> Hcellcomposite;
    std::vector<std::string> Rcellcomposite;
    //collect basexml composite
    const char* tagF = "F";
    const char* tagH = "H";
    const char* tagR = "R";
    XMLElement* compEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("COMPOSITE")->FirstChildElement();
    while(compEle != nullptr) {
        int id = compEle->IntAttribute("ID");
        const char * cstr = compEle->GetText();
        if(tagF[0] == cstr[1]) {
            std::string list = " F[";
            for(auto jt=m_bndComposite[id].begin(); jt!=m_bndComposite[id].end(); ++jt){
                list += std::to_string(*jt)+",";
            }
            list[list.size()-1] = ']';
            list += ' ';
            facecomposite.push_back(list);
        } else if(tagH[0] == cstr[1]) {
            std::string list(cstr);
            Hcellcomposite.push_back(list);
        } else if(tagR[0] == cstr[1]) {
            std::string list(cstr);
            Rcellcomposite.push_back(list);
        }
        compEle = compEle->NextSiblingElement();
    }
    //add additional face composite
    for(auto it=doc.m_bndComposite.begin(); it!=doc.m_bndComposite.end(); ++it) {
        std::set<int> c = it->second;
        if(c.size()==0) {
            std::cout << "composite " << (it->first) << " has no face left." << std::endl;
            continue;
        }
        std::string list=" F[";
        for(auto jt=c.begin(); jt!=c.end(); ++jt) list += std::to_string(faceMap[*jt]) + ",";
        list[list.size()-1] = ']';
        list += ' ';
        facecomposite.push_back(list);
    }
    //reorganize composite
    XMLElement* compRoot = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("COMPOSITE");
    compRoot->DeleteChildren();
    for(int i=0; i<facecomposite.size(); ++i) {
        XMLElement* comp = m_doc.NewElement("C");
        XMLText* text = m_doc.NewText(facecomposite[i].c_str());
        comp->InsertEndChild(text);
        comp->SetAttribute("ID", std::to_string(i).c_str());
        compRoot->InsertEndChild(comp);
    }
    domain.clear();
    if(Hcellcomposite.size()>0 || Hcell.size()>0) {
        //merge base H composite
        std::string list=" H[";
        for(int i=0; i<Hcellcomposite.size(); ++i) {
            int sb = Hcellcomposite[i].find('[');
            int eb = Hcellcomposite[i].find(']');
            if(eb > sb+1) {
                if(list.size()>3) list += ',';
                list += Hcellcomposite[i].substr(sb+1, eb-sb-1);
            }
        }
        for(int i=0; i<Hcell.size(); ++i) {
            if(list.size()>3) list += ',';
            list += std::to_string(Hcell[i]);
        }
        list += "] ";
        XMLElement* comp = m_doc.NewElement("C");
        XMLText* text = m_doc.NewText(list.c_str());
        comp->InsertEndChild(text);
        comp->SetAttribute("ID", std::to_string(facecomposite.size()+1).c_str());
        domain.push_back(facecomposite.size()+1);
        compRoot->InsertEndChild(comp);
    }
    if(Rcellcomposite.size()>0 || Rcell.size()>0) {
        //merge base R composite
        std::string list=" R[";
        for(int i=0; i<Rcellcomposite.size(); ++i) {
            int sb = Rcellcomposite[i].find('[');
            int eb = Rcellcomposite[i].find(']');
            if(eb > sb+1) {
                if(list.size()>3) list += ',';
                list += Rcellcomposite[i].substr(sb+1, eb-sb-1);
            }
        }
        for(int i=0; i<Rcell.size(); ++i) {
            if(list.size()>3) list += ',';
            list += std::to_string(Rcell[i]);
        }
        list += "] ";
        XMLElement* comp = m_doc.NewElement("C");
        XMLText* text = m_doc.NewText(list.c_str());
        comp->InsertEndChild(text);
        comp->SetAttribute("ID", std::to_string(facecomposite.size()+2).c_str());
        domain.push_back(facecomposite.size()+2);
        compRoot->InsertEndChild(comp);
    }
}

void NektarppXml::AddXml(NektarppXml &doc) {
    std::map<int, int> ptsMap;
    std::map<int, int> edgeMap;
    std::map<int, int> faceMap;
    //merge points
    mergePts(doc, ptsMap);
    //merge edges
    mergeEdge(doc, ptsMap, edgeMap);
    //merge face
    mergeFace(doc, edgeMap, faceMap);
    //add cells
    std::vector<int> Rcell;
    std::vector<int> Hcell;
    addCell(doc, faceMap, Rcell, Hcell);
    //modify face composite
    std::vector<int> domain;
    mergeComposite(doc, faceMap, Rcell, Hcell, domain);
    //modify domain
    {
        XMLElement* domainEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("DOMAIN");
        std::string list=" C[";
        for(int i=0; i<domain.size(); ++i) list += std::to_string(domain[i]) + ",";
        list[list.size()-1] = ']';
        list += ' ';
        domainEle->SetText(list.c_str());
    }
    //modify expansion
    {
        const char * tag= "E";
        XMLElement* expansion = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("EXPANSIONS");
        expansion->DeleteChildren();
        for(int i=0;  i<domain.size(); ++i) {
            XMLElement* exp = m_doc.NewElement(tag);
            sprintf(buffer, "C[%d]", domain[i]);
            exp->SetAttribute("COMPOSITE", buffer);
            exp->SetAttribute("NUMMODES", "2");
            exp->SetAttribute("TYPE", "MODIFIED");
            exp->SetAttribute("FIELDS", "u");
            expansion->InsertEndChild(exp);
        }
    }
}

void NektarppXml::OutXml(std::string name) {
    m_doc.SaveFile(name.c_str());
}
