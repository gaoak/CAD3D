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
    m_doc.LoadFile(name.c_str());
}

double NektarppXml::transformz(double z, int nlayers, std::vector<double> &targz){
    int i = round(z*nlayers);
    if(targz.size()>i) {
        return targz[i];
    } else {
        return z;
    }
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
        m_cellsType[id] = (cellEle->Name())[0];
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

void NektarppXml::LoadXml(std::string name, int nlayers, std::vector<double> targz) {
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

void NektarppXml::mergePts(MeshRegion &doc, std::map<int, int> &ptsMap) {
     XMLElement* ptsEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("VERTEX");
    for(auto it=doc.m_pts.begin(); it!=doc.m_pts.end(); ++it) {
        std::vector<double> p = it->second;
        int pId;
        if(pointIsExist(p, pId)) {
            ptsMap[it->first] = pId;
        } else {
            ptsMap[it->first] = ++m_ptsIndexMax;
            m_pts[m_ptsIndexMax] = p;
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

void NektarppXml::mergeEdge(MeshRegion &doc, std::map<int, int> &ptsMap, std::map<int, int> &edgeMap) {
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
            m_edges[m_edgeIndexMax] = e;
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

void NektarppXml::mergeFace(MeshRegion &doc, std::map<int, int> &edgeMap, std::map<int, int> &faceMap) {
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
            m_faces[m_faceIndexMax] = f;
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

void NektarppXml::addCell(MeshRegion &doc, std::map<int, int> &faceMap, std::map<int, int> &cellMap) {
    XMLElement* cellEle = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("ELEMENT");
    for(auto it=doc.m_cells.begin(); it!=doc.m_cells.end(); ++it) {
        std::vector<int> c = it->second;
        for(int i=0; i<c.size(); ++i) {
            c[i] = faceMap[c[i]];
        }
        ++m_cellIndexMax;
        m_cells[m_cellIndexMax] = c;
        m_cellsType[m_cellIndexMax] = doc.m_cellsType[it->first];
        cellMap[it->first] = m_cellIndexMax;
        char tag[2];
        tag[0] = m_cellsType[m_cellIndexMax];
        tag[1] = 0;
        XMLElement* cell = m_doc.NewElement(tag);
        printVector<int>(buffer, "%d ", c);
        XMLText* text = m_doc.NewText(buffer);
        cell->InsertEndChild(text);
        sprintf(buffer, "%d", m_cellIndexMax);
        cell->SetAttribute("ID", buffer);
        cellEle->InsertEndChild(cell);
    }
}

void NektarppXml::mergeComposite(MeshRegion &doc, std::map<int, int> &faceMap, std::map<int, int> &cellMap, std::set<int> &expansion) {
    XMLElement* compRoot = m_doc.FirstChildElement("NEKTAR")->FirstChildElement("GEOMETRY")->FirstChildElement("COMPOSITE");
    compRoot->DeleteChildren();
    //faces
    int index = 0;
    std::map<int, std::set<int> > bndComposite = m_bndComposite;
    m_bndComposite.clear();
    for(auto it=bndComposite.begin(); it!=bndComposite.end(); ++it) {
        if(bndComposite[it->first].size()>0){
            m_bndComposite[index++] = bndComposite[it->first];
        }
    }
    for(auto it=doc.m_bndComposite.begin(); it!=doc.m_bndComposite.end(); ++it) {
        if(doc.m_bndComposite[it->first].size()>0) {
            std::set<int> s;
            for(auto jt=doc.m_bndComposite[it->first].begin(); jt!=doc.m_bndComposite[it->first].end(); ++jt) {
                s.insert(faceMap[*jt]);
            }
            m_bndComposite[index++] = s;
        }
    }
    for(int i=0; i<index; ++i) {
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
    //cells
    int indexs = index;
    std::map<int, std::set<int> > domain = m_domain;
    std::map<int, char > domaintype = m_domainType;
    m_domain.clear();
    m_domainType.clear();
    for(auto it=domain.begin(); it!=domain.end(); ++it) {
        m_domain[index] = it->second;
        m_domainType[index] = domaintype[it->first];
        ++index;
    }
    for(auto it=doc.m_domain.begin(); it!=doc.m_domain.end(); ++it) {
        std::set<int> s;
        for(auto jt=(it->second).begin(); jt!=(it->second).end(); ++jt) {
            s.insert(cellMap[*jt]);
        }
        m_domain[index] = s;
        m_domainType[index] = doc.m_domainType[it->first];
        ++index;
    }
    for(int i=indexs; i<index; ++i) {
        std::string list = " F" + printComposite(m_domain[i]) + " ";;
        list[1] = m_domainType[i];
        XMLElement* comp = m_doc.NewElement("C");
        XMLText* text = m_doc.NewText(list.c_str());
        comp->InsertEndChild(text);
        comp->SetAttribute("ID", std::to_string(i).c_str());
        expansion.insert(i);
        compRoot->InsertEndChild(comp);
    }
}

void NektarppXml::AddMeshRegion(MeshRegion &doc) {
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
    std::map<int, int> cellMap;
    addCell(doc, faceMap, cellMap);
    //modify face composite
    std::set<int> domain;
    mergeComposite(doc, faceMap, cellMap, domain);
    //modify domain
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
    extractBndPts();
}

void NektarppXml::OutXml(std::string name) {
    m_doc.SaveFile(name.c_str());
}

std::string NektarppXml::printComposite(std::set<int> value) {
    std::string str("[");
    int prev = -100, newele, start = -100;
    for(auto it=value.begin(); it!=value.end(); ++it) {
        int v = *it;
        if(v-prev>1) {
            newele = 1;
        } else {
            newele = 0;
        }
        if(newele) {
            if(start>=0) {
                if(prev-start==1) {
                    str += ",";
                    str += std::to_string(prev);
                } else if(prev-start>1) {
                    str += "-";
                    str += std::to_string(prev);
                }
            }
            if(str.size()>1) str += ",";
            str += std::to_string(v);
            start = v;
        }
        prev = v;
    }
    if(newele==0) {
        if(prev-start==1) {
            str += ",";
        } else {
            str += "-";
        }
        str += std::to_string(prev);
    }
    str += "]";
    return str;
}
