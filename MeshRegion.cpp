#include "MeshRegion.h"
#include "Util.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
MeshRegion::MeshRegion(std::string name, double tolerance) {
  m_name = name;
  m_tolerance = tolerance;
  m_ptsIndexMax = -1;
  m_edgeIndexMax = -1;
  m_faceIndexMax = -1;
  m_cellIndexMax = -1;
  for (int i = 0; i < ElementTag.size(); ++i) {
    ElementTypeMap[ElementTag[i]] = i;
  }
  m_minRange = std::vector<double>(3, 1E100);
  m_maxRange = std::vector<double>(3, -1E100);
}

void MeshRegion::RebuildEdgesIndex() {
  m_edgesIndex.clear();
  for (auto it = m_edges.begin(); it != m_edges.end(); ++it) {
    std::set<int> p;
    p.insert((it->second)[0]);
    p.insert((it->second)[1]);
    m_edgesIndex[p] = it->first;
  }
}

void MeshRegion::RebuildFacesIndex() {
  m_facesIndex.clear();
  for (auto it = m_faces.begin(); it != m_faces.end(); ++it) {
    std::set<int> p;
    for (int j = 0; j < (it->second).size(); ++j) {
      p.insert((it->second)[j]);
    }
    m_facesIndex[p] = it->first;
  }
}

int MeshRegion::PointIsExist(std::vector<double> p, int &pId) {
  int index = PointHash(p);
  for (int i = index - 1; i <= index + 1; ++i) {
    if (m_bndPts.find(i) != m_bndPts.end()) {
      for (auto it = m_bndPts[i].begin(); it != m_bndPts[i].end(); ++it) {
        pId = *it;
        if ((fabs(m_pts[pId][0] - p[0]) + fabs(m_pts[pId][1] - p[1]) +
             fabs(m_pts[pId][2] - p[2])) /
                3. <
            m_tolerance) {
          return 1;
        }
      }
    }
  }
  return 0;
}

void MeshRegion::InsertBndPts(int index) {
  if (m_pts.find(index) == m_pts.end())
    return;
  std::vector<double> p = m_pts[index];
  int hash = PointHash(p);
  if (m_bndPts.find(hash) == m_bndPts.end()) {
    std::set<int> ts;
    m_bndPts[hash] = ts;
  }
  m_bndPts[hash].insert(index);
}

void MeshRegion::ExtractBndPts() {
  m_bndPts.clear();
  switch (m_dim) {
  case 3:
    for (auto it = m_bndComposite.begin(); it != m_bndComposite.end(); ++it) {
      std::set<int> comp = it->second;
      for (auto jt = comp.begin(); jt != comp.end(); ++jt) {
        std::vector<int> f = m_faces[*jt];
        for (int k = 0; k < f.size(); ++k) {
          InsertBndPts(m_edges[f[k]][0]);
          InsertBndPts(m_edges[f[k]][1]);
        }
      }
    }
    break;
  case 2:
    for (auto it = m_bndComposite.begin(); it != m_bndComposite.end(); ++it) {
      std::set<int> comp = it->second;
      for (auto jt = comp.begin(); jt != comp.end(); ++jt) {
        InsertBndPts(m_edges[*jt][0]);
        InsertBndPts(m_edges[*jt][1]);
      }
    }
    break;
  case 1:
    for (auto it = m_bndComposite.begin(); it != m_bndComposite.end(); ++it) {
      std::set<int> comp = it->second;
      for (auto jt = comp.begin(); jt != comp.end(); ++jt) {
        InsertBndPts(*jt);
      }
    }
    break;
  default:
    std::cout << "error: unsupported dimension " << m_dim << std::endl;
  }
}

void MeshRegion::OutCompAsGeo3D(std::string name,
                                std::vector<std::vector<double>> center,
                                std::vector<double> radius) {
  std::set<int> usefaces;
  std::set<int> usepts;
  std::set<int> useedges;
  for (auto it = m_bndComposite.begin(); it != m_bndComposite.end(); ++it) {
    std::set<int> pts;
    for (auto jt = (it->second).begin(); jt != (it->second).end(); ++jt) {
      std::vector<int> edge = m_faces[*jt];
      for (int k = 0; k < edge.size(); ++k) {
        pts.insert(m_edges[edge[k]][0]);
        pts.insert(m_edges[edge[k]][1]);
      }
    }
    int select = 0;
    for (auto jt = pts.begin(); jt != pts.end() && select == 0; ++jt) {
      for (int n = 0; n < center.size() && select == 0; ++n) {
        if (fabs(m_pts[*jt][0] - center[n][0]) +
                fabs(m_pts[*jt][1] - center[n][1]) +
                fabs(m_pts[*jt][2] - center[n][2]) <
            3. * radius[n]) {
          select = 1;
          break;
        }
      }
    }
    if (select) {
      for (auto kt = pts.begin(); kt != pts.end(); ++kt)
        usepts.insert(*kt);
      for (auto jt = (it->second).begin(); jt != (it->second).end(); ++jt) {
        usefaces.insert(*jt);
        std::vector<int> edge = m_faces[*jt];
        for (int k = 0; k < edge.size(); ++k) {
          useedges.insert(edge[k]);
        }
      }
    }
  }
  // output
  char buffer[1000];
  std::ofstream geo(name.c_str());
  for (auto it = usepts.begin(); it != usepts.end(); ++it) {
    sprintf(buffer, "Point(%d) = {%20.12e, %20.12e, %20.12e};\n", (*it),
            m_pts[*it][0], m_pts[*it][1], m_pts[*it][2]);
    geo << buffer;
  }
  for (auto it = useedges.begin(); it != useedges.end(); ++it) {
    sprintf(buffer, "Line(%d) = {%d, %d};\n", (*it), m_edges[*it][0],
            m_edges[*it][1]);
    geo << buffer;
  }
  for (auto it = usefaces.begin(); it != usefaces.end(); ++it) {
    std::vector<int> edges = m_faces[*it];
    std::vector<int> sign(edges.size());
    for (int i = 0; i < edges.size(); ++i) {
      if (m_edges[edges[i]][1] == m_edges[edges[(i + 1) % edges.size()]][0] ||
          m_edges[edges[i]][1] == m_edges[edges[(i + 1) % edges.size()]][1]) {
        sign[i] = 1;
      } else {
        sign[i] = -1;
      }
    }
    sprintf(buffer, "Line Loop(%d) = {", (*it));
    geo << buffer;
    for (int i = 0; i < edges.size(); ++i) {
      if (sign[i] < 0)
        geo << "-";
      geo << edges[i];
      if (i < edges.size() - 1)
        geo << ",";
    }
    geo << "};\n";
    sprintf(buffer, "Plane Surface(%d) = {%d};\n", (*it), (*it));
    geo << buffer;
  }
  auto it = usefaces.begin();
  geo << "Surface Loop(0) = {" << (*it);
  for (++it; it != usefaces.end(); ++it) {
    geo << "," << (*it);
  }
  geo << "};\n";
  geo << "Volume(0)={0};\n";
}

void MeshRegion::MergePts(MeshRegion &doc, std::map<int, int> &ptsMap) {
  ptsMap.clear();
  for (auto it = doc.m_pts.begin(); it != doc.m_pts.end(); ++it) {
    std::vector<double> p = it->second;
    int pId;
    if (PointIsExist(p, pId)) {
      ptsMap[it->first] = pId;
      if (m_dim == 1) {
        for (auto jt = m_bndComposite.begin(); jt != m_bndComposite.end();
             ++jt) {
          std::set<int> comp = jt->second;
          if (comp.find(pId) != comp.end()) {
            (jt->second).erase(pId);
          }
        }
        for (auto jt = doc.m_bndComposite.begin();
             jt != doc.m_bndComposite.end(); ++jt) {
          std::set<int> comp = jt->second;
          if (comp.find(it->first) != comp.end()) {
            (jt->second).erase(it->first);
          }
        }
      }
    } else {
      ptsMap[it->first] = ++m_ptsIndexMax;
      m_pts[m_ptsIndexMax] = p;
    }
  }
}

void MeshRegion::MergeEdge(MeshRegion &doc, std::map<int, int> &ptsMap,
                           std::map<int, int> &edgeMap) {
  edgeMap.clear();
  for (auto it = doc.m_edges.begin(); it != doc.m_edges.end(); ++it) {
    std::vector<int> e = it->second;
    std::set<int> es;
    for (int i = 0; i < e.size(); ++i) {
      e[i] = ptsMap[e[i]];
      es.insert(e[i]);
    }
    if (m_edgesIndex.find(es) != m_edgesIndex.end()) {
      edgeMap[it->first] = m_edgesIndex[es];
      if (m_dim == 2) {
        for (auto jt = m_bndComposite.begin(); jt != m_bndComposite.end();
             ++jt) {
          std::set<int> comp = jt->second;
          if (comp.find(m_edgesIndex[es]) != comp.end()) {
            (jt->second).erase(m_edgesIndex[es]);
          }
        }
        for (auto jt = doc.m_bndComposite.begin();
             jt != doc.m_bndComposite.end(); ++jt) {
          std::set<int> comp = jt->second;
          if (comp.find(it->first) != comp.end()) {
            (jt->second).erase(it->first);
          }
        }
      }
    } else {
      edgeMap[it->first] = ++m_edgeIndexMax;
      m_edges[m_edgeIndexMax] = e;
    }
  }
}

void MeshRegion::MergeFace(
    MeshRegion &doc, std::map<int, int> &edgeMap,
    std::map<int, int>
        &faceMap) { // erase a face from boundary if it is not real boundary
  faceMap.clear();
  for (auto it = doc.m_faces.begin(); it != doc.m_faces.end(); ++it) {
    std::vector<int> f = it->second;
    std::set<int> fs;
    for (int i = 0; i < f.size(); ++i) {
      f[i] = edgeMap[f[i]];
      fs.insert(f[i]);
    }
    if (m_facesIndex.find(fs) != m_facesIndex.end()) {
      faceMap[it->first] = m_facesIndex[fs];
      if (m_dim == 3) {
        for (auto jt = m_bndComposite.begin(); jt != m_bndComposite.end();
             ++jt) {
          std::set<int> comp = jt->second;
          if (comp.find(m_facesIndex[fs]) != comp.end()) {
            (jt->second).erase(m_facesIndex[fs]);
          }
        }
        for (auto jt = doc.m_bndComposite.begin();
             jt != doc.m_bndComposite.end(); ++jt) {
          std::set<int> comp = jt->second;
          if (comp.find(it->first) != comp.end()) {
            (jt->second).erase(it->first);
          }
        }
      }
    } else {
      faceMap[it->first] = ++m_faceIndexMax;
      m_faces[m_faceIndexMax] = f;
    }
  }
}

void MeshRegion::AddCell(MeshRegion &doc, std::map<int, int> &faceMap,
                         std::map<int, int> &cellMap) {
  cellMap.clear();
  for (auto it = doc.m_cells.begin(); it != doc.m_cells.end(); ++it) {
    std::vector<int> c = it->second;
    for (int i = 0; i < c.size(); ++i) {
      c[i] = faceMap[c[i]];
    }
    ++m_cellIndexMax;
    m_cells[m_cellIndexMax] = c;
    m_cellsType[m_cellIndexMax] = doc.m_cellsType[it->first];
    cellMap[it->first] = m_cellIndexMax;
  }
}

void MeshRegion::MergeComposite(MeshRegion &doc, std::map<int, int> &faceMap,
                                std::map<int, int> &cellMap) {
  // faces
  int index = 0;
  std::map<int, std::set<int>> bndComposite = m_bndComposite;
  m_bndComposite.clear();
  for (auto it = bndComposite.begin(); it != bndComposite.end(); ++it) {
    if (bndComposite[it->first].size() > 0) {
      m_bndComposite[index++] = bndComposite[it->first];
    }
  }
  for (auto it = doc.m_bndComposite.begin(); it != doc.m_bndComposite.end();
       ++it) {
    if (doc.m_bndComposite[it->first].size() > 0) {
      std::set<int> s;
      for (auto jt = doc.m_bndComposite[it->first].begin();
           jt != doc.m_bndComposite[it->first].end(); ++jt) {
        s.insert(faceMap[*jt]);
      }
      m_bndComposite[index++] = s;
    }
  }
  // cells
  int indexs = index;
  std::map<int, std::set<int>> domain = m_domain;
  std::map<int, char> domaintype = m_domainType;
  m_domain.clear();
  m_domainType.clear();
  for (auto it = domain.begin(); it != domain.end(); ++it) {
    m_domain[index] = it->second;
    m_domainType[index] = domaintype[it->first];
    ++index;
  }
  for (auto it = doc.m_domain.begin(); it != doc.m_domain.end(); ++it) {
    std::set<int> s;
    for (auto jt = (it->second).begin(); jt != (it->second).end(); ++jt) {
      s.insert(cellMap[*jt]);
    }
    m_domain[index] = s;
    m_domainType[index] = doc.m_domainType[it->first];
    ++index;
  }
}

void MeshRegion::AddMeshRegion(MeshRegion &doc) {
  m_name = m_name + doc.m_name;
  std::map<int, int> ptsMap;
  MergePts(doc, ptsMap);
  std::map<int, int> edgeMap;
  MergeEdge(doc, ptsMap, edgeMap);
  RebuildEdgesIndex();
  std::map<int, int> faceMap;
  MergeFace(doc, edgeMap, faceMap);
  RebuildFacesIndex();
  std::map<int, int> cellMap;
  AddCell(doc, faceMap, cellMap);
  MergeComposite(doc, faceMap, cellMap);
  ExtractBndPts();
}

void MeshRegion::GetElementPts(int index, std::vector<int> &pts, char &type) {
  switch (m_dim) {
  case 1:
    GetEdgePts(index, pts, type);
    break;
  case 2:
    GetFacePts(index, pts, type);
    break;
  case 3:
    GetCellPts(index, pts, type);
    break;
  default:
    std::cout << "error: dimension " << m_dim << " not supported." << std::endl;
  }
}

void MeshRegion::GetBndElementPts(int index, std::vector<int> &pts,
                                  char &type) {
  switch (m_dim) {
  case 1:
    GetPts(index, pts, type);
    break;
  case 2:
    GetEdgePts(index, pts, type);
    break;
  case 3:
    GetFacePts(index, pts, type);
    break;
  default:
    std::cout << "error: dimension " << m_dim << " not supported." << std::endl;
  }
}

void MeshRegion::GetElements(std::map<int, std::vector<int>> &elements) {
  switch (m_dim) {
  case 1:
    elements = m_edges;
    break;
  case 2:
    elements = m_faces;
    break;
  case 3:
    elements = m_cells;
    break;
  }
}

void MeshRegion::GetPts(int index, std::vector<int> &pts, char &type) {
  pts.clear();
  type = 0;
  if (m_pts.find(index) == m_pts.end()) {
    return;
  }
  pts.push_back(index);
  type = ElementTag[ePoint];
}

void MeshRegion::GetEdgePts(int index, std::vector<int> &pts, char &type) {
  pts.clear();
  type = 0;
  if (m_edges.find(index) == m_edges.end()) {
    return;
  }
  pts.push_back(m_edges[index][0]);
  pts.push_back(m_edges[index][1]);
  type = ElementTag[eSegment];
}

void MeshRegion::GetFacePts(int index, std::vector<int> &pts, char &type) {
  pts.clear();
  type = 0;
  if (m_faces.find(index) == m_faces.end()) {
    return;
  }
  int ne = m_faces[index].size();
  if (ne == 3)
    type = ElementTag[eTriangle];
  if (ne == 4)
    type = ElementTag[eQuadrilateral];
  for (int i = 0; i < ne; ++i) {
    int e0 = m_faces[index][i];
    int e1 = m_faces[index][(i + 1) % ne];
    if (m_edges[e0][0] == m_edges[e1][0] || m_edges[e0][0] == m_edges[e1][1]) {
      pts.push_back(m_edges[e0][1]);
    } else {
      pts.push_back(m_edges[e0][0]);
    }
  }
}

void MeshRegion::GetFaceUnitNorm(int index, std::vector<double> &norm,
                                 char &type) {
  norm = std::vector<double>(3, 0.);
  if (m_dim <= 2) {
    norm[2] = 1.;
    return;
  }
  std::vector<int> pts;
  GetFacePts(index, pts, type);
  std::vector<double> t1(3);
  std::vector<double> t2(3);
  for (int i = 0; i < 3; ++i) {
    t1[i] = m_pts[pts[1]][i] - m_pts[pts[0]][i];
    t2[i] = m_pts[pts[2]][i] - m_pts[pts[1]][i];
  }
  norm[0] = t1[1] * t2[2] - t1[2] * t2[1];
  norm[1] = t1[2] * t2[0] - t1[0] * t2[2];
  norm[2] = t1[0] * t2[1] - t1[1] * t2[0];
  double l = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
  for (int i = 0; i < 3; ++i) {
    norm[i] /= l;
  }
}

void MeshRegion::GetFaceEdges(int index, std::vector<int> &edges, char &type) {
  edges.clear();
  type = 0;
  if (m_faces.find(index) == m_faces.end()) {
    return;
  }
  int ne = m_faces[index].size();
  if (ne == 3)
    type = ElementTag[eTriangle];
  if (ne == 4)
    type = ElementTag[eQuadrilateral];
  for (int i = 0; i < ne; ++i) {
    int e0 = m_faces[index][i];
    edges.push_back(e0);
  }
}

void MeshRegion::GetCellEdges(int index, std::vector<int> &edges, char &type) {
  edges.clear();
  type = 0;
  if (m_cells.find(index) == m_cells.end()) {
    return;
  }
  type = m_cellsType[index];
  std::vector<int> cell = m_cells[index];
  int ne = cell.size();
  std::set<int> edgeset;
  std::vector<int> faceedges;
  char bottype;
  for (int i = 0; i < ne; ++i) {
    GetFaceEdges(cell[i], faceedges, bottype);
    for (int j = 0; j < faceedges.size(); ++j) {
      edgeset.insert(faceedges[j]);
    }
  }
  for (auto it = edgeset.begin(); it != edgeset.end(); ++it) {
    edges.push_back(*it);
  }
}

double MeshRegion::CheckDirection(std::vector<int> pts, int p1) {
  std::vector<double> t1(3);
  std::vector<double> t2(3);
  for (int i = 0; i < 3; ++i) {
    t1[i] = m_pts[pts[1]][i] - m_pts[pts[0]][i];
    t2[i] = m_pts[pts[2]][i] - m_pts[pts[1]][i];
  }
  std::vector<double> norm(3);
  norm[0] = t1[1] * t2[2] - t1[2] * t2[1];
  norm[1] = t1[2] * t2[0] - t1[0] * t2[2];
  norm[2] = t1[0] * t2[1] - t1[1] * t2[0];
  double res = 0.;
  for (int i = 0; i < 3; ++i) {
    res += norm[i] * (m_pts[p1][i] - m_pts[pts[1]][i]);
  }
  return res;
}

void MeshRegion::GetCellPts(int index, std::vector<int> &pts, char &type) {
  pts.clear();
  type = 0;
  if (m_cells.find(index) == m_cells.end()) {
    return;
  }
  type = m_cellsType[index];
  std::vector<int> cell = m_cells[index];
  int ne = cell.size();
  char bottype;
  if (type == ElementTag[ePrism] || type == ElementTag[eHexahedron]) {
    std::vector<int> bot;
    for (int i = 0; i < ne; ++i) {
      if (m_faces[cell[i]].size() == 3)
        bot.push_back(cell[i]);
    }
    if (bot.size() == 0) {
      bot.push_back(cell[0]);
      bot.push_back(cell[5]);
    }
    GetFacePts(bot[0], pts, bottype);
    std::vector<int> p1;
    GetFacePts(bot[1], p1, bottype);
    if (CheckDirection(pts, p1[0]) < 0.) {
      std::reverse(pts.begin(), pts.end());
    }
    for (int i = 0; i < p1.size(); ++i) {
      for (int j = 0; j < p1.size(); ++j) {
        std::set<int> es;
        es.insert(pts[i]);
        es.insert(p1[j]);
        if (m_edgesIndex.find(es) != m_edgesIndex.end()) {
          pts.push_back(p1[j]);
          break;
        }
      }
    }
  } else if (type == ElementTag[eTetrahedron] || type == ElementTag[ePyramid]) {
    std::vector<int> bot(2);
    bot[0] = cell[0];
    bot[1] = cell[1];
    for (int i = 1; i < ne; ++i) {
      if (m_faces[cell[i]].size() == 4) {
        bot[0] = cell[i];
        bot[1] = cell[(i + 1) % ne];
      }
    }

    GetFacePts(bot[0], pts, bottype);
    std::vector<int> p1;
    GetFacePts(bot[1], p1, bottype);
    if (CheckDirection(pts, p1[0]) < 0.) {
      std::reverse(pts.begin(), pts.end());
    }
    for (int i = 0; i < p1.size(); ++i) {
      int j = 0;
      for (; j < pts.size(); ++j) {
        if (pts[j] == p1[i])
          break;
      }
      if (j == pts.size()) {
        pts.push_back(p1[i]);
        break;
      }
    }
  } else {
    std::cout << "cell not support. " << std::endl;
  }
  std::set<int> temp;
  for (int i = 0; i < m_cells[index].size(); ++i) {
    std::vector<int> face = m_faces[m_cells[index][i]];
  }
}

void MeshRegion::GetCellCenter(int index, std::vector<double> &c) {
  c = std::vector<double>(3, DBL_MAX);
  std::vector<int> pts;
  char type;
  GetCellPts(index, pts, type);
  if (pts.size() < 4) {
    std::cout << "warning cell " << index << " not found." << std::endl;
    return;
  }
  for (int i = 0; i < 3; ++i) {
    c[i] = 0.;
    for (int j = 0; j < pts.size(); ++j) {
      c[i] += m_pts[pts[j]][i];
    }
    c[i] /= pts.size();
  }
}

int MeshRegion::FindSingularElements(double angle) {
  m_singularElements.clear();
  m_singularFace.clear();
  m_singularEdge.clear();
  m_singularPts.clear();
  if (m_dim < 3)
    return 0;
  double cang = cos(angle);
  std::map<int, int> singularEdge;
  char type;
  std::vector<int> edges;
  for (auto comp = m_bndComposite.begin(); comp != m_bndComposite.end();
       ++comp) {
    for (auto face = comp->second.begin(); face != comp->second.end(); ++face) {
      GetFaceEdges(*face, edges, type);
      for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
        if (singularEdge.find(*edge) == singularEdge.end()) {
          singularEdge[*edge] = *face;
        } else {
          std::vector<double> n1(3), n2(3);
          GetFaceUnitNorm(*face, n1, type);
          GetFaceUnitNorm(singularEdge[*edge], n2, type);
          double a = fabs(n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]);
          if (a > cang) {
            singularEdge.erase(*edge);
          } else {
            m_singularFace.insert(*face);
            m_singularFace.insert(singularEdge[*edge]);
          }
        }
      }
    }
  }
  std::vector<int> pts;
  for (auto edge = singularEdge.begin(); edge != singularEdge.end(); ++edge) {
    m_singularEdge.insert(edge->first);
    GetEdgePts(edge->first, pts, type);
    for (auto it = pts.begin(); it != pts.end(); ++it) {
      m_singularPts.insert(*it);
    }
  }
  std::ofstream singularfile("singular.dat");
  singularfile << "variables = x, y, z" << std::endl;
  for (auto it = m_singularPts.begin(); it != m_singularPts.end(); ++it) {
    singularfile << m_pts[*it][0] << " " << m_pts[*it][1] << " "
                 << m_pts[*it][2] << std::endl;
  }
  for (auto ele = m_cells.begin(); ele != m_cells.end(); ++ele) {
    GetCellPts(ele->first, pts, type);
    for (auto it = pts.begin(); it != pts.end(); ++it) {
      if (m_singularPts.find(*it) != m_singularPts.end()) {
        m_singularElements.insert(ele->first);
      }
    }
  }
  return 0;
}

void MeshRegion::ReorgDomain(std::vector<void *> condition,
                             bool detectSingular) {
  if (detectSingular) {
    FindSingularElements();
  }
  int index = 0;
  for (auto it = m_bndComposite.begin(); it != m_bndComposite.end(); ++it) {
    if (index < it->first)
      index = it->first;
  }
  if (index < 100)
    index = 100;
  else if (index < 1000)
    index = 1000;
  else if (index < 1000000)
    index = 1000000;

  void *alwaysTrue;
  condition.push_back(alwaysTrue);
  std::set<int> cells;
  switch (m_dim) {
  case 3:
    for (auto it = m_cells.begin(); it != m_cells.end(); ++it) {
      cells.insert(it->first);
    }
    break;
  case 2:
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it) {
      cells.insert(it->first);
    }
    break;
  case 1:
    for (auto it = m_edges.begin(); it != m_edges.end(); ++it) {
      cells.insert(it->first);
    }
    break;
  default:
    std::cout << "error: dimension " << m_dim << " not supported." << std::endl;
  }
  m_domain.clear();
  m_domainType.clear();
  for (int i = 0; i < condition.size(); ++i) {
    std::vector<std::set<int>> tmpset;
    for (int j = 0; j < ElementTag.size(); ++j)
      tmpset.push_back(std::set<int>());
    std::vector<int> toclear;
    double (*con)(double, double, double) =
        (double (*)(double, double, double))condition[i];
    for (auto it = cells.begin(); it != cells.end(); ++it) {
      std::vector<int> pts;
      char type;
      GetElementPts(*it, pts, type);
      for (int j = 0; j < pts.size(); ++j) {
        std::vector<double> c = m_pts[pts[j]];
        if (i == condition.size() - 1 || con(c[0], c[1], c[2]) > 0.) {
          if (detectSingular && i == 0 &&
              m_singularElements.find(*it) == m_singularElements.end()) {
            break;
          }
          tmpset[ElementTypeMap[type]].insert(*it);
          toclear.push_back(*it);
          break;
        }
      }
    }
    for (int j = 0; j < ElementTag.size(); ++j) {
      if (tmpset[j].size() > 0) {
        m_domain[index + j] = tmpset[j];
        m_domainType[index + j] = ElementTag[j];
      }
    }
    if (toclear.size() == 0)
      continue;
    index += 10;
    for (int j = 0; j < toclear.size(); ++j) {
      cells.erase(toclear[j]);
    }
  }
}

int MeshRegion::FindSharedSingluarEdges(std::set<int> &bnd1,
                                        std::set<int> &bnd2, int &shared,
                                        int &singular) {
  shared = 0;
  singular = 0;
  double nsface = 0;
  if (bnd1.empty() || bnd2.empty())
    return nsface;
  std::set<int> edge1, edge2;
  char type;
  std::vector<int> edges;
  for (auto it = bnd1.begin(); it != bnd1.end(); ++it) {
    GetFaceEdges(*it, edges, type);
    for (auto e = edges.begin(); e != edges.end(); ++e) {
      edge1.insert(*e);
    }
  }
  for (auto it = bnd2.begin(); it != bnd2.end(); ++it) {
    GetFaceEdges(*it, edges, type);
    for (auto e = edges.begin(); e != edges.end(); ++e) {
      edge2.insert(*e);
    }
  }
  for (auto it = edge2.begin(); it != edge2.end(); ++it) {
    if (edge1.find(*it) != edge1.end()) {
      ++shared;
      if (m_singularEdge.find(*it) != m_singularEdge.end()) {
        ++singular;
      }
    }
  }
  for (auto it = bnd2.begin(); it != bnd2.end(); ++it) {
    if (bnd1.find(*it) != bnd1.end()) {
      ++nsface;
    }
  }
  return nsface;
}
// (< angle) (cos > cos angle) no sigular
bool MeshRegion::ReorgBoundary(double angle) {
  bool flag = false;
  FindSingularElements(angle);
  std::vector<std::set<int>> bndComposite;
  std::vector<int> bndId;
  for (auto bnd = m_bndComposite.begin(); bnd != m_bndComposite.end(); ++bnd) {
    bndId.push_back(bnd->first);
    bndComposite.push_back(bnd->second);
  }
  for (int i = 0; i < bndComposite.size() - 1; ++i) {
    for (int j = i + 1; j < bndComposite.size(); ++j) {
      int shared = 0, singluar = 0;
      if (FindSharedSingluarEdges(bndComposite[i], bndComposite[j], shared,
                                  singluar)) {
        std::cout << "error: bnd composits [" << bndId[i] << "] and ["
                  << bndId[j] << "] has shared faces\n";
      }
      if (shared > 0 && singluar == 0) {
        for (auto it = bndComposite[j].begin(); it != bndComposite[j].end();
             ++it) {
          bndComposite[i].insert(*it);
        }
        bndComposite[j].clear();
        flag = true;
      }
    }
  }
  m_bndComposite.clear();
  int index = 0;
  for (int i = 0; i < bndComposite.size(); ++i) {
    if (bndComposite[i].size() > 0) {
      m_bndComposite[index++] = bndComposite[i];
    }
  }
  return flag;
}

void MeshRegion::CheckMesh(double angle) {
  // check points
  std::map<int, int> ptsc;
  for (auto it = m_edges.begin(); it != m_edges.end(); ++it) {
    std::vector<int> e = it->second;
    ptsc[e[0]] += 1;
    ptsc[e[1]] += 1;
  }
  for (auto it = ptsc.begin(); it != ptsc.end(); ++it) {
    if (it->second == 1) {
      std::cout << "error: isolate point found (" << m_pts[it->first][0] << ","
                << m_pts[it->first][1] << "," << m_pts[it->first][2] << ")\n";
    }
  }
  // check edge
  std::map<int, int> edgec;
  char type;
  std::vector<int> edges;
  for (auto it = m_faces.begin(); it != m_faces.end(); ++it) {
    GetFaceEdges(it->first, edges, type);
    for (auto e = edges.begin(); e != edges.end(); ++e) {
      edgec[*e] += 1;
    }
  }
  for (auto it = edgec.begin(); it != edgec.end(); ++it) {
    if (it->second == 1) {
      std::cout << "error: isolate edge found edge[" << it->first << "]\n";
    }
  }
  // check boundary definition
  std::set<int> bndf;
  for (auto it = m_bndComposite.begin(); it != m_bndComposite.end(); ++it) {
    std::set<int> faces = it->second;
    for (auto f = faces.begin(); f != faces.end(); ++f) {
      bndf.insert(*f);
    }
  }
  std::map<int, int> facec;
  for (auto it = m_cells.begin(); it != m_cells.end(); ++it) {
    std::vector<int> faces = it->second;
    for (auto f = faces.begin(); f != faces.end(); ++f) {
      facec[*f] += 1;
    }
  }
  for (auto it = facec.begin(); it != facec.end(); ++it) {
    if (it->second == 1) {
      if (bndf.find(it->first) == bndf.end()) {
        std::cout << "error: undefined boundary face[" << it->first << "]\n";
      } else {
        bndf.erase(it->first);
      }
    }
  }
  for (auto it = bndf.begin(); it != bndf.end(); ++it) {
    std::cout << "error: boundary face [" << *it << "] is inner face\n";
  }
  std::cout << "mesh check finished " << m_name << std::endl;
}

void MeshRegion::OutPutSU2(std::string name) {
  std::ofstream su2(name.c_str());
  su2 << "NDIME= " << m_dim << std::endl;
  su2 << "NPOIN= " << m_pts.size() << std::endl;
  for (auto it = m_pts.begin(); it != m_pts.end(); ++it) {
    su2 << std::setprecision(12) << std::scientific;
    for (int i = 0; i < m_dim; ++i) {
      su2 << std::setw(20) << (it->second)[i] << " ";
    }
    su2 << std::endl;
  }
  std::map<int, std::vector<int>> elements;
  GetElements(elements);
  su2 << "NELEM= " << elements.size() << std::endl;
  for (auto it = elements.begin(); it != elements.end(); ++it) {
    std::vector<int> pts;
    char type;
    GetElementPts(it->first, pts, type);
    su2 << ElementSU2[ElementTypeMap[type]] << " ";
    for (int j = 0; j < pts.size(); ++j)
      su2 << pts[j] << " ";
    su2 << std::endl;
  }
  su2 << "NMARK= " << m_bndComposite.size() << std::endl;
  for (auto it = m_bndComposite.begin(); it != m_bndComposite.end(); ++it) {
    su2 << "MARKER_TAG= C" << it->first << std::endl;
    su2 << "MARKER_ELEMS= " << m_bndComposite[it->first].size() << std::endl;
    for (auto jt = m_bndComposite[it->first].begin();
         jt != m_bndComposite[it->first].end(); ++jt) {
      std::vector<int> pts;
      char type;
      GetBndElementPts(*jt, pts, type);
      su2 << ElementSU2[ElementTypeMap[type]] << " ";
      for (int j = 0; j < pts.size(); ++j)
        su2 << pts[j] << " ";
      su2 << std::endl;
    }
  }
  su2.close();
  std::cout << "write file " << name << std::endl;
}

int MeshRegion::PointHash(std::vector<double> p) {
  double dres = 0., nscatered = 100000.;
  for (int i = 0; i < m_dim; ++i) {
    double len = m_maxRange[i] - m_minRange[i];
    if (len < m_tolerance)
      len = 1.;
    dres = dres + nscatered * (p[i] - m_minRange[i]) / len;
  }
  return round(dres);
}

void MeshRegion::OutGmsh(std::string name) {
  std::ofstream ofile(name);
  ofile << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"; // fixed, real type is 8
  OutGmshNodes(ofile);
  OutGmshElements(ofile);
}

void MeshRegion::OutGmshNodes(std::ofstream &ofile) {
  ofile << "$Nodes\n";
  ofile << m_pts.size() << "\n";
  for (const auto &it : m_pts) {
    ofile << it.first << " " << it.second[0] << " " << it.second[1] << " "
          << it.second[2] << "\n";
  }
  ofile << "$EndNodes\n";
}

void MeshRegion::OutGmshElements(std::ofstream &ofile) {
  std::map<char, int> elmtype;
  for (size_t i = 0; i < ElementTag.size(); ++i) {
    elmtype[ElementTag[i]] = ElementGmsh[i];
  }
  ofile << "$Elements\n";
  ofile << m_cells.size() << "\n";
  for (auto it : m_cells) {
    std::vector<int> pts;
    char type;
    GetCellPts(it.first, pts, type);
    ofile << it.first << " " << elmtype[type] << " 2 0 0 ";
    for (auto jt : pts) {
      ofile << jt << " ";
    }
    ofile << '\n';
  }
  ofile << "$EndElements\n";
}
