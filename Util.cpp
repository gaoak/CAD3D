#include"Util.h"
#include<iostream>
#include<cmath>
std::string printComposite(std::set<int> value) {
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

void parserDouble(const char * cstr, std::vector<double> & value) {
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

void parserUInt(const char * cstr, std::vector<int> & value) {
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

double transformz(double z, int nlayers, std::vector<double> &targz){
    int i = round(z*nlayers);
    if(targz.size()>i) {
        return targz[i];
    } else {
        return z;
    }
}

void transform(std::vector<double> &p, double AoA) {
    double x = p[0], y = p[1];
    p[0] = x*cos(AoA) + y*sin(AoA);
    p[1] =-x*sin(AoA) + y*cos(AoA);
}
