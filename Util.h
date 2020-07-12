#ifndef UTIL_H
#define UTIL_H
#include<vector>
#include<string>
#include<set>

double transformz(double z, int nlayers, std::vector<double> &targz);
void parserUInt(const char * cstr, std::vector<int> & value);
void parserDouble(const char * cstr, std::vector<double> & value);
template<typename T>
void printVector(char * str, const char * format, std::vector<T> values) {
    int pos = 0;
    for(int i=0; i<values.size(); ++i) {
        int len = sprintf(str+pos, format, values[i]);
        pos += len;
    }
}
std::string printComposite(std::set<int> value);
void transform(std::vector<double> &p, double AoA);
#endif
