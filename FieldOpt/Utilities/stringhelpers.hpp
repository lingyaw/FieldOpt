/// This file contains helper methods for working with time.
#ifndef STRINGHELPERS_FUNCTIONS_H
#define STRINGHELPERS_FUNCTIONS_H

#include <string>
#include <vector>

using namespace std;

template <typename T>
inline string vec_to_str(vector<T> vec) {
    string str = "";
    if (vec.size() == 0) return str;
    str = to_string(vec[0]);
    if (vec.size() > 1) {
        for (int i = 1; i < vec.size(); ++i) {
            str = str + ", " + to_string(vec[i]);
        }
    }
    return str;
}

#endif // STRINGHELPERS_FUNCTIONS_H
