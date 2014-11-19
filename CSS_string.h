#include <string>
#include <algorithm>
#include <functional>


using namespace std;

static inline string &ltrim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
    return s;
}

// trim from end
static inline string &rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline string &trim(string &s) {
    return ltrim(rtrim(s));
}


void StringReplaceAll(std::string& subject, const string& search, const string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
}

vector<string> StringSplit(string s, const char d){
    vector<string> e;
    string item;

    for (stringstream ss(s); getline(ss, item, d); ){
        e.push_back( item );
    }
    return e;
}

bool is_valid_int(const string &s){
    for (unsigned int i=0; i<s.length(); ++i){
        if (s.at(i) < 48 || s.at(i) > 57){
            return false;
        }
    }

    return true;
}
