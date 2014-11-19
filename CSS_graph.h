#include <fstream>
#include <string>
#include <set>


using namespace std;



bool read_graph(unordered_map< string, unordered_map<string, double> > &G,
                set<string> &node_set,
                string path,
                int is_undirected){
    
    bool   is_valid = true;

    string node1;
    string node2;
    double weight =  1.0;

    int    max_index = -1;


    vector<string> tok;

    string::size_type sz;


    ifstream f_in( path.c_str() );

    for (string l; getline(f_in, l); ){
        StringReplaceAll(l, "\t", " ");
        
        tok = StringSplit(trim(l), ' ');
        
        if ( !(tok.size() == 2 || tok.size() == 3) ){
            is_valid = false;
            break;
        }

        if (tok.size() == 3){
            weight = stod(tok[2], &sz);
        }else{
            weight = 1.0;
        }

        G[tok[0]][tok[1]] = weight;
        if (is_undirected){
            G[tok[1]][tok[0]] = weight;
        }
        
        node_set.insert(tok[0]);
        node_set.insert(tok[1]);
    }

    f_in.close();

    return is_valid;
}

