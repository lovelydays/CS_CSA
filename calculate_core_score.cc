
#include <cmath>
#include <cctype>

#include <fstream>
#include <sstream>
#include <iostream>

#include <string>
#include <utility>
#include <algorithm>
#include <functional>

#include <vector>
#include <unordered_map>
#include <unordered_set>


#include "CSS_string.h"
#include "CSS_graph.h"
#include "CSS_make_vec.h"
#include "CSS_CSA.h"


using namespace std;


void write_core_score(vector<string> &idx2name,
                      vector<double> &core_score_vec, 
                      string path){

    ofstream f_out( path.c_str() );

    for (int i=0; i<idx2name.size(); ++i){
        f_out << idx2name[i] << " " << core_score_vec[i] << endl;
    }

    f_out.close();

}



int main(int argc, char **argv){
    
    if (argc != 5){
        cout << endl;
        cout << "Usage: " << argv[0] << " <graph file> <d alpha> <d beta> <is undirected?> " << endl;
        cout << endl;
        cout << "   Node indices should be one-based. " << endl;
        cout << "   An undirected link should be represented as two edges. " << endl;
        cout << "       ex) 1 2 and 2 1 " <<endl;
        cout << endl;
        cout << "   Is a input graph undirected?" << endl;
        cout << "       Yes -> <is_undirected?> = 1" << endl;
        cout << "       No  -> <is_undirected?> = 1" << endl;
        cout << endl;
        cout << "   Core score out file path." << endl;
        cout << "       CoreScore-<graph file>.txt" << endl;
        cout << endl;
        cout << "   Core score out file format." << endl;
        cout << "       node 'i'<space>CS(i)<line break>" << endl;
        cout << endl;
        return 0;
    }

    string::size_type sz;
    
    string graph_path  = argv[1];
    string CS_out_path = "CoreScore-"+graph_path+".txt";
    
    double d_alpha = stod(argv[2], &sz);
    double d_beta  = stod(argv[3], &sz);

    int    is_undirected = stoi(argv[4], &sz);

    unordered_map< string, unordered_map<string, double> > G;
    unordered_map< int, unordered_map<int, double> > G_int;
    
    set< string > node_set;

    vector< pair<int,int> > edge_vec;

    // Read a graph and check its validity
    if ( !read_graph(G, node_set, graph_path, is_undirected) ){
        cout << endl;
        cout << "  Graph format is incorrect." << endl;
        cout << endl;
        return 0;
    }

    vector<string>             idx2name( node_set.begin(), node_set.end() );
    unordered_map<string, int> name2idx;

    for (int i=0; i<idx2name.size(); ++i){
        name2idx[ idx2name[i] ] = i;
    }

    vector<int> node_vec(idx2name.size());

    for (int i=0; i<node_vec.size(); ++i){
        node_vec[i] = i+1;
    }

    make_G_int_and_edge_vec(G, G_int, name2idx, edge_vec, is_undirected);


    vector<double> core_score_vec;

    CSA_iter(G_int, node_vec, edge_vec, d_alpha, d_beta, is_undirected, core_score_vec);


    write_core_score(idx2name, core_score_vec, CS_out_path);
    

    return 1;
}
