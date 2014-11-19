#include <unordered_map>


using namespace std;


void make_node_vec(unordered_map< int, unordered_map<int, double> > &G, 
                   vector<int> &node_vec){

    for (unordered_map< int, unordered_map<int, double> >::iterator it=G.begin(); it != G.end(); ++it){
        node_vec.push_back( (*it).first );
    }

    sort( node_vec.begin(), node_vec.end() );
}

/*
void make_edge_vec(unordered_map< int, unordered_map<int, double> > &G, 
                   vector< pair<int,int> > &edge_vec,
                   int is_undirected){

    int node1 = -1;
    int node2 = -1;

    unordered_map< int, double > nb_map;

    for (unordered_map< int, unordered_map<int, double> >::iterator it=G.begin(); it != G.end(); ++it){
        node1  = (*it).first;
        nb_map = (*it).second;
        
        for (unordered_map< int, double >::iterator it2=nb_map.begin(); it2 != nb_map.end(); ++it2){
            node2 = (*it2).first;

            if ( is_undirected && node1 >= node2 ){
                continue;
            }

            edge_vec.push_back( make_pair(node1, node2) );
        }
    }

}*/


void make_G_int_and_edge_vec(unordered_map< string, unordered_map<string, double> > &G,
                             unordered_map< int,    unordered_map<int,    double> > &G_int,
                             unordered_map< string, int > &name2idx,
                             vector< pair<int,int> > &edge_vec,
                             int is_undirected){

    string name1;
    string name2;

    int idx1;
    int idx2;

    unordered_map< string, double > nb_map;

    for (unordered_map< string, unordered_map<string, double> >::iterator it=G.begin(); it != G.end(); ++it){
        name1  = (*it).first;
        nb_map = (*it).second;
        
        for (unordered_map< string, double >::iterator it2=nb_map.begin(); it2 != nb_map.end(); ++it2){
            name2 = (*it2).first;

            idx1 = name2idx[name1];
            idx2 = name2idx[name2];

            G_int[idx1][idx2] = (*it2).second;

            if ( is_undirected && idx1 > idx2 ){
                continue;
            }

            edge_vec.push_back( make_pair(idx1, idx2) );

        }
    }
    
}
