#include <cmath>
#include <utility>

#include <vector>
#include <unordered_map>


using namespace std;

void calculate_C(vector<int> &node_vec, 
                 double alpha,
                 double beta,
                 vector<double> &C){

    C.resize(node_vec.size());
    
    double N           = (double)( (int)node_vec.size() );
    double one_alpha   = 1.0 - alpha;
    double floor_betaN = floor( beta*N );

    double c1 = one_alpha / ( 2.0*floor_betaN );
    
    double c2 = (one_alpha) / ( 2.0*(N - floor_betaN) );
    double c3 = ( 1.0 + alpha ) * 0.5;

    for (int t=0; t<N; ++t){
        int i = node_vec[t]; // 0 -> 1

        C[ t ] = (i<=floor_betaN) ? i * c1 : (i - floor_betaN)*c2 + c3;
    }

}


double calculate_R(double alpha, 
                   double beta, 
                   unordered_map< int, unordered_map<int, double> > &G,
                   vector<double> &C,
                   vector< int > &node_vec,
                   vector< pair<int,int> > &edge_vec,
                   int is_undirected){

    double R = 0.0;
    double d = 1.0;

    if (is_undirected){
        d = 2.0;
    }

    for (int t=0; t<edge_vec.size(); ++t){
        int i = edge_vec[t].first;
        int j = edge_vec[t].second;

        R += G[i][j] * C[i] * C[j] * d;
    }

    return R;
}

/*
void calculate_CS(double d_alpha,
                  double d_beta,
                  unordered_map< int, unordered_map<int, double> > &G,
                  vector< int > &node_vec,
                  vector< pair<int,int> > &edge_vec,
                  unordered_map<int, double> &CS){

    for (double alpha=0.0; alpha<=1.0; alpha += d_alpha){
        for (double beta=0.0; beta<1.0; beta += d_beta){
            vector<double> C;

            calculate_C(node_vec, alpha, beta, C);

            double R = calculate_R(alpha, beta, G, C, node_vec, edge_vec, is_undirected);

            for (int t=0; t<node_vec.size(); ++t){
                int i = node_vec[t];
                CS[i] += C[t]*R;
            }
        }
    }

}
*/

