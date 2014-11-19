
#include <cmath>
#include <chrono>
#include <iostream>
#include <algorithm>

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <omp.h>

#include "CSS_calculate.h"

#include "dSFMT-2.2.3/dSFMT.c"


using namespace std;



using namespace std;

    
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}


void make_double_range01(vector<double> &v, double d){
    for (double a=0.0; a <= 1.0; a += d){
        v.push_back(a);
    }
}

inline double get_rnd_01co_real(){
    return dsfmt_gv_genrand_close_open();
}

inline int get_rnd_0V(int V){
    return (int)(dsfmt_gv_genrand_close_open()*V);
}


inline pair<int,int> get_rnd_int_pair_0V(int V){
    int i1 = (int)(dsfmt_gv_genrand_close_open()*V);
    int i2 = (int)(dsfmt_gv_genrand_close_open()*V);
    while (i1 == i2){
        i2 = (int)(dsfmt_gv_genrand_close_open()*V);
    }
    return make_pair(i1, i2);
}


void sample_int_pairs_range_0V(set< pair<int,int> > &sample_set, int V1, int V2, int N_rnd){
    int d1;
    int d2;

    while (sample_set.size() != N_rnd){
        d1 = get_rnd_0V(V1);
        d2 = get_rnd_0V(V2);

        sample_set.insert( make_pair(d1, d2) );
    }
}


void shuffle_node_vec(vector<int> &node_vec2){

    random_shuffle( node_vec2.begin(), node_vec2.end(), get_rnd_0V);
}


void fill_idx2node_vec(vector<int> &node_vec, vector< vector<int> > &idx2node_vec){

    idx2node_vec[0] = node_vec;

    shuffle_node_vec(idx2node_vec[0]);

    for (int i=1; i<idx2node_vec.size(); ++i){
        idx2node_vec[i] = idx2node_vec[i-1];

        shuffle_node_vec(idx2node_vec[i]);
    }
}


void shuffle_node_vec_nb(vector<int> &node_vec_nb, int idx1, int idx2){
    int t = -1;

    t = node_vec_nb[idx1];
    node_vec_nb[idx1] = node_vec_nb[idx2];
    node_vec_nb[idx2] = t;
}


void fill_idx2node_vec_nb(vector<int>               &node_vec_t, 
                          vector< vector<int> >     &idx2node_vec_nb,
                          int idx1,
                          unordered_map<int,double> &t_nb){

    vector<int> nb_vec;

    for (unordered_map<int,double>::iterator it=t_nb.begin(); it != t_nb.end(); ++it){
        nb_vec.push_back( (*it).first );
    }

    for (int i=0; i<nb_vec.size(); ++i){
        idx2node_vec_nb[i] = node_vec_t;
        
        shuffle_node_vec_nb(idx2node_vec_nb[i], idx1, nb_vec[i]);
    }

}


void calculate_R_iter(unordered_map< int, unordered_map<int, double> > &G,
                      vector< pair<int,int> > &edge_vec,
                      const double alpha,
                      const double beta,
                      vector< vector<int> > &idx2node_vec,
                      vector<double> &idx2R,
                      vector< vector<double> > &idx2C,
                      int is_undirected){
    
    for (int i=0; i<idx2node_vec.size(); ++i){
        vector<double> C;
        calculate_C(idx2node_vec[i], alpha, beta, C);

        idx2R[i] = calculate_R(alpha, beta, G, C, idx2node_vec[i], edge_vec, is_undirected);

        idx2C[i] = C;
    }

}


void local_optimization(unordered_map< int, unordered_map<int, double> > &G_int,
                        vector< pair<int,int> > &edge_vec,
                        const double alpha,
                        const double beta,
                        vector< vector<int> > &idx2node_vec,
                        vector<double> &idx2R,
                        vector< vector<double> > &idx2C,
                        int is_undirected){
    for (int t=0; t<idx2node_vec.size(); ++t){

        for (int t2=0; t2<(int)(idx2node_vec.size()*1); ++t2){
            int rnd_idx = get_rnd_0V( idx2node_vec[t].size() );
            vector< vector<int> > idx2node_vec_nb;
            idx2node_vec_nb.resize( G_int[rnd_idx].size() );
            
            fill_idx2node_vec_nb(idx2node_vec[t], idx2node_vec_nb, rnd_idx, G_int[rnd_idx]);

            vector< double >         idx2R_nb( idx2node_vec_nb.size() );
            vector< vector<double> > idx2C_nb( idx2node_vec_nb.size() );
            
            calculate_R_iter(G_int, edge_vec, alpha, beta, idx2node_vec_nb, idx2R_nb, idx2C_nb, is_undirected);

            // Find max R index
            double tmp_R     = idx2R[t];
            int    tmp_R_idx = -1;

            for (int t2=0; t2<idx2R_nb.size(); ++t2){
                if (idx2R_nb[t2] > tmp_R){
                    tmp_R     = idx2R_nb[t2];
                    tmp_R_idx = t2;
                }
            }

            if (tmp_R > idx2R[t]){
                idx2R[t] = tmp_R;
                idx2C[t] = idx2C_nb[tmp_R_idx];

                idx2node_vec[t] = idx2node_vec_nb[tmp_R_idx];
            }else{
                break;
            }
        }
    }
}



double calculate_distance_between_core_vectors(vector<double> &v1, vector<double> &v2){
    double sol_dist = 0.0;
    for (int i=0; i<v1.size(); ++i){
        sol_dist += (v1[i]-v2[i])*(v1[i]-v2[i]);
    }
    return sol_dist;
}

double calculate_distance_between_core_vectors_iter(vector< vector<double> > &idx2C,
                                                    unordered_map< int, unordered_map<int, double> > &idx2idx2sol_dist){
    double D_avr;
    
    int N = idx2C.size();
    
    for (int i=0; i<N-1; ++i){
        for (int j=i+1; j<N; ++j){
            double sol_dist = calculate_distance_between_core_vectors(idx2C[i], idx2C[j]);
            
            idx2idx2sol_dist[i][j] = sol_dist;

            D_avr += sol_dist;
        }
    }
    D_avr = D_avr / (double)(N*(N-1))*2;

    return D_avr;
}


int sum_seed_flags(vector<int> &seed_flags){
    int s = 0;
    for (int i=0; i<seed_flags.size(); ++i){
        s += seed_flags[i];
    }
    return s;
}


void make_pair_set(vector<double> &idx2R, set< pair<double,int> > &s){
    for (int i=0; i<idx2R.size(); ++i){
        s.insert( make_pair(idx2R[i], i) );
    }
}


void select_solutions(vector<double> &idx2R,
                      set< pair<double, int> > &R_idx_set,
                      int N_s){

    int idx;

    while ( R_idx_set.size() < N_s ){
        idx = get_rnd_0V( idx2R.size() );

        R_idx_set.insert( make_pair(idx2R[idx], idx) );
    }
}


void crossover(vector< vector<int> > &idx2node_vec,
               vector< pair<double, int> > &R_idx_vec,
               set< vector<int> > &T_c,
               int N_c){

    int idx;
    pair<double,int> tmp_pair;

    while (T_c.size() < N_c){
        idx = get_rnd_0V( R_idx_vec.size() );

        vector<int> tmp_node_vec = idx2node_vec[idx];

        rotate( tmp_node_vec.begin(), tmp_node_vec.begin()+(int)(tmp_node_vec.size()*0.2), tmp_node_vec.end() );
        //rotate( tmp_node_vec.begin(), tmp_node_vec.begin()+1, tmp_node_vec.end() );
        
        T_c.insert( tmp_node_vec );
    }
}


void mutation(vector< vector<int> > &idx2node_vec,
              vector< pair<double, int> > &R_idx_vec,
              set< vector<int> > &T_m,
              int N_m){

    int idx1, idx2;
    int rnd_idx;
    int t;

    pair<int,int> idx_pair;

    while ( T_m.size() < N_m ){
        idx1 = get_rnd_0V( R_idx_vec.size() );

        idx2 = R_idx_vec[idx1].second;

        vector<int> tmp_node_vec = idx2node_vec[idx2];

        for (int i=0; i<(int)(tmp_node_vec.size()*0.1); ++i){
        //for (int i=0; i<1; ++i){
            idx_pair = get_rnd_int_pair_0V( tmp_node_vec.size() );

            t = tmp_node_vec[idx_pair.first];
            tmp_node_vec[idx_pair.first]  = tmp_node_vec[idx_pair.second];
            tmp_node_vec[idx_pair.second] = t;
        }
        T_m.insert( tmp_node_vec );
    }
}


void find_nearest_solution_in_P(vector<double> &w,
                                vector< vector<double> > &idx2C,
                                set< pair<double, int> > &nearest_sols){

    for (int i=0; i<idx2C.size(); ++i){
        nearest_sols.insert( make_pair(calculate_distance_between_core_vectors(w, idx2C[i]), i) );
    }
}


void replace_solution(int t,
                      pair<double, int> w,
                      vector< vector<int> > &idx2node_vec,
                      vector<double> &idx2R,
                      vector< vector<double> > &idx2C,
                      vector< vector<int> > &Tc_Tm,
                      vector<double> &Tc_Tm2R,
                      vector< vector<double> > &Tc_Tm2C,
                      unordered_map< int, unordered_map<int, double> > idx2idx2sol_dist,
                      set< pair<double,int> > &nearest_sols){

    int idx1, idx2;
    
    // Replace solution distances
    vector< pair<double, int> > nearest_sols_vec( nearest_sols.begin(), nearest_sols.end() );

    for (int i=1; i<nearest_sols_vec.size(); ++i){
        idx1 = w.second;
        idx2 = nearest_sols_vec[i].second;

        if (idx1 < idx2){
            idx2idx2sol_dist[idx1][idx2] = w.first;
        }else{
            idx2idx2sol_dist[idx2][idx1] = w.first;
        }
    }

    // Replace idx2R
    idx2R[w.second] = Tc_Tm2R[t];

    // Replace idx2C
    idx2C[w.second] = Tc_Tm2C[t];

    // Replace idx_pair_vec
    idx2node_vec[w.second] = Tc_Tm[t];
}


void CSA(unordered_map< int, unordered_map<int, double> > &G_int,
         vector< int > &node_vec,
         vector< pair<int,int> > &edge_vec,
         const double alpha, 
         const double beta,
         const int    is_undirected,
         vector<double> &max_idx2C,
         double         &max_idx2R){
    
    const double CSA_alpha = 0.995;
    
    const int N_ran_sols = 50;
    const int N_s        = 30;
    const int N_c        = 10;
    const int N_m        = 5;


    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();

    dsfmt_gv_init_gen_rand(seed);

    // Initialize the first bank, P, with N_ran_sols
    //     by randomize node_vec
    vector< vector<int> > idx2node_vec;
    idx2node_vec.resize(N_ran_sols);

    fill_idx2node_vec(node_vec, idx2node_vec);

    
    vector< double >         idx2R( idx2node_vec.size() );
    vector< vector<double> > idx2C( idx2node_vec.size() );

    for (int out_loop=0; out_loop<2; ++out_loop){
            
        //     Calculate R for each pair
        calculate_R_iter(G_int, edge_vec, alpha, beta, idx2node_vec, idx2R, idx2C, is_undirected);   

        // Optimize P by a local maximization
        local_optimization(G_int, edge_vec, alpha, beta, idx2node_vec, idx2R, idx2C, is_undirected);


        // Initialize seed flags
        vector<int> seed_flags( idx2R.size(), 0 );

        // Get D_avr
        unordered_map< int, unordered_map<int, double> > idx2idx2sol_dist;
            
        double D_avr = calculate_distance_between_core_vectors_iter(idx2C, idx2idx2sol_dist);
        double D_cut = D_avr * 0.5;

        
        while ( sum_seed_flags(seed_flags) < seed_flags.size() ){
            
            set< pair<double, int> > R_idx_set;

            select_solutions(idx2R, R_idx_set, N_s);
            
            vector< pair<double, int> > R_idx_vec( R_idx_set.begin(), R_idx_set.end() );

            for (int i=0; i<R_idx_vec.size(); ++i){
                seed_flags[R_idx_vec[i].second] = 1;
            }

            // Crossover
            set< vector<int> > T_c;

            crossover(idx2node_vec, R_idx_vec, T_c, N_c);

            // Mutation
            set< vector<int> > T_m;

            mutation(idx2node_vec, R_idx_vec, T_m, N_m);


            set< vector<int> > Tc_Tm_set;

            Tc_Tm_set.insert(T_c.begin(), T_c.end());
            Tc_Tm_set.insert(T_m.begin(), T_m.end());


            vector< vector<int> > Tc_Tm( Tc_Tm_set.begin(), Tc_Tm_set.end() );
            vector< vector<int> >::iterator it_Tc_Tm;


            vector< double > Tc_Tm2R( Tc_Tm.size() );
            vector< vector<double> > Tc_Tm2C( Tc_Tm.size() );

            calculate_R_iter(G_int, edge_vec, alpha, beta, Tc_Tm, Tc_Tm2R, Tc_Tm2C, is_undirected);
            local_optimization(G_int, edge_vec, alpha, beta, Tc_Tm, Tc_Tm2R, Tc_Tm2C, is_undirected);

            pair<double, int> u;
            pair<double, int> w;

            for (int t=0; t<Tc_Tm.size(); ++t){
                R_idx_set.clear();

                make_pair_set(idx2R, R_idx_set);

                u = *(R_idx_set.begin());

                if ( Tc_Tm2R[t] < idx2R[u.second] ){
                    continue;
                }
                
                set< pair<double,int> > nearest_sols;

                find_nearest_solution_in_P(Tc_Tm2C[t], idx2C, nearest_sols);
                
                w = *(nearest_sols.begin());

                if ( w.first > D_cut ){
                    replace_solution(t, u, idx2node_vec, idx2R, idx2C, Tc_Tm, Tc_Tm2R, Tc_Tm2C, idx2idx2sol_dist, nearest_sols);   
                    seed_flags[u.second] = 0;
                }else if( Tc_Tm2R[t] > idx2R[w.second] ){
                    replace_solution(t, w, idx2node_vec, idx2R, idx2C, Tc_Tm, Tc_Tm2R, Tc_Tm2C, idx2idx2sol_dist, nearest_sols);   
                    seed_flags[w.second] = 0;
                }
            }
            D_cut = max( CSA_alpha*D_cut, D_avr*0.2 );
            
            if (D_cut < 0.01){
                break;
            }
            //cout << D_cut << " " << sum_seed_flags(seed_flags) << endl;
        }
    }

    set< pair<double, int> > tmp_set;

    for (int i=0; i<idx2R.size(); ++i){
        tmp_set.insert( make_pair(idx2R[i], i) );
    }

    vector< pair<double, int> > tmp_vec(tmp_set.begin(), tmp_set.end());

    map< int, double > core_score;

    max_idx2C = idx2C[tmp_vec[tmp_vec.size()-1].second];
    max_idx2R = idx2R[tmp_vec[tmp_vec.size()-1].second];
    /*
    for (int j=0; j<idx2C[tmp_vec[tmp_vec.size()-1].second].size(); ++j){
        core_score[j] += idx2C[tmp_vec[tmp_vec.size()-1].second][j]*idx2R[tmp_vec[tmp_vec.size()-1].second];
    }*/

}


void CSA_iter(unordered_map< int, unordered_map<int, double> > &G_int, 
              vector<int> &node_vec,
              vector< pair<int,int> > &edge_vec,
              const double d_alpha, 
              const double d_beta,
              int is_undirected,
              vector<double> &core_score_vec){

    vector<double> range_alpha;
    vector<double> range_beta;

    // Initialize range vectors
    make_double_range01( range_alpha, d_alpha );
    make_double_range01( range_beta , d_beta  );

    
    // Initialize SFMT random number generator
    
    //double alpha = range_alpha[get_rnd_0V( range_alpha.size() )];
    //double beta  = range_beta[get_rnd_0V( range_beta.size() )];

    vector<double> core_score( node_vec.size(), 0.0 );
    /*
    double *core_score = new double[node_vec.size()];
    for (int i=0; i<node_vec.size(); ++i){
        core_score[i] = 0.0;
    }
    */

    //omp_set_num_threads(2);

//#pragma omp parallel for shared(core_score)
    for (int t1=0; t1<range_alpha.size(); ++t1){
        for (int t2=0; t2<range_beta.size(); ++t2){
            double alpha = range_alpha[t1];
            double beta  = range_beta[t2];
            
            vector<double> max_idx2C;
            double         max_idx2R;

            CSA(G_int, node_vec, edge_vec, alpha, beta, is_undirected, max_idx2C, max_idx2R);
            
            for (int i=0; i<max_idx2C.size(); ++i){
                core_score[i] += max_idx2C[i]*max_idx2R;
            }
            cout << "alpha=" << alpha << ", beta=" << beta << endl;
        }
    }
    
    double mn = ( *min_element( core_score.begin(), core_score.end() ) );
    double mx = ( *max_element( core_score.begin(), core_score.end() ) );

    double mn_inv = 1.0/mn;
    double mx_inv = 1.0/mx;

    double tmp1 = (mx-mn);
    double tmp2 = 1.0/tmp1;
    double tmp3 = mn*tmp2;

    for (int i=0; i<node_vec.size(); ++i){
        //core_score[i] = core_score[i]*tmp2-tmp3;
        core_score[i] = core_score[i]*mx_inv;
        
        if (core_score[i] < 0 || core_score[i] < 1E-12){ 
            core_score[i] = 0; 
        }
        
        core_score_vec.push_back( core_score[i] );
    }

    //delete [] core_score;
}
