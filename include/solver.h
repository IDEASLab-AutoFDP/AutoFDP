
#ifndef SOLVER_H
#define SOLVER_H

#include "time.h"
#include <sys/time.h>
#include <iomanip>
#pragma warning(disable: 4702)
#pragma warning(pop)
#include <ctime>
#include <iterator>
#include <numeric>
#include <random>
#include <regex>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>
#include "metric.h"
#include <algorithm>
#include <cassert>
#include "PivotMDS.h"
#include <omp.h> 
#include <time.h>
#include <sstream>
#include <set>
#include "shuffle.h"
#include "mkpath.h"


#define randf() (1.0f*random()/RAND_MAX)


#define INF 0x3f3f3f3f

#define DEFAULT_WEIGHT 1
#define velocityDecay  (0.8f)

int str2int(std::string s);
float str2float(std::string s);

class graph {
    public:
        int n, m;
        float k1,k2,a1,a2,b1,b2;
        
        std::vector<Node> nodes;
        std::vector<Link> links;

        std::vector<int> pivots;
        std::vector<int> flag;

        float loss_scale;
        std::vector<std::vector<double>> sp;

        std::vector<std::pair<int,int>> constraints_sgd;

        std::vector<float> schedule;

        rk_state rstate;


        std::vector<std::vector<double>> P;

        graph();

        int get_n(){return n;}
        void set_n(int _n){n=_n;}
        int get_m(){return m;}
        void set_m(int _m){n=_m;}
        void get_degree();
        

        void initgraph(std::string filename);
        void initRandomPosition();
        void solveDij();

        void preSolveSGD(double eps = 0.1, int t_max = 200, int seed = 42, float eta_max = 100.0f);
        bool solve_bilevel(int iter);
        std::vector<Node> solve_k_with_noupdate(int iter,float Delta_a1=0.0f,float Delta_b1=0.0f,float Delta_a2=0.0f,float Delta_b2=0.0f,float Delta_k1=0.0f,float Delta_k2=0.0f);
    
        float opti_SE(int iter,double Delta=0.01f,double k=1.0f);
        float opti_IL(int iter,double Delta=0.01f,double k=1.0f);
        float opti_IL_mult(int iter, double Delta, double k=1.0f);
        float opti_NP(int iter,double Delta=0.01f,int neigh=2,double k=1.0f);
        float opti_CL(int iter,double Delta=0.01f,double k=1.0f);
        float opti_CA(int iter,double Delta=0.01f,double k=1.0f);
        float opti_AR(int iter,double Delta=0.01f,double k=1.0f);
        float opti_MA(int iter,double Delta=0.01f,double k=1.0f);
        float opti_NR(int iter,double Delta=0.01f,double k=1.0f);
        float opti_GP(int iter,double Delta=0.01f,double k=1.0f);

        void calP_DeepGD();
        float opti_MA_DeepGD(int iter,double Delta=0.01f,double k=1.0f);
        float opti_NR_DeepGD(int iter, double Delta=0.01f,double k=1.0f);
        float opti_tSNE_DeepGD(int iter, double Delta, double k);

        void drawPathSVG(std::string filename, std::string path,bool output_pos=false);
    };

#endif
