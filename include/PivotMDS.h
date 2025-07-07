#ifndef _PIVOTMDS_H_
#define _PIVOTMDS_H_
#include "Node.h"
#include "Edge.h"
#include <vector>
#include <queue>
#include <string.h>
#include "randomkit.h"
#include <map>
#include <stack>
#include "Link.h"
#include <set>
#include <limits>
#include <cmath>
#include <algorithm>
#define Pair pair<int,int>
#define MAX 200000+100

   std::vector<std::vector<int> > build_graph_unweighted(int n, int m, int *I, int *J);

   void PivotMDS(std::vector<Node> &nodes, int k, int N, std::vector<std::vector<double>> &shortPat, std::vector<int> &pivot);

   void Dij(int V, int E, std::vector<Edge> e[MAX], int st, std::vector<float> &dist);

   std::vector<int> maxmin_random_sp_unweighted(const std::vector<std::vector<int> > &graph, int n_pivots, int p0, int seed);
   float Jaccard(std::set<int> set1,std::set<int> set2);
   void DFS2(int n,int m,std::vector<Edge> e[MAX],int st,std::vector<int> &kn, int neighbor);
   void maxmin_bfs_unweighted(const std::vector<std::vector<int> > &graph, const int p, std::vector<int> &mins, std::vector<int> &argmins);
   float norm(std::vector<float>& b_k, int k);
   void mult(std::vector<float>& b_k, std::vector<std::vector<float>>& B, int k);
   void DFS1(int n, int m, const std::vector<std::vector<int>> &e, int st, std::vector<int> &kn, int neighbor);
   void BFS(int V,int E,std::vector<Edge> e[MAX],int st,std::vector<int> &ind,std::vector<int>& count_com);
#endif