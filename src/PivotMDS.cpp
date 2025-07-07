#include "PivotMDS.h"
void Dij(int V,int E,std::vector<Edge> e[MAX],int st,std::vector<float> &dist){
    std::priority_queue<std::Pair,std::vector<std::Pair>,std::greater<std::Pair>> q;
    int visited[V];
    memset(visited,0,sizeof(visited));
    for(int i=0;i<V;i++) dist[i]=MAX;
    dist[st]=0;
    q.push(std::make_pair(0,st));
    while(!q.empty()){
        std::Pair t=q.top();
        q.pop();
        if(visited[t.second]) continue;
        visited[t.second]=1;

        for(int i=0;i<e[t.second].size();i++){
            int son=e[t.second][i].end;
            if(dist[son]>dist[t.second]+e[t.second][i].w){
                dist[son]=dist[t.second]+e[t.second][i].w;
                if(!visited[son])
                    q.push(std::make_pair(dist[son],son));
            }
        }
    }
}
void maxmin_bfs_unweighted(const std::vector<std::vector<int> >& graph, const int p, std::vector<int>& mins, std::vector<int>& argmins)
{
    int n = graph.size();
    std::queue<int> q;
    std::vector<int> d(n, -1);

    q.push(p);
    d[p] = 0;
    while (!q.empty())
    {
        int current = q.front();
        q.pop();

        for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
        {
            const int &next = graph[current][i_edge];
            if (d[next] == -1)
            {
                q.push(next);
                d[next] = d[current] + 1;
                if (d[next] < mins[next])
                {
                    mins[next] = d[next];
                    argmins[next] = p;
                }
            }
        }
    }
}
std::vector<std::vector<int> > build_graph_unweighted(int n, int m, int* I, int* J)
{
    // used to make graph undirected, in case it is not already
    std::vector<std::set<int>> undirected(n);
    std::vector<std::vector<int> > graph(n);

    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw std::invalid_argument("i or j bigger than n");

        if (i != j && undirected[j].find(i) == undirected[j].end()) // if edge not seen
        {
            undirected[i].insert(j);
            undirected[j].insert(i);
            graph[i].push_back(j);
            graph[j].push_back(i);
        }
    }
    return graph;
}
std::vector<int> maxmin_random_sp_unweighted(const std::vector<std::vector<int> >& graph, int n_pivots, int p0, int seed)
{
    int n = graph.size();

    std::vector<int> mins(n, std::numeric_limits<int>::max());
    std::vector<int> argmins(n, -1);

    // first pivot
    mins[p0] = 0;
    argmins[p0] = p0;
    maxmin_bfs_unweighted(graph, p0, mins, argmins);
    for (int i = 0; i < n; i++)
    {
        if (argmins[i] == -1)
            throw std::invalid_argument("graph has multiple connected components");
    }

    // remaining pivots
    // std::mt19937 rng(seed);
    // std::uniform_real_distribution<> uniform(0, 1);
    rk_state rstate;
    rk_seed(seed, &rstate);
    for (int ij = 1; ij < n_pivots; ij++)
    {
        // choose pivots with probability min
        int min_total = 0;
        for (int i = 0; i < n; i++)
        {
            min_total += mins[i];
        }
        int sample = rk_interval(min_total, &rstate);
        int cumul = 0;
        int argmax = -1;
        for (int i = 0; i < n; i++)
        {
            cumul += mins[i];
            if (cumul >= sample)
            {
                argmax = i;
                break;
            }
        }
        if (argmax == -1)
            throw std::invalid_argument("unweighted pivot sampling failed");

        mins[argmax] = 0;
        argmins[argmax] = argmax;
        maxmin_bfs_unweighted(graph, argmax, mins, argmins);
    }
    return argmins;
}

void DFS1(int n, int m, const std::vector<std::vector<int>> &e, int st, std::vector<int> &kn, int neighbor) {
    bool visited[n];
    int layer[n];
    memset(layer,0,n);
    memset(visited,0,n);
    visited[st]=1;
    std::queue<int> q;
    q.push(st);
    layer[st]=0;
    kn.pop_back();
    while(!q.empty()){
        int u=q.front();
        q.pop();
        for(int i=0;i<e[u].size();i++){
            int v=e[u][i];
            if(!visited[v]) {
                layer[v]=layer[u]+1;
                visited[v]=1;
                if(layer[v]<=neighbor){
                    q.push(v);
                    kn.push_back(v);
                }
            }
        }
    }
}

void BFS(int V,int E,std::vector<Edge> e[MAX],int st,std::vector<int> &ind,std::vector<int>& count_com){
    //int visited[V];
    int *visited=new int[V];
    memset(visited,0,sizeof(visited));
    std::queue<int> q1;
    q1.push(st);
    //count_com[ind[st]]++;
    while(!q1.empty()){
        int u=q1.front();
        q1.pop();
        if(visited[u]) continue;
        visited[u]=1;
        ind[u]=ind[st];
        //cout<<"indu="<<ind[u]<<"  indst"<<ind[st]<<endl;
        count_com[ind[u]]++;
        for(int i=0;i<e[u].size();i++){
            int v=e[u][i].end;
            if(!visited[v]){
                q1.push(v);

            }
        }
    }
}

float Jaccard(std::set<int> set1,std::set<int> set2){
    std::set<int> inter_set;
    std::set<int> union_set;
    //sort(set1.begin(),set1.end());
    //sort(set2.begin(),set2.end());

    std::set_intersection(set1.begin(),set1.end(),set2.begin(),set2.end(), inserter(inter_set,inter_set.begin()));
    std::set_union(set1.begin(),set1.end(),set2.begin(),set2.end(), inserter(union_set,union_set.begin()));
    return (float)inter_set.size()/(float)union_set.size();
}

void DFS2(int n,int m,std::vector<Edge> e[MAX],int st,std::vector<int> &kn, int neighbor){
//void DFS_K(int n, int m, const vector<vector<int>> &e, int st, vector<int> &kn, int neighbor) {
    bool visited[n];
    int layer[n];
    memset(layer,0,n);
    memset(visited,0,n);
    visited[st]=1;
    std::queue<int> q;
    q.push(st);
    layer[st]=0;
    kn.pop_back();
    while(!q.empty()){
        int u=q.front();
        q.pop();
        for(int i=0;i<e[u].size();i++){
            int v=e[u][i].end;
            if(!visited[v]) {
                layer[v]=layer[u]+1;
                visited[v]=1;
                if(layer[v]<=neighbor){
                    q.push(v);
                    kn.push_back(v);
                }
            }
        }

    }
}

void mult(std::vector<float>& b_k, std::vector<std::vector<float>>& B, int k) {
    std::vector<float> b_k1(k);
    for (int i = 0; i < k; i++) {
        b_k1[i] = b_k[i];
        b_k[i] = 0;
    }
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            b_k[i] += B[i][j] * b_k1[j];
        }
    }
}
float norm(std::vector<float>& b_k, int k) {
    float a = 0;
    for (int i = 0; i < k; i++) {
        a += (b_k[i] * b_k[i]);
    }
    return a;
}
void PivotMDS(std::vector<Node> &nodes,int  k,int N,std::vector<std::vector<double>> &shortPat,std::vector<int> &pivot){
    int iteration = 100;
    std::vector<std::vector<float>> C(N, std::vector<float>(k));
    std::vector<std::vector<int>> dkj(N, std::vector<int>(k));
    //select k pivots
    //rand
    //set ShortestPath
    //calculate B
    //calculate CC^T
    //calculate first second

    //Init dik
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            dkj[i][j] = shortPat[i][pivot[j]];

        }
    }
    float sum = 0;
    std::vector<float> sumn(N);
    std::vector<float> sumk(k);
    for (int i = 0; i < N; i++) {
        sumn[i] = 0;
        // sumk[i] = 0;
    }
    for (int i = 0; i < k; i++) {
        // sumn[i] = 0;
        sumk[i] = 0;
    }
    //calculate Cij
    for (int i = 0; i < N; i++) {
        //sumk[i] = 0;
        for (int j = 0; j < k; j++) {
            //dkj[i][j]=1;//shortPat[pivot[i]][j]
            sum += (dkj[i][j] * dkj[i][j]);
            sumn[i] += (dkj[i][j] * dkj[i][j]);
        }
        sumn[i] /= k;
    }
    sum /= (N * k);

    for (int j = 0; j < k; j++) {
        //sumn[j] = 0;
        for (int i = 0; i < N; i++) {
            //dkj[i][j]=1;//shortPat[pivot[i]][j]
            //sum+=(dkj[i][j]*dkj[i][j]);
            sumk[j] += (dkj[i][j] * dkj[i][j]);
        }
        sumk[j] /= N;
    }
    //calculate C
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            C[i][j] = -0.5f * ((dkj[i][j] * dkj[i][j]) - sumk[j] - sumn[i] + sum);
        }
    }
    std::vector<std::vector<float>> B(k, std::vector<float>(k));
    // B.resize(N,vector<float>(N));

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            for (int _ = 0; _ < N; _++) {
                B[i][j] += C[_][i] * C[_][j];
            }
        }
    }
    std::vector<float> b_k(k);
    for (int i = 0; i < k; i++) {
        if (i % 2 == 0) {
            b_k[i] = 1;
        }
        else {
            b_k[i] = -1;
        }
    }
    //?????????????
    for (int i = 0; i < iteration; i++) {
        mult(b_k, B, k);
        //float norm1=norm(b_k,k);
        float a = 0;
        for (int i = 0; i < k; i++) {
            a += (b_k[i] * b_k[i]);
        }
        for (int i = 0; i < k; i++) {
            b_k[i] /=sqrtf(a);
        }
    }
    //?????
    float normV = norm(b_k, k);
    float a = 0;
    for (int i = 0; i < k; i++) {
        a += (b_k[i] * b_k[i]);
    }
    for (int i = 0; i < k; i++) {
        b_k[i] /= sqrtf(a);
    }
    //???????
    float lamdab = 0;
    std::vector<float> temp(k);
    for (int i = 0; i < k; i++) {
        temp[i] = b_k[i];
    }
    mult(temp, B, k);

    for (int i = 0; i < k; i++) {
        lamdab += temp[i] * b_k[i];
        // cout << "temp" << i << "=   b_k="<<b_k[i] << " ";
    }
    std::vector<std::vector<float>> A(k, std::vector<float>(k));
    std::vector<float> a_k(k);
    float lamdaa = 0;
    // float normV = norm(b_k, k);

    //??A
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            A[i][j] = B[i][j] - lamdab * b_k[i] * b_k[j] / normV;
        }
    }
    //?????A_k
    for (int i = 0; i < k; i++) {
        if (i % 2 == 0) {
            a_k[i] = 1;
        }
        else {
            a_k[i] = -1;
        }
    }

    //??��?????????
    for (int i = 0; i < iteration; i++) {
        mult(a_k, A, k);
        //float norm1=norm(b_k,k);
        float a2 = 0;
        for (int i = 0; i < k; i++) {
            a2 += (a_k[i] * a_k[i]);
        }
        for (int i = 0; i < k; i++) {
            a_k[i] /= sqrtf(a2);
        }
    }
    //?????
    float a2 = 0;
    for (int i = 0; i < k; i++) {
        a2 += (a_k[i] * a_k[i]);
    }
    //norm
    //ak
    for (int i = 0; i < k; i++) {
        a_k[i] /= sqrtf(a2);
    }

    // vector<float> x(N);
    // vector<float> y(N);

    for (int i = 0; i < N; i++) {
        nodes[i].x = 0;
        nodes[i].y = 0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            nodes[i].x  += C[i][j] * b_k[j];
            nodes[i].y += C[i][j] * a_k[j];
        }
    }
    /*for (int i = 0; i < N; i++) {
        cout<<nodes[i].x<<" "<<nodes[i].y<<" ";
    }*/


}