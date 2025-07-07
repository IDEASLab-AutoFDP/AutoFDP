#include "solver.h"

int str2int(std::string s)
{
    std::stringstream stream(s);
    int i;
    stream >> i;
    stream.clear();
    return i;
}

float str2float(std::string s)
{
    std::stringstream stream(s);
    float i;
    stream >> i;
    stream.clear();
    return i;
}

void graph::get_degree()
{
    int V = nodes.size();
    std::vector<int> deg(V, 0);
    for (int i = 0; i < links.size(); i++)
    {
        deg[links[i].source]++;
        deg[links[i].target]++;
    }
    std::cout << "deg=[";
    for (int i = 0; i < V; i++)
    {
        std::cout << deg[i] << ",";
        if (i == V / 2)
        {
            std::cout << std::endl;
        }
    }
    std::cout << "]" << std::endl;
}

std::vector<int> split(const std::string &str, const std::string &delim)
{
    std::vector<int> res;
    if ("" == str)
        return res;
    char *strs = new char[str.length() + 1];
    strcpy(strs, str.c_str());

    char *d = new char[delim.length() + 1];
    strcpy(d, delim.c_str());

    char *p = strtok(strs, d);
    while (p)
    {
        std::string s = p;
        int y = atoi(s.c_str());
        res.push_back(y);
        p = strtok(NULL, d);
    }
    return res;
}

std::vector<std::string> split_string(const std::string &str, const std::string &delim)
{
    std::vector<std::string> res;
    if ("" == str)
        return res;
    char *strs = new char[str.length() + 1];
    strcpy(strs, str.c_str());

    char *d = new char[delim.length() + 1];
    strcpy(d, delim.c_str());

    char *p = strtok(strs, d);
    while (p)
    {
        std::string s = p;
        // int y = atoi(s.c_str());
        res.push_back(s);
        p = strtok(NULL, d);
    }
    return res;
}

graph::graph()
{
    n = 0;
    m = 0;

    a1 = 1;
    b1 = -2;
    a2 = 0;
    b2 = -1.5;
    k1 = 1.0;
    k2 = 1.0;
}

void graph::initgraph(std::string filename)
{
    std::ifstream input;
    input.open(filename.data());
    assert(input.is_open());
    std::string s;
    std::string b = "break";
    while (getline(input, s))
    {
        if (s.at(0) != '%')
        {
            break;
        }
    }
    std::vector<int> I, J;

    std::vector<int> di = split(s, " ");
    // if the first row is num of node and edge
    if (di.size() == 3)
    {
        n = di[1];
        m = di[2];
    }
    else
    {
        I.push_back(di[0]);
        J.push_back(di[1]);
    }

    bool cluster_file = false;
    while (getline(input, s))
    {
        if (s.at(0) == 'b')
        {
            cluster_file = true;
            break;
        }
        di = split(s, " ");
        I.push_back(di[0]);
        J.push_back(di[1]);
    }
    if (cluster_file)
    {
        while (getline(input, s))
        {
            flag.push_back(str2int(s));
        }
    }
    input.close();

    int max_i = *max_element(I.begin(), I.end());
    int max_j = *max_element(J.begin(), J.end());
    int min_i = *min_element(I.begin(), I.end());
    int min_j = *min_element(J.begin(), J.end());
    int max_ind = max_i > max_j ? max_i : max_j;
    int min_ind = min_i < min_j ? min_i : min_j;
    n = max_ind - min_ind + 1;
    m = I.size();
    std::vector<std::set<int>> undirected(n);

    for (int ij = 0; ij < m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw std::invalid_argument("i or j bigger than n");

        if (i != j && undirected[j].find(i) == undirected[j].end()) // if edge not seen
        {
            undirected[i].insert(j);
            undirected[j].insert(i);
            Link temp;
            temp.source = i - min_ind;
            temp.target = j - min_ind;
            links.push_back(temp);
        }
    }

    m = links.size();

    nodes.resize(n);

    std::cout << "build undirected graph n=" << n << " m=" << m << std::endl;
}

void graph::initRandomPosition()
{
    if (nodes.size() != n)
    {
        nodes.resize(n);
    }
    for (int i = 0; i < this->n; i++)
    {
        double radius = 10.0 * sqrt(i), angle = i * 3.141592653589793f * (3.0 - sqrt(5.0));
        Node temp;
        temp.id = i;
        temp.x = radius * cos(angle);
        temp.y = radius * sin(angle);
        nodes[i] = temp;
    }
    for (int i = 0; i < this->m; i++)
    {
        nodes[links[i].source].deg++;
        nodes[links[i].target].deg++;
    }
}

void graph::solveDij()
{

    std::vector<Edge> e[n];
    int st, end, weight;
    for (int i = 0; i < this->m; i++)
    {
        st = links[i].source;
        end = links[i].target;
        Edge tmp;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
        st = links[i].target;
        end = links[i].source;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
    }
    for (int i = 0; i < n; i++)
    {
        std::vector<float> dist(n);
        Dij(n, m, e, i, dist);
        std::vector<double> isp(n);
        for (int j = 0; j < n; j++)
        {
            if (dist[j] == MAX)
            {
                isp[j] = 0;
            }
            else
            {
                isp[j] = dist[j];
            }
        }
        sp.push_back(isp);
    }
}

void graph::preSolveSGD(double eps, int t_max, int seed, float eta_max)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (i != j)
                constraints_sgd.push_back({i, j});
        }
    }

    rk_seed(seed, &rstate);
    fisheryates_shuffle(constraints_sgd, rstate);


    float eta_min = eps;
    float lambd = log(eta_min / eta_max) / (t_max - 1);
    for (int i = 0; i < t_max; i++)
    {
        schedule.push_back(eta_max * exp(lambd * i));
    }
}

bool graph::solve_bilevel(int iter)
{
    // SGD
    double Delta = 0.001;
    double delta_avg = 0;
    double sum_delta = 0;
    float alpha = 0;
    if (iter % 10 == 1)
        fisheryates_shuffle(constraints_sgd, rstate);

    double Delta_max = 0;

    alpha = schedule[iter];
    for (int con = 0; con < constraints_sgd.size(); con++)
    {
        float rx = 0;
        float ry = 0;

        const std::pair<int,int> &t = constraints_sgd[con];
        const int &i = t.first, &j = t.second;
        const double d_ij = sp[i][j];
        alpha = schedule[iter];

        if (d_ij == 0)
            continue;

        float mvx = (nodes[i].x - nodes[j].x);
        float mvy = (nodes[i].y - nodes[j].y);
        float dist = sqrtf(mvx * mvx + mvy * mvy + 0.001f) + 0.001f;
        if (mvx != mvx)
        {
            printf("mvx error\n");
            std::cout << nodes[i].x << " " << nodes[i].y << " " << nodes[j].x << " " << nodes[j].y << std::endl;
            return 1;
        }
        if (mvy != mvy)
        {
            printf("mvy error\n");
            return 1;
        }
        if (dist < 0.01f)
            dist = 0.01f;
        float fccur = 1 / dist;
        float wcur = 1 / d_ij;
        float f = 0;
        bool isSping = true;
        float ka = k1, kr = k2, fa = 0, fr = 0;
        if (ka != 0)
        {
            fa = ka;
            fa *= pow(dist, a1);
            fa *= pow(d_ij, b1);
        }
        // a
        fr = 0;
        if (kr != 0)
        {
            fr = kr;

            fr *= pow(dist, a2);
            fr *= pow(d_ij, b2);
        }
        // r
        f = alpha * (fa - fr);
        if (f != f)
        {
            fprintf(stderr, "solve fffff error\n");
            fprintf(stderr, "ka=%f kr=%f a1=%f a2=%f b1=%f b2=%f\n", ka, kr, a1, a2, b1, b2);
            f = 1e4;
        }

        if (abs(f) > abs(dist - d_ij))
        {
            f = dist - d_ij;
        }
        if (abs(f) > 1e3 * alpha)
        {
            if (f > 0)
            {
                f = 1e3 * alpha;
            }
            else
            {
                f = -1.0 * 1e3 * alpha;
            }
        }
        rx += 1.0f * mvx / dist * f * 0.5;
        ry += 1.0f * mvy / dist * f * 0.5;

        nodes[i].x -= rx;
        nodes[i].y -= ry;
        nodes[j].x += rx;
        nodes[j].y += ry;
    }

    return 1;
}

std::vector<Node> graph::solve_k_with_noupdate(int iter, float Delta_a1, float Delta_b1, float Delta_a2, float Delta_b2, float Delta_k1, float Delta_k2)
{
    // SGD
    float alpha = 0;
    std::vector<Node> _nodes(nodes);
    float _a1 = a1, _a2 = a2, _b1 = b1, _b2 = b2, _k1 = k1, _k2 = k2;
    if (Delta_a1)
    {
        _a1 = a1 + Delta_a1;
    }
    if (Delta_a2)
    {
        _a2 = a2 + Delta_a2;
    }
    if (Delta_b1)
    {
        _b1 = b1 + Delta_b1;
    }
    if (Delta_b2)
    {
        _b2 = b2 + Delta_b2;
    }
    if (Delta_k1)
    {
        _k1 = k1 + Delta_k1;
    }
    if (Delta_k2)
    {
        _k2 = k2 + Delta_k2;
    }

    for (int con = 0; con < constraints_sgd.size(); con++)
    {
        float rx = 0;
        float ry = 0;

        const std::pair<int,int> &t = constraints_sgd[con];
        const int &i = t.first, &j = t.second;
        const double d_ij = sp[i][j];
        alpha = schedule[iter];

        if (d_ij == 0)
            continue;

        float mvx = (_nodes[i].x - _nodes[j].x);
        float mvy = (_nodes[i].y - _nodes[j].y);
        float dist = sqrtf(mvx * mvx + mvy * mvy + 0.001f) + 0.001f;
        if (mvx != mvx)
        {
            printf("update mvx error\n");
            std::cout << "mvx=" << mvx << " mvy=" << mvy << std::endl;
            return _nodes;
        }
        if (mvy != mvy)
        {
            printf("update mvy error\n");
            return _nodes;
        }
        if (dist < 0.01f)
            dist = 0.01f;
        float fccur = 1 / dist;
        float wcur = 1 / d_ij;
        float f = 0;
        bool isSping = true;
        float ka = _k1, kr = _k2, fa = 0, fr = 0;
        if (ka != 0)
        {
            fa = ka;
            fa *= pow(dist, _a1);
            fa *= pow(d_ij, _b1);
        }
        // a
        if (kr != 0)
        {
            fr = kr;
            fr *= pow(dist, _a2);
            fr *= pow(d_ij, _b2);
        }
        // r
        f = alpha * (fa - fr);
        if (f != f)
        {
            fprintf(stderr, "notupdate fffff error\n");
            fprintf(stderr, "ka=%f kr=%f a1=%f a2=%f b1=%f b2=%f i=%d j=%d iter=%d\n", ka, kr, _a1, _a2, _b1, _b2, i, j, iter);
            fprintf(stderr, "fa=%f fr=%f f=%f dist=%f mvx=%f mvy=%f _nodes[i].x=%f _nodes[i].y=%f _nodes[j].x=%f _nodes[j].y=%f\n", fa, fr, f, dist, _nodes[i].x, _nodes[i].y, _nodes[j].x, _nodes[j].x);
        }

        if (abs(f) > abs(dist - d_ij))
        {
            f = dist - d_ij;
        }
        if (abs(f) > 1e3 * alpha)
        {
            if (f > 0)
            {
                f = 1e3 * alpha;
            }
            else
            {
                f = -1.0 * 1e3 * alpha;
            }
        }

        rx += 1.0f * mvx / dist * f * 0.5;
        ry += 1.0f * mvy / dist * f * 0.5;

        _nodes[i].x -= rx;
        _nodes[i].y -= ry;
        _nodes[j].x += rx;
        _nodes[j].y += ry;
    }
    return std::move(_nodes);
}

float graph::opti_NP(int iter, double Delta, int neigh, double k)
{
    float NP = testNP(nodes, links, sp, neigh);
    if (NP == 1.0)
    {
        return false;
    }
    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float NP_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            NP_d = testNP(_nodes, links, sp, neigh);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            NP_d = testNP(_nodes, links, sp, neigh);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            NP_d = testNP(_nodes, links, sp, neigh);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            NP_d = testNP(_nodes, links, sp, neigh);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            NP_d = testNP(_nodes, links, sp, neigh);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            NP_d = testNP(_nodes, links, sp, neigh);
        }

        float delta_d = (NP_d - NP) / Delta;
        delta[i] = delta_d;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 += (k * delta_a1);
    a2 += (k * delta_a2);
    b1 += (k * delta_b1);
    b2 += (k * delta_b2);
    k1 += (k * delta_k1);
    k2 += (k * delta_k2);
    return NP;
}

float graph::opti_IL_mult(int iter, double Delta, double k)
{
    float IL = testIL(nodes, links);
    if (iter == 0)
    {
        loss_scale = 0.5 / (IL + 0.01f);
    }

    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float IL_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            IL_d = testIL(_nodes, links);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            IL_d = testIL(_nodes, links);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            IL_d = testIL(_nodes, links);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            IL_d = testIL(_nodes, links);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            IL_d = testIL(_nodes, links);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            IL_d = testIL(_nodes, links);
        }

        float delta_d = (IL_d - IL) / Delta;
        delta[i] = delta_d * loss_scale;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 -= std::max(-0.1, std::min(0.1, k * delta_a1));
    a2 -= std::max(-0.1, std::min(0.1, k * delta_a2));
    b1 -= std::max(-0.1, std::min(0.1, k * delta_b1));
    b2 -= std::max(-0.1, std::min(0.1, k * delta_b2));
    k1 -= std::max(-0.1, std::min(0.1, k * delta_k1));
    k2 -= std::max(-0.1, std::min(0.1, k * delta_k2));
    return IL;
}

float graph::opti_IL(int iter, double Delta, double k)
{
    float IL = testIL(nodes, links);
    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float IL_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            IL_d = testIL(_nodes, links);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            IL_d = testIL(_nodes, links);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            IL_d = testIL(_nodes, links);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            IL_d = testIL(_nodes, links);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            IL_d = testIL(_nodes, links);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            IL_d = testIL(_nodes, links);
        }

        float delta_d = (IL_d - IL) / Delta;
        delta[i] = delta_d;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    std::cout<<"IL="<<IL<<std::endl;
    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 -= std::max(-0.2, std::min(0.2, k * delta_a1));
    a2 -= std::max(-0.2, std::min(0.2, k * delta_a2));
    b1 -= std::max(-0.2, std::min(0.2, k * delta_b1));
    b2 -= std::max(-0.2, std::min(0.2, k * delta_b2));
    k1 -= std::max(-0.2, std::min(0.2, k * delta_k1));
    k2 -= std::max(-0.2, std::min(0.2, k * delta_k2));
    return IL;
}

float graph::opti_SE(int iter, double Delta, double k)
{

    float SE = testSE(nodes, links, sp);

    float delta[6];
#pragma omp parallel for num_threads(4)
    for (int i = 0; i < 6; i++)
    {

        float SE_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            SE_d = testSE(_nodes, links, sp);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            SE_d = testSE(_nodes, links, sp);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            SE_d = testSE(_nodes, links, sp);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            SE_d = testSE(_nodes, links, sp);
        }
        else if(i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            SE_d = testSE(_nodes, links, sp);
        }
        else if(i== 5 )
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            SE_d = testSE(_nodes, links, sp);
        }

        float delta_d = (SE_d - SE) / Delta;
        delta[i] = delta_d;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;
    a1 -= k * delta_a1;
    a2 -= k * delta_a2;
    b1 -= k * delta_b1;
    b2 -= k * delta_b2;
    k1 -= k * delta_k1;
    k2 -= k * delta_k2;

    return SE;
}

float graph::opti_AR(int iter, double Delta, double k)
{
    float AR = testAR(nodes);
    if (iter == 0)
    {
        loss_scale = 0.1 / (AR + 0.01f);
    }

    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float AR_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            AR_d = testAR(_nodes);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            AR_d = testAR(_nodes);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            AR_d = testAR(_nodes);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            AR_d = testAR(_nodes);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            AR_d = testAR(_nodes);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            AR_d = testAR(_nodes);
        }
        float delta_d = (AR_d - AR) / Delta;
        delta[i] = delta_d * loss_scale;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 -= k * delta_a1;
    a2 -= k * delta_a2;
    b1 -= k * delta_b1;
    b2 -= k * delta_b2;
    k1 -= k * delta_k1;
    k2 -= k * delta_k2;

    return AR;
}

float graph::opti_CA(int iter, double Delta, double k)
{
    float CA = testCA(nodes, links);
    if (iter == 0)
    {
        loss_scale = 0.5 / (CA + 0.01f);
    }
    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float CA_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            CA_d = testCA(_nodes, links);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            CA_d = testCA(_nodes, links);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            CA_d = testCA(_nodes, links);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            CA_d = testCA(_nodes, links);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            CA_d = testCA(_nodes, links);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            CA_d = testCA(_nodes, links);
        }

        float delta_d = (CA_d - CA) / Delta;
        delta[i] = delta_d * loss_scale;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 -= (k * delta_a1);
    a2 -= (k * delta_a2);
    b1 -= (k * delta_b1);
    b2 -= (k * delta_b2);
    k1 -= (k * delta_k1);
    k2 -= (k * delta_k2);
    return CA;
}

float graph::opti_GP(int iter, double Delta, double k)
{
    float GA = testGP(nodes, links);
    if (iter == 0)
    {
        loss_scale = 0.5 / (GA + 0.01f);
    }
    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float GA_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            GA_d = testGP(_nodes, links);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            GA_d = testGP(_nodes, links);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            GA_d = testGP(_nodes, links);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            GA_d = testGP(_nodes, links);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            GA_d = testGP(_nodes, links);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            GA_d = testGP(_nodes, links);
        }

        float delta_d = (GA_d - GA) / Delta;
        delta[i] = delta_d * loss_scale;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    delta_a1 *= (k * 1e-2);
    delta_b1 *= (k * 1e-2);
    delta_a2 *= (k * 1e-2);
    delta_b2 *= (k * 1e-2);
    delta_k1 *= (k * 1e-2);

    delta_k2 *= (k * 1e-2);

    delta_a1 = delta_a1 > 0.2 ? 0.2 : delta_a1;
    delta_a1 = delta_a1 < -0.2 ? -0.2 : delta_a1;
    delta_b1 = delta_b1 > 0.2 ? 0.2 : delta_b1;
    delta_b1 = delta_b1 < -0.2 ? -0.2 : delta_b1;
    delta_b2 = delta_b2 > 0.2 ? 0.2 : delta_b2;
    delta_b2 = delta_b2 < -0.2 ? -0.2 : delta_b2;
    delta_a2 = delta_a2 > 0.2 ? 0.2 : delta_a2;
    delta_a2 = delta_a2 < -0.2 ? -0.2 : delta_a2;
    delta_k1 = delta_k1 > 0.1 ? 0.1 : delta_k1;
    delta_k1 = delta_k1 < -0.1 ? -0.1 : delta_k1;
    delta_k2 = delta_k2 > 0.1 ? 0.1 : delta_k2;
    delta_k2 = delta_k2 < -0.1 ? -0.1 : delta_k2;

    a1 -= delta_a1;
    a2 -= delta_a2;
    b1 -= delta_b1;
    b2 -= delta_b2;
    k1 -= delta_k1;
    k2 -= delta_k2;

    return GA;
}

float graph::opti_NR(int iter, double Delta, double k)
{
    float NR = testNR(nodes);
    if (iter == 0)
    {
        loss_scale = 0.1 / (NR + 0.01f);
    }

    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float NR_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            NR_d = testNR(_nodes);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            NR_d = testNR(_nodes);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            NR_d = testNR(_nodes);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            NR_d = testNR(_nodes);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            NR_d = testNR(_nodes);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            NR_d = testNR(_nodes);
        }

        float delta_d = (NR_d - NR) / Delta;
        delta[i] = delta_d * loss_scale;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 -= (k * delta_a1);
    a2 -= (k * delta_a2);
    b1 -= (k * delta_b1);
    b2 -= (k * delta_b2);
    k1 -= (k * delta_k1);
    k2 -= (k * delta_k2);

    return NR;
}

float graph::opti_NR_DeepGD(int iter, double Delta, double k)
{
    float NR = testNR_DeepGD(nodes);
    if (iter == 0)
    {
        loss_scale = 0.1 / (NR + 0.01f);
    }

    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float NR_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            NR_d = testNR_DeepGD(_nodes);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            NR_d = testNR_DeepGD(_nodes);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            NR_d = testNR_DeepGD(_nodes);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            NR_d = testNR_DeepGD(_nodes);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            NR_d = testNR_DeepGD(_nodes);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            NR_d = testNR_DeepGD(_nodes);
        }

        float delta_d = (NR_d - NR) / Delta;
        delta[i] = delta_d * loss_scale;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];


    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 -= (k * delta_a1);
    a2 -= (k * delta_a2);
    b1 -= (k * delta_b1);
    b2 -= (k * delta_b2);
    k1 -= (k * delta_k1);
    k2 -= (k * delta_k2);

    return NR;
}

float graph::opti_MA(int iter, double Delta, double k)
{
    float MA = testMA(nodes, links);
    if (iter == 0)
    {
        loss_scale = 0.5 / (MA + 0.01f);
    }
    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float MA_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            MA_d = testMA(_nodes, links);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            MA_d = testMA(_nodes, links);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            MA_d = testMA(_nodes, links);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            MA_d = testMA(_nodes, links);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            MA_d = testMA(_nodes, links);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            MA_d = testMA(_nodes, links);
        }

        float delta_d = (MA_d - MA) / Delta;
        delta[i] = delta_d * loss_scale;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    std::cout << "MA=" << MA << std::endl;
    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 += std::max(-0.2, std::min(0.2, k * delta_a1));
    a2 += std::max(-0.2, std::min(0.2, k * delta_a2));
    b1 += std::max(-0.2, std::min(0.2, k * delta_b1));
    b2 += std::max(-0.2, std::min(0.2, k * delta_b2));
    k1 += std::max(-0.2, std::min(0.2, k * delta_k1));
    k2 += std::max(-0.2, std::min(0.2, k * delta_k2));

    return MA;
}

float graph::opti_MA_DeepGD(int iter, double Delta, double k)
{
    float MA = testMA_DeepGD(nodes, links);
    if (iter == 0)
    {
        loss_scale = 0.5 / (MA + 0.01f);
    }
    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float MA_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            MA_d = testMA_DeepGD(_nodes, links);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            MA_d = testMA_DeepGD(_nodes, links);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            MA_d = testMA_DeepGD(_nodes, links);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            MA_d = testMA_DeepGD(_nodes, links);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            MA_d = testMA_DeepGD(_nodes, links);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            MA_d = testMA_DeepGD(_nodes, links);
        }

        float delta_d = (MA_d - MA) / Delta;
        delta[i] = delta_d * loss_scale;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    std::cout << "MA=" << MA << std::endl;
    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 -= std::max(-0.2, std::min(0.2, k * delta_a1));
    a2 -= std::max(-0.2, std::min(0.2, k * delta_a2));
    b1 -= std::max(-0.2, std::min(0.2, k * delta_b1));
    b2 -= std::max(-0.2, std::min(0.2, k * delta_b2));
    k1 -= std::max(-0.2, std::min(0.2, k * delta_k1));
    k2 -= std::max(-0.2, std::min(0.2, k * delta_k2));


    return MA;
}

float graph::opti_CL(int iter, double Delta, double k)
{
    float CR = testCL(nodes, links);
    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {

        float CR_d=0;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            CR_d = testCL(_nodes, links);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            CR_d = testCL(_nodes, links);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            CR_d = testCL(_nodes, links);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            CR_d = testCL(_nodes, links);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            CR_d = testCL(_nodes, links);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            CR_d = testCL(_nodes, links);
        }

        float delta_d = (CR_d - CR) / Delta;
        delta[i] = delta_d;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 += (k * delta_a1);
    a2 += (k * delta_a2);
    b1 += (k * delta_b1);
    b2 += (k * delta_b2);
    k1 += (k * delta_k1);
    k2 += (k * delta_k2);
    return CR;
}


void graph::drawPathSVG(std::string filename, std::string path, bool outpos)
{
    std::string filename1 = path + filename.substr(0, filename.length() - 4) + "_layout.svg";
    std::string newpath = path;
    // foder not exist
    if (IsFolderExist(newpath.c_str()) != -1)
    {
        int isCreate = light::mkpath(newpath);
        std::cout << "folder create" << isCreate << std::endl;
    }

    float maxx = -1e9, minx = 1e9, maxy = -1e9, miny = 1e9;
    std::vector<Node> _nodes = nodes;

    for (int i = 0; i < _nodes.size(); i++)
    {
        if (_nodes[i].x > maxx)
            maxx = _nodes[i].x;
        if (_nodes[i].x < minx)
            minx = _nodes[i].x;
        if (_nodes[i].y > maxy)
            maxy = _nodes[i].y;
        if (_nodes[i].y < miny)
            miny = _nodes[i].y;
    }
    float lenth = std::max(maxx - minx, maxy - miny);
    for (int i = 0; i < _nodes.size(); i++)
    {
        _nodes[i].x -= minx;
        _nodes[i].y -= miny;
    }
    for (int i = 0; i < _nodes.size(); i++)
    {
        _nodes[i].x /= lenth;
        _nodes[i].x *= 1000;
        _nodes[i].x += 50;
        _nodes[i].y /= lenth;
        // _nodes[i].y *= 1000 + 100;
        _nodes[i].y *= 1000;
        _nodes[i].y += 50;
    }

    fprintf(stderr, "new drawSVG\n");
    std::ofstream f;
    f.open(filename1, std::ios::out);
    // FILE *fp1=freopen((char*) filename1.data(),"w",stdout);
    f << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n"
         "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
         "\n"
         "<svg width=\"1100\" height=\"1100\" version=\"1.1\"\n"
         "xmlns=\"http://www.w3.org/2000/svg\">"
      << std::endl;

    for (int i = 0; i < this->m; i++)
    {
        f << "<line x1=\"" << _nodes[links[i].source].x << "\" y1=\"" << _nodes[links[i].source].y << "\" x2=\"" << _nodes[links[i].target].x << "\" y2=\"" << _nodes[links[i].target].y << "\"\n"                                                                                                                                                                      "style=\"stroke:rgb(99,99,99);stroke-width:4\"/>"
          << std::endl;
    }

    for (int i = 0; i < this->n; i++)
    {
        f << "<circle id=\"" << i <<"\" cx=\"" << _nodes[i].x << "\" cy=\"" << _nodes[i].y << "\" r=\"6\" stroke=\"RGB(31,120,180)\"\n"
                                                                            "stroke-width=\"4\" fill=\"RGB(31,120,180)\"/>"
          << std::endl;
    }

    f << "</svg>" << std::endl;
    f.close();
}

float graph::opti_tSNE_DeepGD(int iter, double Delta, double k)
{
    float tsne = testtSNE_DeepGD(nodes,P);
    float delta[6];
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < 6; i++)
    {
        float tsne_d;
        if (i == 0)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, Delta, 0, 0, 0, 0, 0);
            tsne_d = testtSNE_DeepGD(_nodes,P);
        }
        else if (i == 1)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, Delta, 0, 0, 0, 0);
            tsne_d = testtSNE_DeepGD(_nodes,P);
        }
        else if (i == 2)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, Delta, 0, 0, 0);
            tsne_d = testtSNE_DeepGD(_nodes,P);
        }
        else if (i == 3)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, Delta, 0, 0);
            tsne_d = testtSNE_DeepGD(_nodes,P);
        }
        else if (i == 4)
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, Delta, 0);
            tsne_d = testtSNE_DeepGD(_nodes,P);
        }
        else
        {
            std::vector<Node> _nodes = solve_k_with_noupdate(iter, 0, 0, 0, 0, 0, Delta);
            tsne_d = testtSNE_DeepGD(_nodes,P);
        }

        float delta_d = (tsne_d - tsne) / Delta;
        delta[i] = delta_d;
    }

    float delta_a1 = delta[0];
    float delta_b1 = delta[1];
    float delta_a2 = delta[2];
    float delta_b2 = delta[3];
    float delta_k1 = delta[4];
    float delta_k2 = delta[5];

    std::cout<<"tsne= "<<tsne<<std::endl;
    
    delta_a1 *= 1e-2;
    delta_b1 *= 1e-2;
    delta_a2 *= 1e-2;
    delta_b2 *= 1e-2;
    delta_k1 *= 1e-2;
    delta_k2 *= 1e-2;

    a1 -= (k * delta_a1);
    a2 -= (k * delta_a2);
    b1 -= (k * delta_b1);
    b2 -= (k * delta_b2);
    k1 -= (k * delta_k1);
    k2 -= (k * delta_k2);
    return tsne;
}

void graph::calP_DeepGD()
{
    P.clear();
    for (int i = 0; i < n; i++)
    {
        std::vector<double> tmp(n, 0);
        for (int j = 0; j < n; j++)
        {
            if(j==i)
                continue;
            tmp[j] = exp(-2 * sp[i][j] / 2);
        }
        
        double _sum = std::accumulate(tmp.begin(), tmp.end(), 0.0);
        for (int j = 0; j < n; j++)
            tmp[j] /= _sum;
        P.push_back(tmp);
    }
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
        {
            double v1=P[i][j];
            double v2=P[j][i];
            P[i][j]=(v1+v2)/(2*n);
            P[j][i]=(v1+v2)/(2*n);
        }
}