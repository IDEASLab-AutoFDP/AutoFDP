#include "metric.h"

#define PAI acos(-1)

int judgeCross(Node x, Node y, Node z, Node h)
{
    float m = (x.x - z.x) * (h.y - z.y) - (x.y - z.y) * (h.x - z.x);
    float n = (y.x - z.x) * (h.y - z.y) - (y.y - z.y) * (h.x - z.x);
    if (m * n <= 0)
        return 1;
    return 0;
}
float multi(Node p1, Node p2, Node p0)
{
    return (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
}
Node intersect(Node a, Node b, Node c, Node d)
{
    Node p;
    p.x = (multi(a, d, c) * b.x - multi(b, d, c) * a.x) / (multi(a, d, c) - multi(b, d, c));
    p.y = (multi(a, d, c) * b.y - multi(b, d, c) * a.y) / (multi(a, d, c) - multi(b, d, c));
    return p;
}
double testCL(std::vector<Node> &nodes, std::vector<Link> &links)
{
    int V = nodes.size(), E = links.size();
    float n_all = E * (E - 1) / 2;
    float n_imp = 0, c_max = 0, cn = 0;
    float crosslessness = 1;
    // compute degree
    std::vector<int> count(V, 0);
    for (int i = 0; i < E; i++)
    {
        count[links[i].source]++;
        count[links[i].target]++;
    }

    for (int i = 0; i < V; i++)
    {
        n_imp += count[i] * (count[i] - 1);
    }
    n_imp /= 2;
    c_max = n_all - n_imp;
    if (c_max == 0)
        std::cout << "cross error" << n_all << " " << E << std::endl;
    for (int i = 0; i < E; i++)
    {
        int src = links[i].source, tgt = links[i].target;
        if (src == tgt)
            continue;

        float src_x = nodes[src].x, src_y = nodes[src].y;
        float tgt_x = nodes[tgt].x, tgt_y = nodes[tgt].y;
        for (int j = i + 1; j < E; j++)
        {
            int temp_src = links[j].source, temp_tgt = links[j].target;
            if (src == temp_src || tgt == temp_src || src == temp_tgt || tgt == temp_tgt)
                continue;
            if (temp_tgt == temp_src)
                continue;
            // compute if cross
            Node a, b, c, d;
            a.x = src_x, a.y = src_y;
            b.x = tgt_x, b.y = tgt_y;
            c.x = nodes[temp_src].x, c.y = nodes[temp_src].y;
            d.x = nodes[temp_tgt].x, d.y = nodes[temp_tgt].y;

            int ans1 = judgeCross(a, b, c, d);
            int ans2 = judgeCross(c, d, a, b);

            if (ans1 == 1 && ans2 == 1)
            {
                cn++;
            }
        }
    }
    crosslessness = 1 - sqrtf(cn / c_max);
    // std::cout<<"cross="<<crosslessness<<std::endl;
    return crosslessness;
}

double testCA(std::vector<Node> &nodes,std::vector<Link> &links){
    int V=nodes.size(),E=links.size();
    float n_all=E*(E-1)/2;
    float n_imp=0,c_max=0,cnt=0;
    float crosslessness=1;
    //compute degree
    std::vector<int> count(V,0);
    for(int i=0;i<E;i++){
        count[links[i].source]++;
        count[links[i].target]++;
    }

    for(int i=0;i<V;i++){
        n_imp+=count[i]*(count[i]-1);
    }
    n_imp/=2;
    c_max=n_all-n_imp;
    float sum_cr=0;
    if(c_max==0) std::cout<<"cross error"<<n_all<<" "<<E<<std::endl;
    for(int i=0;i<E;i++){
        int src=links[i].source,tgt=links[i].target;
        if(src==tgt) continue;

        float src_x=nodes[src].x,src_y=nodes[src].y;
        float tgt_x=nodes[tgt].x,tgt_y=nodes[tgt].y;
        for(int j=i+1;j<E;j++){
            int temp_src=links[j].source,temp_tgt=links[j].target;
            if(src==temp_src||tgt==temp_src||src==temp_tgt||tgt==temp_tgt) continue;
            if(temp_tgt==temp_src) continue;
            //compute if cross
            Node a,b,c,d;
            a.x=src_x,a.y=src_y;
            b.x=tgt_x,b.y=tgt_y;
            
           // fprintf(stderr,"xx1=%f yy1=%f xx2=%f yy2=%f \n",nodes[temp_src].x,nodes[temp_src].y,nodes[temp_tgt].x,nodes[temp_tgt].y);
            c.x=nodes[temp_src].x,c.y=nodes[temp_src].y;
            d.x=nodes[temp_tgt].x,d.y=nodes[temp_tgt].y;
            //fprintf(stderr,"c.x=%f c.y=%f d.x=%f d.y=%f \n",c.x,c.y,d.x,d.y);
            int ans1=judgeCross(a,b,c,d);
            int ans2=judgeCross(c,d,a,b);

            float ag=0;
            if(ans1==1&&ans2==1){
                //cross angle
               
                double vex1=tgt_x-src_x,vex2=nodes[temp_tgt].x-nodes[temp_src].x;
                double vey1=tgt_y-src_y,vey2=nodes[temp_tgt].y-nodes[temp_src].y;
                double mul1= (vex1*vex2+vey1*vey2);
                double norm1=sqrtf(vex1*vex1+vey1*vey1),norm2=sqrtf(vex2*vex2+vey2*vey2);
                if(norm1==0||norm2==0){continue;}
                double norm=mul1/(norm2*norm1);
                if(norm>1){norm=0.99;}
                if(norm<-1){norm=-0.99;}
                double cos_sim=acos(norm);
                if(cos_sim!=cos_sim){
                    fprintf(stderr,"x1=%f y1=%f x2=%f y2=%f  xx1=%f yy1=%f xx2=%f yy2=%f \n"
                    ,src_x,src_y,tgt_x,tgt_y,c.x,c.y,d.x,d.y);
                    fprintf(stderr,"x1=%f y1=%f x2=%f y2=%f  xx1=%f yy1=%f xx2=%f yy2=%f \n"
                    ,a.x,a.y,b.x,b.y,nodes[temp_src].x,nodes[temp_src].y,nodes[temp_tgt].x,nodes[temp_tgt].y);
                     fprintf(stderr,"vex1=%f vex2=%f vey1=%f vey2=%f \n",vex1,vex2,vey1,vey2);
                    fprintf(stderr,"mult=%f norm1=%f norm2=%f cos_sim=%f\n",mul1,norm1,norm2,cos_sim);
                    fprintf(stderr,"cos_sim error\n");
                    cos_sim=PI;
                    continue;
                }
                cnt++;
                sum_cr+=abs((cos_sim-PI/2)/(PI/2));
            }

        }
    }
    if(cnt==0){return 0;}
    float cross_angle=sum_cr;//cnt;
    return cross_angle;
}

float testSE(std::vector<Node> &nodes, std::vector<Link> &links, std::vector<std::vector<double>> &sp)
{
    double edgelen = 0;
    int edgenum = 0;
    float sum_ddij = 0, sum_d = 0, scale = 1;
    for (int i = 0; i < nodes.size(); i++)
    {
        for (int j = 0; j < nodes.size(); j++)
        {
            if (sp[i][j] == 0)
                continue;
            float disx = nodes[i].x - nodes[j].x, disy = nodes[i].y - nodes[j].y;
            sum_ddij += (sqrt(disx * disx + disy * disy) * sp[i][j]) / (sp[i][j] * sp[i][j]);
            sum_d += (disx * disx + disy * disy) / (sp[i][j] * sp[i][j]);
        }
    }
    scale = sum_ddij / sum_d;
    for (int i = 0; i < nodes.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (sp[i][j] == 0)
                continue;
            double d = sqrt(pow(scale * (nodes[i].x - nodes[j].x), 2) + pow(scale * (nodes[i].y - nodes[j].y), 2));
            double dis = pow(d - sp[i][j], 2);
            edgelen += dis / pow(sp[i][j], 2);
            edgenum++;
        }
    }
    return edgelen / edgenum;
}

double testIL(std::vector<Node> &nodes, std::vector<Link> &links)
{
    double edgelen = 0;
    int E = links.size();
    double sum_edgelen = 0, scale = 1;

    for (int e = 0; e < links.size(); e++)
    {
        int i = links[e].source, j = links[e].target;
        double disx = nodes[i].x - nodes[j].x, disy = nodes[i].y - nodes[j].y;
        sum_edgelen += (sqrt(disx * disx + disy * disy));
    }
    scale = sum_edgelen / (float)E;

    for (int e = 0; e < links.size(); e++)
    {
        int i = links[e].source, j = links[e].target;
        double d = sqrtf(pow((nodes[i].x - nodes[j].x), 2) + pow((nodes[i].y - nodes[j].y), 2)) / scale;
        double dis = pow(d - 1, 2);
        edgelen += dis;
    }
    return edgelen / (double)E;
}

float testNP(std::vector<Node> &nodes,std::vector<Link> &links,std::vector<std::vector<double>> &sp,int neigh)
{
    double NP = 0;
    int leafnum = 0;
    int number=nodes.size();
    int countt=0;
    for (int i = 0; i < number; i++) {
        std::vector<int> dijknn;
        std::vector<double> dknn;
        for (int j = 0; j < number; j++) {
            if (j != i) {
                if (sp[i][j] <= neigh&&sp[i][j]>0) {
                    dijknn.push_back(j);
                }
            }
        }
        for (int j = 0; j < dijknn.size(); j++) {
            dknn.push_back(1000000000000);
        }
        for (int j = 0; j < number; j++) {
            if (j != i) {
                double curdis = sqrt(pow(nodes[i].x - nodes[j].x, 2) + pow(nodes[i].y - nodes[j].y, 2));
                for (int k = 0; k < dijknn.size(); k++) {
                    if (curdis < dknn[k]) {
                        for (int l = dijknn.size() - 1; l > k; l--) {
                            dknn[l] = dknn[l - 1];
                        }
                        dknn[k] = curdis;
                        break;
                    }
                }
            }
        }
        int intersect = 0;
        int uni = 0;
        for (int j = 0; j < dijknn.size(); j++) {
            double curdis = sqrt(pow(nodes[i].x - nodes[dijknn[j]].x, 2) + pow(nodes[i].y - nodes[dijknn[j]].y, 2));
            if (curdis <= dknn[dknn.size() - 1]) {
                intersect++;
            }
            else {
                uni += 2;
            }
        }
        if (dijknn.size() == 0) { NP += 0; leafnum += 1; }
        else { NP += (double)intersect / (uni + intersect); }
        countt+=intersect;
    }
    //cout << "NP" <<neigh<<"="<<NP / number << endl;
    return NP / number;
}

double testMA(std::vector<Node> &nodes, std::vector<Link> &links)
{
    int V = nodes.size(), E = links.size();
    double min_angle = 0, sum_angle = 0;

    std::vector<Edge> e[V];
    int st, end, weight;
    for (int i = 0; i < E; i++)
    {
        // cin>>st>>end>>weight;
        st = links[i].source;
        end = links[i].target;
        Edge tmp;
        // tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp); // vector
        st = links[i].target;
        end = links[i].source;
        // tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
    }
    // cout<<"111"<<endl;
    for (int i = 0; i < V; i++)
    {
        if (e[i].size() < 2)
            continue;
        double theta = (double)360 / (double)e[i].size();
        double min_v_angle = 360;
        // cout<<count[i]<<" "<<e[i].size()<<" ";
        for (int j = 0; j < e[i].size() - 1; j++)
        {
            int v0 = i;
            int v1 = e[i][j].end;
            double v1_x = nodes[v0].x - nodes[v1].x;
            double v1_y = nodes[v0].y - nodes[v1].y;
            for (int k = j + 1; k < e[i].size(); k++)
            {
                int v2 = e[i][k].end;
                // compute angel
                double v2_x = nodes[v0].x - nodes[v2].x;
                double v2_y = nodes[v0].y - nodes[v2].y;
                if ((sqrt(v1_x * v1_x + v1_y * v1_y) * sqrt(v2_x * v2_x + v2_y * v2_y)) == 0)
                    continue;
                double angle = ((v1_x * v2_x) + (v1_y * v2_y)) / (sqrt(v1_x * v1_x + v1_y * v1_y) * sqrt(v2_x * v2_x + v2_y * v2_y));
                if (angle <= 1 && angle >= -1)
                {
                    angle = acos(angle) * 180 / PI;
                }

                if (min_v_angle > angle)
                    min_v_angle = angle;
            }
        }
        double temp_angle = abs((theta - min_v_angle) / theta);
        sum_angle += temp_angle;
    }
    min_angle = 1 - sum_angle / (double)V;
    return min_angle;
}

float calTheta(float x1,float y1,float x2,float y2)
{
    float x = x2 - x1;
    float y = y2 - y1;
    float r = sqrt(x * x + y * y);
    // if(r == 0)
    //     return 0.0001;
    float theta = acos(x / r);
    if(y < 0)
        theta = 2 * PAI - theta;
    return theta;
}


// deepgd_loss
float testMA_DeepGD(std::vector<Node> &nodes, std::vector<Link> &links)
{
    float L_MA = 0;
    int n = nodes.size();
    std::vector<std::vector<int>> adjMtx(n);
    for(auto it:links)
    {
        int u = it.target;
        int v = it.source;
        adjMtx[u].push_back(v);
        adjMtx[v].push_back(u);
    }
    for (int i = 0; i < n; i++)
    {
        int degree = adjMtx[i].size();//该点度数
        if (degree == 1)
            continue;
        std::vector<float> thetas;
        for (auto it : adjMtx[i])
            thetas.push_back(calTheta(nodes[i].x, nodes[i].y, nodes[it].x, nodes[it].y));
        std::sort(thetas.begin(), thetas.end());
        
        float thetas_adjust = thetas[0];
        for (int j = 0; j < thetas.size(); j++)
            thetas[j] -= thetas_adjust;
        
        thetas.push_back(2 * PAI);
        std::vector<float> loss;
        for (int j = 1; j < thetas.size(); j++)
        {
            float angle = thetas[j] - thetas[j - 1];
            if(angle > PAI)
                angle = 2 * PAI - angle;
            loss.push_back(angle);
        }
        float optimal = 2 * PAI / degree;
        for (int j = 0; j < loss.size(); j++)
            L_MA += abs(optimal - loss[j]);
    }
    return L_MA;
}

float testAR(std::vector<Node> &nodes)
{
    int n = nodes.size();
    Eigen::MatrixXf pos(n, 2);
    for (int i = 0; i < n; i++)
    {
        pos(i, 0) = nodes[i].x;
        pos(i, 1) = nodes[i].y;
    }
    Eigen::VectorXf pos_mean = pos.colwise().mean();
    Eigen::MatrixX2f pos_centered = pos.rowwise() - pos_mean.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(pos_centered, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double ratio = svd.singularValues()[1] / svd.singularValues()[0];
    return -std::log(ratio);
}


double testNR(std::vector<Node> &nodes)
{
    int V = nodes.size();
    double r = 1 / sqrtf((double)V);

    double sum_vr = 0;
    double target = 0.1;
    double d_max = 0, d_min = 1e9;
    for (int i = 0; i < V; i++)
    {
        for (int j = i + 1; j < V; j++)
        {
            double mvx = nodes[i].x - nodes[j].x, mvy = nodes[i].y - nodes[j].y;
            double dist = sqrtf(mvx * mvx + mvy * mvy);
            if (dist > d_max)
            {
                d_max = dist;
            }
            if (d_min > dist)
            {
                d_min = dist;
            }
        }
    }
    d_min = d_max + 1;
    float cnt = 0;
    for (int i = 0; i < V; i++)
    {
        for (int j = i + 1; j < V; j++)
        {
            float mvx = nodes[i].x - nodes[j].x, mvy = nodes[i].y - nodes[j].y;
            float dist = sqrtf(mvx * mvx + mvy * mvy);
            float nr = (1 - dist / (r * d_max)) > 0 ? (1 - dist / (r * d_max)) : 0;
            if (nr != 0)
            {
                cnt++;
            }
            sum_vr += nr * nr;
        }
    }
    if (cnt == 0)
    {
        return 0;
    }
    cnt *= pow((1 - 1 / r), 2);
    return sum_vr;
}

double testNR_DeepGD(std::vector<Node> &nodes)
{
    int n = nodes.size();
    double NR = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double mvx = nodes[i].x - nodes[j].x, mvy = nodes[i].y - nodes[j].y;
            double dist = sqrt(mvx * mvx + mvy * mvy);
            NR += exp(-dist);
        }
    }
    NR *= 2;
    NR /= n * (n - 1);
    return NR;
}

double testGP(std::vector<Node> &nodes, std::vector<Link> &links)
{
    int V = nodes.size(), E = links.size();
    std::vector<Node> _nodes = nodes;
    float maxx = 0, minx = 100000000000, maxy = 0, miny = 10000000000;
    // fprintf(stderr,"_nodes0.x=%f _nodes1.x=%f nodes0.x=%f nodes1.x=%f \n",_nodes[0].x,_nodes[1].x,nodes[0].x,nodes[1].x);
    for (int i = 0; i < V; i++)
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
    float xlenth = maxx - minx;
    float ylenth = maxy - miny;
    float lenth = std::max(xlenth, ylenth);

    double min_ga = 1e12 + 1;
    double sum_ga = 0, sum_r = 0;
    for (int i = 0; i < E; i++)
    {
        int src = links[i].source, tgt = links[i].target;
        double x1 = nodes[src].x, y1 = nodes[src].y;
        double x2 = nodes[tgt].x, y2 = nodes[tgt].y;

        double c_x = (x1 + x2) / 2, c_y = (y1 + y2) / 2;
        double radii = 1e-6 + sqrtf((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / 2;

        for (int j = 0; j < V; j++)
        {
            if (j == src || j == tgt)
                continue;
            double k_x = nodes[j].x, k_y = nodes[j].y;
            double dist_to_min = sqrtf((k_x - c_x) * (k_x - c_x) + (k_y - c_y) * (k_y - c_y));
            double d_m = 0;
            if (radii - dist_to_min > 0)
            {
                d_m = radii - dist_to_min;
            }
            d_m /= lenth;
            sum_ga += d_m * d_m;
            sum_r += radii * radii;
            // if(min_ga>dist_to_min/radii){min_ga=dist_to_min/radii;}
            // sum_ga=
            // if(r-rad>0){
            //     sum_ga+=(r-rad)*(r-rad);
            // }
        }
    }
    return sum_ga; // sum_r;
}

std::vector<std::vector<double>> calQ(std::vector<Node> nodes)
{
    int n=nodes.size();
    std::vector<std::vector<double>> Q(n, std::vector<double>(n, 0.0));
    double sum=0.0;
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
        {
            double dx=nodes[i].x-nodes[j].x;
            double dy=nodes[i].y-nodes[j].y;
            double value = pow(1+dx*dx+dy*dy,-1);
            Q[i][j]=value;
            Q[j][i]=value;
            sum+=value;
        }
    sum*=2;
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
        {
            Q[i][j] /= sum;
            Q[j][i] /= sum;
        }
    return Q;
}

double testtSNE_DeepGD(std::vector<Node> nodes, std::vector<std::vector<double>> P)
{
    int n=nodes.size();
    std::vector<std::vector<double>> Q=calQ(nodes);
    double tsne=0.0;
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
            tsne += P[i][j]*log(P[i][j]/Q[i][j]);
    tsne*=2;
    return tsne;
}
