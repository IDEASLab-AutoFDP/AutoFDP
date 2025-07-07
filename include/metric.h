#ifndef _METRIC_H_
#define _METRIC_H_
#include <Node.h>
#include <Link.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <Edge.h>
#include <algorithm>
#include <Eigen/Dense>

#define PI 3.141592653

int judgeCross(Node x,Node y,Node z,Node h);
float multi(Node p1,Node p2,Node p0);
Node intersect(Node a,Node b,Node c,Node d);


double testCL(std::vector<Node> &nodes,std::vector<Link> &links);
double testMA(std::vector<Node> &nodes,std::vector<Link> &links);
float testSE(std::vector<Node> &nodes,std::vector<Link> &links,std::vector<std::vector<double>> &sp);
double testIL(std::vector<Node> &nodes, std::vector<Link> &links);
float testAR(std::vector<Node> &nodes);
double testCA(std::vector<Node> &nodes,std::vector<Link> &links);
double testGP(std::vector<Node> &nodes, std::vector<Link> &links);
double testNR(std::vector<Node> &nodes);
float testNP(std::vector<Node> &nodes,std::vector<Link> &links,std::vector<std::vector<double>> &sp,int neigh);


float testMA_DeepGD(std::vector<Node> &nodes, std::vector<Link> &links);
double testNR_DeepGD(std::vector<Node> &nodes);
double testtSNE_DeepGD(std::vector<Node> nodes, std::vector<std::vector<double>> P);

#endif