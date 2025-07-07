
#ifndef _NODE_H_
#define _NODE_H_
#include <utility>
#include <string>
typedef std::pair<int,int> node_pair;
struct Node {
    int id;
    std::string _id;
    float x, y;
    int deg = 0;
};
#endif