#include "layout.h"
using namespace std;

int main(int argc, char **argv)
{
    string filename = "visbrazil.mtx";

    graph g;
    g.initgraph("alldata/" + filename);
    g.initRandomPosition();

    int tmax = 50;
    int choose_metric = 4;

    AutoFDP(g, tmax, choose_metric);

    g.drawPathSVG(filename, "result/");

   return 0;
}