

#include<iostream>
#include<fstream>
#include<string>
#include<random>
#include<algorithm>
#include<math.h>
#include "stdlib.h"
#include "time.h"
#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>

using namespace lemon;
using namespace std;


void generate_random_tree(ListGraph& g, int n)
{
    ListGraph::Node* nodes = new ListGraph::Node[n];
    //ListGraph::Node u;
    for(int i=0; i<n; i++) nodes[i] = g.addNode();
    srand(time(NULL));
    for(int i=1; i<n; i++)
    {
        int j = (rand() % i);
        int k = (rand() % (n-i)) + i;
        g.addEdge(nodes[j],nodes[k]);
        ListGraph::Node temp = nodes[i];
        nodes[i] = nodes[k];
        nodes[k] = temp;
    }
    delete[] nodes;
    return;
}

void generate_random_connected_graph(ListGraph& g, int n, int m)
{

    generate_random_tree(g,n);
    srand(time(NULL));
    for(int i=0; i<m-n+1; i++)
    {
        if(countEdges(g) == n*(n-1)/2) break;
        int j = (rand() % n);
        while(countIncEdges(g,g.nodeFromId(j)) == n-1) j = (rand() % n);
        int k = (rand() % n);
        ListGraph::Edge e = findEdge(g,g.nodeFromId(j),g.nodeFromId(k));
        while(j==k || e != INVALID)
        {
            k = (rand() % n);
            e = findEdge(g,g.nodeFromId(j),g.nodeFromId(k));
        }
        g.addEdge(g.nodeFromId(j),g.nodeFromId(k));
    }
    return;
}

void generate_random_groups(vector<vector<int>>& groups, int n, int k)
{
    int nV0;
    float f;
    do
    {
        f = (float)(rand() % (RAND_MAX>>1)) / RAND_MAX;
        nV0 = (int) (f*n);
    } while( nV0==0 || n-nV0 <= k);
    int nVG = n - nV0 - k;
    int nVperGroup = (int)(nVG/k);

    vector<int> nodeInd;
    for(int i=0; i<n; i++) nodeInd.push_back(i);
    random_shuffle(nodeInd.begin(),nodeInd.end());
    for(int i=0; i<nV0; i++) nodeInd.pop_back();
    for(int i=0; i<k; i++){ groups[i].push_back(nodeInd.back()); nodeInd.pop_back();}
    srand(time(NULL));
    for(int i = 0; i<k; i++)
    {
        for(int j=0; j<nVperGroup; j++)
        {
            groups[i].push_back(nodeInd[rand() % nVG]);
        }
    }

    return;
}

double euclidean_dist(pair<int,int>& P, pair<int,int>& Q)
{
    int fr = P.first - Q.first;
    int sc = P.second - Q.second;
    return sqrt(fr*fr + sc*sc);
}

int manhattan_dist(pair<int,int>& P, pair<int,int>& Q)
{
    return abs(P.first - Q.first) + abs(P.second - Q.second);
}




void random_weights(int m, vector<int>& wght)
{
    for(int i=0; i<m; i++) wght.push_back(rand() % 100 +1);
    return;
}


void random_euclidean_dist(vector<pair<int,int>>& edges, vector<int>& wght, int n, int m)
{
    vector<pair<int,int>> points;
    srand(time(NULL));
    for(int i=0; i<n; i++) points.push_back(make_pair(rand()%100 + 1,rand()%100 + 1));
    for(vector<pair<int,int>>::iterator it = edges.begin(); it!=edges.end(); ++it) wght.push_back(ceil(euclidean_dist(points[it->first],points[it->second])));
    return;

}


void random_manhattan_dist(vector<pair<int,int>>& edges, vector<int>& wght, int n, int m)
{
    vector<pair<int,int>> points;
    srand(time(NULL));
    for(int i=0; i<n; i++) points.push_back(make_pair(rand()%100 + 1,rand()%100 + 1));
    for(vector<pair<int,int>>::iterator it = edges.begin(); it!=edges.end(); ++it) wght.push_back(manhattan_dist(points[it->first],points[it->second]));
    return;

}



void generate_random_gst_instance(ofstream& inFile, string instance_name, string remark, string type, int m, int n, int k)
{

    ListGraph g;
    vector<vector<int>> groups(k);
    generate_random_connected_graph(g,n,m);
    generate_random_groups(groups,n,k);
    vector<pair<int,int>> edges;
    for(ListGraph::EdgeIt e(g); e!=INVALID; ++e) edges.push_back(make_pair(min(g.id(g.u(e)),g.id(g.v(e))),max(g.id(g.u(e)),g.id(g.v(e)))));
    sort(edges.begin(),edges.end());
    vector<int> wght;

    if(type=="rand") random_weights(m,wght);
    if(type=="euclid") random_euclidean_dist(edges,wght,n,m);
    if(type=="grid") random_manhattan_dist(edges,wght,n,m);
    //for(ListGraph::EdgeIt e(g); e!=INVALID; ++e) cout<<g.id(g.u(e))<<"\t"<<g.id(g.v(e))<<endl;

    inFile<<"0 string 33D32945  STP Steiner Tree Problem File"<<endl;
    inFile<<"SECTION Comment"<<endl;
    inFile<<"Name '"<<instance_name<<"'"<<endl;
    inFile<<"Remark 'T.-D. Nguyen and P.-T. Do: An ant colony optimization algorithm for solving Group Steiner Proble/home/sjelic/Dokumenti/stxxlProjektim (2013)'"<<endl;
    inFile<<"END"<<endl<<endl;
    inFile<<"SECTION Graph"<<endl;
    inFile<<"Nodes "<<n<<endl;
    inFile<<"Edges "<<m<<endl;
    int i=0;
    for(vector<pair<int,int>>::iterator it = edges.begin(); it!=edges.end(); ++it) inFile<<"E "<<it->first<<" "<<it->second<<" "<<wght[i++]<<endl;
    inFile<<"END"<<endl<<endl;
    inFile<<"SECTION Terminals"<<endl;
    inFile<<"Terminals "<<k<<endl;
    for(int i=0; i<k; i++)
    {
        inFile<<"T";
        for(unsigned int j=0; j<groups[i].size(); j++) inFile<<" "<<groups[i][j];
        inFile<<endl;
    }
    inFile<<"END"<<endl<<endl;
    inFile<<"EOF";
    return;
}
