#include <iostream>
#include "gst_gen.h"

using namespace std;

int main()
{
    string type, name, remark;
    cin>>type;
    if(type=="rand") name = "RANDOM";
    if(type=="euclid") name = "EUCLIDEAN";
    if(type=="grid") name = "GRID";
    remark="T.-D. Nguyen and P.-T. Do: An ant colony optimization algorithm for solving Group Steiner Problem (2013)";
    int j=1;
    //for(int n=60; n<=300; n+=10)
    {

        int m,n,K;
        cin >> n >> m >> K;
        //int m = ( rand() % ((int)(n*(n-1)/4) - (int)(1.5*n)) ) + ((int)(1.5*n));
        //int K = ( rand() % ((int)(n/4) - (int)(log2(n))) ) + ((int)(log2(n)));
        stringstream ss;
        if(j<=9) ss<<"0";
        ss<<j<<"_"<<type<<"_"<<"n="<<n<<"_"<<"m="<<m<<"_"<<"k="<<K<<".stp";
        ofstream inFile(ss.str());
        generate_random_gst_instance(inFile, name, remark, type, m, n, K);
        j++;
    }

    return 0;
}
