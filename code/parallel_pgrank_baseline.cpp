#include<bits/stdc++.h>
#include<iostream>
#include<omp.h>
using namespace std;

vector <int> adj[1000000];
vector <int> rev[1000000];
int outdeg[1000000];
int N,M;
double d=0.15;
double threshold=0.0000000001;
double pr[1000000];
double prev1[1000000];
void compute(){
    double error=10000000000.0;
    int count=0;

    #pragma omp parallel for
    for(int i=1;i<=N;i++)
        prev1[i]=d/N;
    while(error>threshold){
        count+=1;
        #pragma omp parallel for  
        for(int i=1;i<=N;i++){
            pr[i]=d/(N*1.0);
            for(int v=0;v<rev[i].size();v++){
                int ver=rev[i][v];
                pr[i]+=((prev1[ver]*1.0)/(outdeg[ver]*1.0))*(1.0-d);
            }
        }
        double nerr=0.0;
        #pragma omp parallel for reduction(max : nerr)
        for(int i=1;i<=N;i++){
            nerr=max(nerr,abs(prev1[i]-pr[i]));
            prev1[i]=pr[i];
        }
        error=nerr;

    }
}

void print_pagerank(){
    cout<<"Final Page rank"<<endl;
    for(int i = 1; i <= N; i++){
        cout<<pr[i]<<endl;
    }
}

int main(){
    double time;
    time = omp_get_wtime();
    cin>>N>>M;
    for(int i=0;i<M;i++)
    {
        int u,v;
        cin>>u>>v;
        adj[u].push_back(v);
        rev[v].push_back(u);
        outdeg[u]++;
    }
    compute();
    print_pagerank();
    cout<<"Time taken to compute page rank of "<<N <<" nodes is"<< omp_get_wtime()-time<<endl;
    return 0;
}