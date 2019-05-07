#include<bits/stdc++.h>
#include<iostream>
#include<omp.h>
using namespace std;
vector <int> adj[1000000];
vector <int> rev[1000000];
int vis[1000000],scc[1000000];
int level[1000000];
vector <int> LEVEL[200];
vector <int> SCC[1000000];
deque <int> L;
int outdeg[1000000];
int N,M;
double d=0.15;
double threshold=0.00000001;
int SCC_MAX;
double pr[1000000];
double prev1[1000000];
int dead[1000000];
int block_indeg[1000000];
int MAX_LL;
vector <int> final_adj[1000000];
deque <int> block_zero;
vector <int> block_topo;
vector <int> block_adj[1000000];

void dfs1(int u){
    if(vis[u]==0){
        vis[u]=1;
        for(int i=0;i<adj[u].size();i++)
        {
            dfs1(adj[u][i]);
        }
        L.push_front(u);
    }
}
void print_pagerank(){
    cout<<"Final Page rank"<<endl;
    for(int i = 1; i <= N; i++){
        cout<<pr[i]<<endl;
    }
}
void assign(int u,int root){
    if(scc[u]==0){
        scc[u]=root;
        SCC[root].push_back(u);
        for(int i=0;i<rev[u].size();i++){
            assign(rev[u][i],root);
        }
    }
}
void calculate_scc(){
    for(int i=1;i<=N;i++){
        if(!vis[i])
            dfs1(i);
    }
    
    int SCC_count=1;
    deque<int>::iterator it = L.begin();
    while(it!=L.end()){
        if(!scc[*it]){
         
            assign(*it,SCC_count++);
        }
        it++;
    }
    SCC_MAX=SCC_count;
}

void dfs2(int u){
    if(vis[u]==2)
        return ;
    vis[u]=2;
    for(int i=0;i<adj[u].size();i++){
        if(scc[u]!=scc[adj[u][i]]){
            block_adj[scc[u]].push_back(scc[adj[u][i]]);
            block_indeg[scc[adj[u][i]]]++;
        }
        if(vis[adj[u][i]]!=2)
            dfs2(adj[u][i]);
    }
}


void topological_ordering_scc(){
    for(int i=1;i<=N;i++){
        if(vis[i]!=2)
            dfs2(i);
    }
    int count=0;
    for(int i=1;i<SCC_MAX;i++)
    {
        if(block_indeg[i]==0){
            block_zero.push_back(i);
            level[i]=1;
            ++count;
            LEVEL[1].push_back(i);
        }
    }

    int max_lev=1;
    while(!block_zero.empty()){
        int u=block_zero.front();
        block_zero.pop_front();
        block_topo.push_back(u);
        for(int i=0;i<block_adj[u].size();i++){
            block_indeg[block_adj[u][i]]--;
            if(block_indeg[block_adj[u][i]]==0){
                block_zero.push_back(block_adj[u][i]);
                level[block_adj[u][i]]=1+level[u];
                LEVEL[1+level[u]].push_back(block_adj[u][i]);
                max_lev=max(max_lev,1+level[u]);
            }
        }
    }
    MAX_LL=max_lev;
}
void calculate_pgrank_SCC(int sc){
    double error=10000000000.0;
    double threshold2=threshold/N;
    // #pragma omp parallel for
    for(int i=0;i<SCC[sc].size();i++){
        prev1[SCC[sc][i]]=d/N;
    }
    while(error>threshold){
        //#pragma omp parallel for
        for(int i=0;i<SCC[sc].size();i++){
            if(dead[SCC[sc][i]]==1)
               continue;
            pr[SCC[sc][i]]=d/(N*1.0);
            for(int v=0;v<rev[SCC[sc][i]].size();v++){
                int ver=rev[SCC[sc][i]][v];
                pr[SCC[sc][i]]+=((prev1[ver]*1.0)/(outdeg[ver]*1.0))*(1.0-d);
            }

        }
        double nerr=0.0;
        for(int i=0;i<SCC[sc].size();i++){
            nerr=max(nerr,abs(prev1[SCC[sc][i]]-pr[SCC[sc][i]]));
            if(abs(prev1[SCC[sc][i]]-pr[SCC[sc][i]]<=threshold2))
                dead[SCC[sc][i]]=1;
            prev1[SCC[sc][i]]=pr[SCC[sc][i]];
        }
        error=nerr;
        
    }
}


void compute_pagerank(){
    for(int level=1;level<=MAX_LL;level++){
       //#pragma omp parallel for
       for(int i=0;i<LEVEL[level].size();i++){
         calculate_pgrank_SCC(LEVEL[level][i]);
        } 
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
        if(u!=v){
            adj[u].push_back(v);
            rev[v].push_back(u);
            outdeg[u]++;
        }
    }
    calculate_scc();
    topological_ordering_scc();
    compute_pagerank();
    print_pagerank();
    cout<<"Time taken to compute page rank of "<<N <<" nodes is"<< omp_get_wtime()-time<<endl;
    return 0;
}