#include<cstdio>
#include<iostream>
#include<sys/time.h>
#include<istream>
#include<ostream>
#include<fstream>
#include<chrono>
#include<iomanip>
#include<windows.h>
#define ll long long
using namespace std;
const int maxn=1e4+5;
ll sum[maxn],b[maxn][maxn],a[maxn];
int n;
int main(){
    freopen("1.in","r",stdin);
    //freopen("1.out")
    cin>>n;
    for(int i=0;i<n;++i){
        a[i]=i+1;
        for(int j=0;j<n;++j){
            b[i][j]=i+j+1;
        }
    }
    ll head,freq,tail;
    QueryPerformanceFrequency((LARGE_INTEGER *) &freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int i=0;i<n;++i){
        sum[i]=0.0;
    }
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            sum[i]+=b[j][i]*a[j];
        }
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    ofstream out;
    out.open("1.out",ios::app);
    out<<fixed<<setprecision(8)<<(tail-head)*1000.0/freq<<" "; 
    out.close();
    //printf("%ld ",time_usec);
    //printf("Elapsed time: %lf microseconds\n", time_usec);
}
