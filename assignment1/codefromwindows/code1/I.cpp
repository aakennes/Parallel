#include<cstdio>
#include<iostream>
#include<sys/time.h>
#include<istream>
#include<ostream>
#include<fstream>
#include<chrono>
#include<iomanip>
#include<cmath>
#include<windows.h>
#define ll long long
using namespace std;
const int maxn=1e4+5;
ll sum[maxn],b[maxn][maxn],a[maxn];
int n,cnt;
long double ans1,ans2;
void Init(){
    for(int i=0;i<n;++i){
        a[i]=i+1;
        for(int j=0;j<n;++j){
            b[i][j]=i+j+1;
        }
    }
    for(int i=0;i<n;++i){
        sum[i]=0.0;
    }
}
void Run1(){
    Init();
    ll head,freq,tail;
    QueryPerformanceFrequency((LARGE_INTEGER *) &freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            sum[i]+=b[j][i]*a[j];
        }
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    ans1+=(tail-head)*1000.0/freq;
}
void Run2(){
    Init();
    ll head,freq,tail;
    QueryPerformanceFrequency((LARGE_INTEGER *) &freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int j=0;j<n;++j){
        for(int i=0;i<n;++i){
            sum[i]+=b[j][i]*a[j];
        }
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    ans2+=(tail-head)*1000.0/freq;
}
int main(){
    freopen("1.in","r",stdin);
    freopen("1.out","w",stdout);
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    cin>>n>>cnt;
    for(int i=1;i<=cnt;++i){
        Run1();
        Run2();
    }
    cout<<fixed<<setprecision(8)<<ans1/cnt<<" "<<ans2/cnt<<endl;
    //printf("%ld ",time_usec);
    //printf("Elapsed time: %lf microseconds\n", time_usec);
}
