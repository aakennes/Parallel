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
const int maxn=1e8+5;
ll sum,a[maxn],sum1,sum2;
long double ans1,ans2,ans3;
int n,m,cnt;
void Run1(){
    sum=0;
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    QueryPerformanceFrequency((LARGE_INTEGER *) &freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int i=0;i<n;++i){
        sum+=a[i];
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    ans1+=(tail-head)*1000.0/freq;
}
void Run2(){
    sum=sum1=sum2=0;
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    QueryPerformanceFrequency((LARGE_INTEGER *) &freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int i=0;i<n;i+=2){
        sum1+=a[i];
        sum2+=a[i+1];
    }
    sum=sum1+sum2;
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    ans2+=(tail-head)*1000.0/freq;
}
void Run3(){
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    QueryPerformanceFrequency((LARGE_INTEGER *) &freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int j=n;j>1;j/=2){
        for(int i=0;i<j/2;++i){
            a[i]=a[i*2]+a[i*2+1];
        }
    }
    sum=a[0];
    QueryPerformanceCounter((LARGE_INTEGER *)&tail);
    ans3+=(tail-head)*1000.0/freq;
}
int main(){
    freopen("1.in","r",stdin);
    freopen("1.out","w",stdout);
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    cin>>m>>cnt;
    n=pow(2,m);
    for(int i=1;i<=cnt;++i){
        Run1();
        Run2();
        Run3();
    }
    cout<<fixed<<setprecision(8)<<ans1/cnt<<" "<<ans2/cnt<<" "<<ans3/cnt<<endl;
    //printf("%ld ",time_usec);
    //printf("Elapsed time: %lf microseconds\n", time_usec);
}
