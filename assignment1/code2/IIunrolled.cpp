#include<cstdio>
#include<iostream>
#include<sys/time.h>
#include<istream>
#include<ostream>
#include<fstream>
#include<chrono>
#include<iomanip>
#include<cmath>
#define ll long long
using namespace std;
const int maxn=1e8+5;
ll sum,a[maxn],sum1,sum2,sum3,sum4,sum5;
long double ans1,ans2,ans3,ans4,ans5,ans6,ans7,ans8;
int n,m,cnt;
void Run1(){//1*1
    sum=0;
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    auto Start=chrono::high_resolution_clock::now();
    for(int i=0;i<n;++i){
        sum+=a[i];
    }
    auto End=chrono::high_resolution_clock::now();
    chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
    ans1+=elapsed.count();
}
void Run2(){//2*2
    sum=sum1=sum2=0;
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    auto Start=chrono::high_resolution_clock::now();
    for(int i=0;i<n;i+=2){
        sum1+=a[i];
        sum2+=a[i+1];
    }
    sum=sum1+sum2;
    auto End=chrono::high_resolution_clock::now();
    chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
    ans2+=elapsed.count();
}
void Run3(){
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    auto Start=chrono::high_resolution_clock::now();
    for(int j=n;j>1;j/=2){
        for(int i=0;i<j/2;++i){
            a[i]=a[i*2]+a[i*2+1];
        }
    }
    sum=a[0];
    auto End=chrono::high_resolution_clock::now();
    chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
    ans3+=elapsed.count();
}
void Run4(){//3*3
    sum=sum1=sum2=sum3=0;
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    auto Start=chrono::high_resolution_clock::now();
    for(int i=0;i<n;i+=3){
        sum1+=a[i];
        sum2+=a[i+1];
        sum3+=a[i+2];
    }
    sum=sum1+sum2+sum3;
    auto End=chrono::high_resolution_clock::now();
    chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
    ans4+=elapsed.count();
}
void Run5(){//2*1
    sum=0;
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    auto Start=chrono::high_resolution_clock::now();
    for(int i=0;i<n;i+=2){
        sum+=a[i];
        sum+=a[i+1];
    }
    auto End=chrono::high_resolution_clock::now();
    chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
    ans5+=elapsed.count();
}
void Run6(){//3*1
    sum=0;
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    auto Start=chrono::high_resolution_clock::now();
    for(int i=0;i<n;i+=3){
        sum+=a[i];
        sum+=a[i+1];
        sum+=a[i+2];
    }
    auto End=chrono::high_resolution_clock::now();
    chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
    ans6+=elapsed.count();
}
void Run7(){//4*4
    sum=sum1=sum2=sum3=sum4=0;
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    auto Start=chrono::high_resolution_clock::now();
    for(int i=0;i<n;i+=4){
        sum1+=a[i];
        sum2+=a[i+1];
        sum3+=a[i+2];
        sum4+=a[i+3];
    }
    sum=sum1+sum2+sum3+sum4;
    auto End=chrono::high_resolution_clock::now();
    chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
    ans7+=elapsed.count();
}
void Run8(){//3*3
    sum=sum1=sum2=sum3=sum4=sum5=0;
    for(int i=0;i<n;++i){
        a[i]=i*10+1;
    }
    ll head,freq,tail;
    auto Start=chrono::high_resolution_clock::now();
    for(int i=0;i<n;i+=5){
        sum1+=a[i];
        sum2+=a[i+1];
        sum3+=a[i+2];
        sum4+=a[i+3];
        sum5+=a[i+4];
    }
    sum=sum1+sum2+sum3+sum4+sum5;
    auto End=chrono::high_resolution_clock::now();
    chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
    ans8+=elapsed.count();
}
int main(){
    freopen("1.in","r",stdin);
    freopen("1.out","w",stdout);
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    cin>>m>>cnt;
    n=pow(2,m);
    for(int i=1;i<=cnt;++i){
        //Run1();//1*1
        //Run2();//2*2
        //Run3();
        Run4();//3*3
        Run5();//2*1
        Run6();//3*1
        Run7();//4*4
        Run8();//5*5
    }//1*1 2*1 3*1 2*2 3*3 4*4 5*5
    cout<<fixed<<setprecision(8)<<ans5/cnt<<" "<<ans6/cnt<<" "<<ans4/cnt<<" "<<ans7/cnt<<" "<<ans8/cnt<<endl;
    //printf("%ld ",time_usec);
    //printf("Elapsed time: %lf microseconds\n", time_usec);
}
