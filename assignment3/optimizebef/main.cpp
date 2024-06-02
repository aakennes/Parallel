#include"../params.h"
// #include"readwrite.h"
#include"ntt.h"

#include<cstdio>
#include<iostream>
#include<cstring>

const int maxn=1e6+5;

int p=Original_Q;
int n=Original_N;

int countw,wn[maxn];
int r[maxn];
int rt[maxn],irt[maxn];
Mint RT[maxn],IRT[maxn];

int a[maxn],b[maxn],ab[maxn];
Mint A[maxn],B[maxn],AB[maxn];
MMint aa[maxn],bb[maxn],aabb[maxn];

void poly_mul(){
    int aa[maxn],bb[maxn],aabb[maxn];
    memcpy(aa,a,sizeof(a));
    memcpy(bb,b,sizeof(b));
    memcpy(aabb,ab,sizeof(ab));
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            aabb[i+j]+=1LL*aa[i]*bb[j]%p;
            aabb[i+j]%=p;
        }
    } 
    for(int i=0;i<2*n-1;++i){
        std::cout<<aabb[i]<<" ";
    }
}
int nn[11]={256,512,1024,2048,4096,8192,16384,32768,65536,131072};
int qq[5]={1409,3329,7681,12289};
#include<sys/time.h>
#include<chrono>
void time_test(){
    freopen("1.out","w",stdout);
    for(int i=0;i<=9;++i){
        for(int j=0;j<=3;++j){
            long double ans=0;
            int cnt=100;
            for(int k=1;k<=100;++k){
                FILE* fp;
                // fRead(fp,a,b,nn[i],qq[j]);
                fRead(fp,A,B,nn[i],qq[j]);
                // std::cout<<nn[i]<<" "<<qq[j]<<'\n';
                auto Start=std::chrono::high_resolution_clock::now();
                // ntt_common(a,b,ab,r,nn[i]);
                // ntt_Montgomery(a,b,ab,r,nn[i]);
                // ntt_Montgomery_Mint(A,B,AB,r,nn[i]);
                // ntt_dif(a,b,ab,rt,irt,nn[i]);
                // ntt_dif_Mint(A,B,AB,RT,IRT,nn[i]);
                // ntt_dif_x4(a,b,ab,rt,irt,nn[i]);
                ntt_dif_x4_Mint(A,B,AB,RT,IRT,nn[i]);
                auto End=std::chrono::high_resolution_clock::now();
                std::chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
                ans+=elapsed.count();
                // std::cout<<elapsed.count()<<std::endl;
            }   
            std::cout<<ans/cnt<<std::endl;
        }
    }
    
}
int main(){
    time_test();
    // FILE* fp;
    // fRead(fp,A,B,nn[0],qq[0]);
    // ntt_dif_x4_Mint(A,B,AB,RT,IRT,nn[0]);
    // std::cout<<"asd\n";
    // ntt_dif_x4(a,b,ab,rt,irt,nn[0]);
    // ntt_dif(a,b,ab,rt,irt,nn[0]);
    // ntt_Montgomery_Mint(A,B,AB,r,nn[0]);
    // ntt_dif_Mint(A,B,AB,RT,IRT,nn[0]);
    // for(int i=0;i<10;++i)std::cout<<AB[i]<<" ";
    // // poly_mul();
    // // ntt_common(a,b,ab,r);
    // // ntt_dif(a,b,ab,rt,irt);
    // // ntt_dif_x4(a,b,ab,rt,irt);
    // // for(int i=0;i<nn[0];++i)std::cout<<a[i]<<" ";

    // for(int i=0;i<nn[0];++i)aa[i]=In((u32)a[i]);
    // for(int i=0;i<nn[0];++i)bb[i]=In((u32)b[i]);
    // std::cout<<aa[0]<<"\n";
    // ntt_Montgomery_MMint(aa,bb,aabb,r,nn[0]);
    // std::cout<<get(aabb[1])<<'\n';
    // for(int i=0;i<nn[0]*2-1;++i)std::cout<<get(aabb[i])<<" ";
    // fWrite(fp,ab,n);
    // puts("a--------------");
    // fRead(fp,A,B,nn[0],qq[0]);
    // std::cout<<A[0].getv()<<'\n';
    // ntt_Montgomery_Mint(A,B,AB,r,nn[0]);
    // std::cout<<AB[0]<<'\n';
    // ntt_dif_Mint(A,B,AB,RT,IRT);
    // ntt_dif_x4_Mint(A,B,AB,RT,IRT);
    // time_test();
    // fWrite(fp,AB,n);

    // fclose(fp);
    return 0;
}