#include"params.h"
#include"readwrite.h"
// #include"ntt.h"

#include<cstdio>
#include<iostream>
#include<cstring>

const int maxn=1e6+5;

int p=Original_Q;
int n=Original_N;

int countw,wn[maxn];
int r[maxn];
int rt[maxn],irt[maxn];



int a[maxn],b[maxn],ab[maxn];
Mint A[maxn],B[maxn],AB[maxn];
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
int main(){
    // findw(countw,wn);
    FILE* fp;
    // fRead(fp,a,b,n);
    // // poly_mul();
    // // ntt_common(a,b,ab,r);
    // // ntt_dif(a,b,ab,rt,irt);
    // // ntt_dif_x4(a,b,ab,rt,irt);
    // // ntt_Montgomery(a,b,ab,r);
    // fWrite(fp,ab,n);


    fRead(fp,A,B,n);
    ntt_Montgomery_Mint(A,B,AB,r);
    fWrite(fp,AB,n);
    // fclose(fp);
    return 0;
}