#include"params.h"
#include"readwrite.h"
#include"Math.h"
#include"ntt.h"

#include<cstdio>
#include<iostream>
#include<cstring>

const int maxn=9e3+5;

int p=Original_Q;
int n=Original_N;

int countw,wn[maxn];
int r[maxn];



int a[maxn],b[maxn],ab[maxn];
void poly_mul(){
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            ab[i+j]+=1LL*a[i]*b[j]%p;
            ab[i+j]%=p;
        }
    }
    for(int i=0;i<2*n;++i){
        std::cout<<ab[i]<<" ";
    }
}
int main(){
    findw(countw,wn);
    FILE* fp;
    fRead(fp,a,b,n);
    // poly_mul();
    ntt_common(a,b,ab,r);
    fWrite(fp,ab,n);
    
    // fclose(fp);
    return 0;
}