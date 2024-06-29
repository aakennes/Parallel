#include"Math.h"

int qpow(int a,int x,int p){
    int ans=1;
    while(x){
        if(x&1)ans=1LL*ans*a%p;
        a=1LL*a*a%p;
        x>>=1;
    }
    return ans%p;
}

void findw(int countw,int wn[])
{
    int w, i, q = Original_Q, n = Pt3_N;
    countw = 0;
    for (w = 2; w < q; w++){
        for (i = 0; i <= n; i++){
            if (qpow(w, i, q) == q - 1)
                if (i != n) continue;
                else //printf("%4d,",wn[countw++] = w);
                    wn[countw++] = w;
        }
    }
}

int barrett_mod(int a,int b){
    u128 y = (((u128)1 << 64) / b)>> 64;
    u64 x = a - ((u128)a * y) * b;
    // std::cout<<a<<" "<<b<<" "<<x<<'\n';
    while(x >= Original_Q)x-=Original_Q;
    // std::cout<<a<<" "<<b<<" "<<x<<'\n';
    return x;
}

int barrett_add(int a,int b,int p){
    i128 aa = (i64)a+b;
    // u128 y = (((u128)1 << 64) / p) >> 64;
    u64 x = aa - (aa * barret_y >> 64) * p;
    if(x >= p)x -= p;
    return x;
}

int barrett_sub(int a,int b,int p){
    i64 aa = a-b;
    if(aa < 0)aa += p;
    return aa;
}

int barrett_mul(int a,int b,int p){
    i128 aa = (i64)a*b;
    u64 x = aa - (aa * barret_y >> 64) * p;
    if(x >= p)x -= p;
    return x;
}

#define BARRETT_MUL(a, b, p) ({ \
    i128 aa = (i64)(a) * (b); \
    u64 x = aa - ((aa * (barret_y)) >> 64) * (p); \
    if (x >= (p)) x -= (p); \
    x; \
})
int barrett_qpow(int a,int x,int p){
    int ans=1;
    while(x){
        if(x&1)ans=BARRETT_MUL(ans,a,p);
        // a=1LL*a*a%p;
        a=BARRETT_MUL(a,a,p);
        x>>=1;
    }
    return ans;
}
#undef BARRETT_MUL