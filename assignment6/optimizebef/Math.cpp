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