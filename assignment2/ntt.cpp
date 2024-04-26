#include "ntt.h"

#define p Original_Q
#define n Original_N

using Mint = MontgomeryModInt32<p>; // 定义一个模 P 的模整数类型

void ntt_common(int *a,int *r,int limit,int type){
    for(int i = 0; i < limit; i++) {
		if(i < r[i]){
            std::swap(a[i], a[r[i]]);
        }
    }
    int wnn = 3;
    int invwnn = qpow(wnn, p-2, p);
	for(int mid = 1; mid < limit; mid <<= 1) {	
		int Wn = qpow( type == 1 ? wnn : invwnn , (p - 1) / (mid << 1),p);
		for(int j = 0; j < limit; j += (mid << 1)) {
			int w = 1;
			for(int k = 0; k < mid; k++, w = (1LL*w * Wn) % p) {
                int x = a[j + k], y = 1LL * w * a[j + k + mid] % p;             
                a[j + k] = (1LL*x + y) % p,
                a[j + k + mid] = (1LL* x - y + p) % p;
			}
		}
	}
}

void ntt_common(int *a,int *b,int *ab,int *r){
    
    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
	for(i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    ntt_common(a,r,limit, 1);
    ntt_common(b,r,limit, 1);
	for(i = 0; i < limit; i++){
        // printf("%d\n",b[i]);
        ab[i] = 1LL * a[i] * b[i] % p;
    } 
	ntt_common(ab,r,limit, -1);
    int invn=qpow(limit,p-2,p);

    // std::cout<<invn<<'\n';
    for(i = 0; i < 2*n; i++){
        ab[i] = (1LL * ab[i] * invn) % p;
    } 

}

void ntt_Montgomery(int *a,int *r,int limit,int type){

    for(int i = 0; i < limit; i++) {
		if(i < r[i]){
            std::swap(a[i], a[r[i]]);
        }
    }
    int wnn = 3;
    Mint Mwnn = wnn;
    Mint Minvwnn = 1/Mwnn;
	for(int mid = 1; mid < limit; mid <<= 1) {	
        Mint Moperator1 = type == 1 ? Mwnn : Minvwnn;
        int64_t Moperator2 = (p - 1) / (mid << 1);
		Mint MWn = Moperator1.pow(Moperator2);
		for(int j = 0; j < limit; j += (mid << 1)) {
			Mint Mw = 1;
			for(int k = 0; k < mid; k++, Mw = Mw * MWn) {                   
                Mint x = a[j + k];
                Mint Meven = a[j + k + mid];
				Mint y = Mw * Meven;                    
				a[j + k] = (x + y).get(),
				a[j + k + mid] = (x - y).get();
			}
		}
	}
}

void ntt_Montgomery(int *a,int *b,int *ab,int *r){

    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
	for(i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    ntt_Montgomery(a,r,limit, 1);
    ntt_Montgomery(b,r,limit, 1);
	for(i = 0; i < limit; i++){
        Mint Mai = a[i];
        Mint Mbi = b[i];
        ab[i] = (Mai * Mbi).get();
    } 
	ntt_Montgomery(ab,r,limit, -1);
    Mint Mlimit = limit;
    int invn = (1 / Mlimit).get();
    for(i = 0; i < 2*n; i++){
        Mint Mabi = ab[i];
        Mint Minvn = invn;
        ab[i] = (Mabi * Minvn).get();
    } 
}

void ntt_dit(int *a,int *r,int limit){
    for (int i = 1, l = limit >> 1; i < limit; i <<= 1, l >>= 1) {
        for (int j = 0, q = 0; j < i; j++, q += l << 1) {
            int w = r[i + j];
            for (int k = q; k < q + l; k++) {
                int x = a[k] , y = 1LL * a[k + l] * w % p;
                a[k + l] = (1LL * x - y + p)%p;
                a[k] = (1LL * x + y)%p;
            }
        }
    }
}

void ntt_dif(int *a,int *r,int limit){
    for (int l = 1, i = limit >> 1; i; l <<= 1, i >>= 1) {
        for (int j = 0, q = 0; j < i; j++, q += l << 1) {
            int w = r[i + j];
            for (int k = q; k < q + l; k++) {
                int x = a[k] , y = a[k + l];
                a[k + l] = (1LL * x - y + p)%p * w %p;
                a[k] = (1LL * x + y)%p;
            }
        }
    }
}

void ntt_dif(int *a,int *b,int *ab,int *rt,int *irt){
    int L=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    int wnn=3;
    int invwnn=qpow(wnn,p-2,p);
    rt[0] = 1;
    rt[1] = qpow(wnn, (p - 1) / limit / 2,p);
    for (int i = 2; i < limit; i++) rt[i] = 1LL * rt[i-1] * rt[1] % p;
    irt[0] = 1;
    irt[1] = qpow(rt[1], p - 2,p);
    for (int i = 2; i < limit; i++) irt[i] = 1LL * irt[i-1] * irt[1] % p;
    for (int i = 0, j = 0; i < limit; i++) {
        if (i > j) {
            std::swap(rt[i], rt[j]);
            std::swap(irt[i], irt[j]);
        }
        printf("%d ",rt[i]);
        for (int t = limit >> 1; (j ^= t) < t; t >>= 1);
    }
    ntt_dit(a,rt,limit);
    ntt_dit(b,rt,limit);
	for(int i = 0;i<limit;i++){
        // printf("%d ",b[i]);
        ab[i]=1LL*a[i]*b[i]%p;
    }
	ntt_dif(ab,irt,limit);
    int invn=qpow(limit,p-2,p);
    // std::cout<<invn<<'\n';
    for(int i = 0; i < 2 * n; i++){
        ab[i] = (1LL * ab[i] * invn) % p;
        printf("%d\n",ab[i]);
    } 
}