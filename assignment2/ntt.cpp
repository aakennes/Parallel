#include "ntt.h"

#define p Original_Q
#define n Original_N

void ntt_common(int *a,int *r,int limit,int type){
    for(int i = 0; i < limit; i++) {
		if(i < r[i]){
            std::swap(a[i], a[r[i]]);
        }
    }
    int wnn=3;
    int invwnn=qpow(wnn,p-2,p);
    // std::cout<<wnn<<" "<<invwnn<<"\n";
	for(int mid = 1; mid < limit; mid <<= 1) {	
		int Wn = qpow( type == 1 ? wnn : invwnn , (p - 1) / (mid << 1),p);
		for(int j = 0; j < limit; j += (mid << 1)) {
			int w = 1;
			for(int k = 0; k < mid; k++, w = (w * Wn) % p) {
				 int x = a[j + k], y = 1LL * w * a[j + k + mid] % p;
				 a[j + k] = (x + y) % p,
				 a[j + k + mid] = (x - y + p) % p;
			}
		}
	}
    // std::cout<<wnn<<"\n";
}

void ntt_common(int *a,int *b,int *ab,int *r){
    
    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    // std::cout<<L<<" "<<limit<<'\n';
	for(i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    ntt_common(a,r,limit, 1);
    ntt_common(b,r,limit, 1);
	for(i = 0; i < limit; i++){
        // printf("%d\n",a[i]);
        ab[i] = (1LL * a[i] * b[i]) % p;
    } 
	ntt_common(ab,r,limit, -1);
    int invn=qpow(limit,p-2,p);
    for(i = 0; i < 2*n; i++){
        // printf("%d\n",a[i]);
        ab[i] = (1LL * ab[i] * invn) % p;
    } 
}