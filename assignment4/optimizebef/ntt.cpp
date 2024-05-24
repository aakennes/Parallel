#include"Math.h"

#define p Original_Q
int W[100005];
const int thread_count = 2;

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
void ntt_common(int *a,int *b,int *ab,int *r,int n){
    
    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    // printf("%d\n",limit);
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
void ntt_common_openmp(int *a,int *r,int limit,int type){
    for(int i = 0; i < limit; i++) {
		if(i < r[i]){
            std::swap(a[i], a[r[i]]);
        }
    }
    int wnn = 3;
    int invwnn = qpow(wnn, p-2, p);
    #pragma omp parallel num_threads(thread_count)

    for(int mid = 1; mid < limit; mid <<= 1) {	
        int Wn = qpow(type == 1 ? wnn : invwnn, (p - 1) / (mid << 1), p);
        // init
        W[0] = 1;
        for(int k = 1; k < mid; ++k){
            W[k] = 1LL * W[k-1] * Wn % p;
        }

        #pragma omp for schedule(guided)
        for(int j = 0; j < limit; j += (mid << 1)) {
            // #pragma omp for schedule(guided)
            for(int k = 0; k < mid; k++) {
                int x = a[j + k], y = 1LL * W[k] * a[j + k + mid] % p;             
                a[j + k] = (1LL*x + y) % p,
                a[j + k + mid] = (1LL* x - y + p) % p;
            }
        }
    }
    
}
void ntt_common_openmp(int *a,int *b,int *ab,int *r,int n){
    
    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    // printf("%d\n",limit);
	for(i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    ntt_common_openmp(a,r,limit, 1);
    ntt_common_openmp(b,r,limit, 1);
	for(i = 0; i < limit; i++){
        // printf("%d\n",b[i]);
        ab[i] = 1LL * a[i] * b[i] % p;
    } 
	ntt_common_openmp(ab,r,limit, -1);
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
void ntt_Montgomery(int *a,int *b,int *ab,int *r,int n){

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
void ntt_Montgomery_openmp(int *a,int *r,int limit,int type){
    
    for(int i = 0; i < limit; i++) {
		if(i < r[i]){
            std::swap(a[i], a[r[i]]);
        }
    }
    int wnn = 3;
    Mint Mwnn = wnn;
    Mint Minvwnn = 1/Mwnn;

    #pragma omp parallel num_threads(thread_count)

	for(int mid = 1; mid < limit; mid <<= 1) {	
        Mint Moperator1 = type == 1 ? Mwnn : Minvwnn;
        int64_t Moperator2 = (p - 1) / (mid << 1);
		Mint MWn = Moperator1.pow(Moperator2);
        #pragma omp for schedule(guided)
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
void ntt_Montgomery_openmp(int *a,int *b,int *ab,int *r,int n){

    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
	for(i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    ntt_Montgomery_openmp(a,r,limit, 1);
    ntt_Montgomery_openmp(b,r,limit, 1);
	for(i = 0; i < limit; i++){
        Mint Mai = a[i];
        Mint Mbi = b[i];
        ab[i] = (Mai * Mbi).get();
    } 
	ntt_Montgomery_openmp(ab,r,limit, -1);
    Mint Mlimit = limit;
    int invn = (1 / Mlimit).get();
    for(i = 0; i < 2*n; i++){
        Mint Mabi = ab[i];
        Mint Minvn = invn;
        ab[i] = (Mabi * Minvn).get();
    } 
}

void ntt_Montgomery_Mint(Mint *a,int *r,int limit,int type){
    
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
				Mint y = Mw * a[j + k + mid];                    
				a[j + k] = (x + y),
				a[j + k + mid] = (x - y);
			}
		}
	}
}
void ntt_Montgomery_Mint(Mint *a,Mint *b,Mint *ab,int *r,int n){

    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
	for(i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    ntt_Montgomery_Mint(a,r,limit, 1);
    ntt_Montgomery_Mint(b,r,limit, 1);
	for(i = 0; i < limit; i++){
        ab[i]= a[i] * b[i];
    } 
	ntt_Montgomery_Mint(ab,r,limit, -1);
    Mint Mlimit = limit;
    Mint invn = (1 / Mlimit);
    for(i = 0; i < 2*n; i++){
        ab[i] = (ab[i] * invn);
    } 
}
void ntt_Montgomery_Mint_openmp(Mint *a,int *r,int limit,int type){
    
    for(int i = 0; i < limit; i++) {
		if(i < r[i]){
            std::swap(a[i], a[r[i]]);
        }
    }
    int wnn = 3;
    Mint Mwnn = wnn;
    Mint Minvwnn = 1/Mwnn;
    #pragma omp parallel num_threads(thread_count)
	for(int mid = 1; mid < limit; mid <<= 1) {	
        Mint Moperator1 = type == 1 ? Mwnn : Minvwnn;
        int64_t Moperator2 = (p - 1) / (mid << 1);
		Mint MWn = Moperator1.pow(Moperator2);

        #pragma omp for schedule(guided)
		for(int j = 0; j < limit; j += (mid << 1)) {
			Mint Mw = 1;
			for(int k = 0; k < mid; k++, Mw = Mw * MWn) {                   
                Mint x = a[j + k];
				Mint y = Mw * a[j + k + mid];                    
				a[j + k] = (x + y),
				a[j + k + mid] = (x - y);
			}
		}
	}
}
void ntt_Montgomery_Mint_openmp(Mint *a,Mint *b,Mint *ab,int *r,int n){

    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
	for(i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    ntt_Montgomery_Mint_openmp(a,r,limit, 1);
    ntt_Montgomery_Mint_openmp(b,r,limit, 1);
	for(i = 0; i < limit; i++){
        ab[i]= a[i] * b[i];
    } 
	ntt_Montgomery_Mint_openmp(ab,r,limit, -1);
    Mint Mlimit = limit;
    Mint invn = (1 / Mlimit);
    for(i = 0; i < 2*n; i++){
        ab[i] = (ab[i] * invn);
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
void ntt_dif(int *a,int *b,int *ab,int *rt,int *irt,int n){
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
    // for (int i = 2; i < limit; i++) printf("%d ", rt[i]);
    irt[0] = 1;
    irt[1] = qpow(rt[1], p - 2,p);
    for (int i = 2; i < limit; i++) irt[i] = 1LL * irt[i-1] * irt[1] % p;
    // for (int i = 2; i < limit; i++) printf("%d ", irt[i]);
    for (int i = 0, j = 0; i < limit; i++) {
        if (i > j) {
            std::swap(rt[i], rt[j]);
            std::swap(irt[i], irt[j]);
        }
        // printf("%d ",rt[i]);
        for (int t = limit >> 1; (j ^= t) < t; t >>= 1);
    }
    ntt_dit(a,rt,limit);
    ntt_dit(b,rt,limit);
	for(int i = 0;i<limit;i++){
        ab[i]=1LL*a[i]*b[i]%p;
        // printf("%d ",ab[i]);
    }
	ntt_dif(ab,irt,limit);
    int invn=qpow(limit,p-2,p);
    // std::cout<<invn<<'\n';
    for(int i = 0; i < 2 * n; i++){
        ab[i] = (1LL * ab[i] * invn) % p;
        // printf("%d\n",ab[i]);
    } 
}
void ntt_dit_openmp(int *a, int *r, int limit) {
    for (int i = 1, l = limit >> 1; i < limit; i <<= 1, l >>= 1) {
        #pragma omp parallel for num_threads(thread_count) schedule(static)
        for (int j = 0; j < i; j++) {
            int w = r[i + j];
            for (int k = j * (l << 1); k < j * (l << 1) + l; k++) {
                int x = a[k], y = 1LL * a[k + l] * w % p;
                a[k + l] = (1LL * x - y + p) % p;
                a[k] = (1LL * x + y) % p;
            }
        }
    }
}
void ntt_dif_openmp(int *a, int *r, int limit) {
    for (int l = 1, i = limit >> 1; i; l <<= 1, i >>= 1) {
        #pragma omp parallel for num_threads(thread_count) schedule(static)
        for (int j = 0; j < i; j++) {
            int w = r[i + j];
            for (int k = j * (l << 1); k < j * (l << 1) + l; k++) {
                int x = a[k], y = a[k + l];
                a[k + l] = 1LL * (x - y + p) * w % p;
                a[k] = (x + y) % p;
            }
        }
    }
}
void ntt_dif_openmp(int *a,int *b,int *ab,int *rt,int *irt,int n){
    int limit = 1;
    while (limit <= 2 * n - 2) {
        limit <<= 1;
    }

    int wnn = 3;
    int invwnn = qpow(wnn, p-2, p);

    rt[0] = 1;
    rt[1] = qpow(wnn, (p - 1) / limit / 2, p);
    for (int i = 2; i < limit; i++) rt[i] = 1LL * rt[i-1] * rt[1] % p;

    irt[0] = 1;
    irt[1] = qpow(rt[1], p - 2, p);
    for (int i = 2; i < limit; i++) irt[i] = 1LL * irt[i-1] * irt[1] % p;

    for (int i = 0, j = 0; i < limit; i++) {
        if (i > j) {
            std::swap(rt[i], rt[j]);
            std::swap(irt[i], irt[j]);
        }
        for (int t = limit >> 1; (j ^= t) < t; t >>= 1);
    }

    ntt_dit_openmp(a, rt, limit);
    ntt_dit_openmp(b, rt, limit);

    for (int i = 0; i < limit; i++) {
        ab[i] = 1LL * a[i] * b[i] % p;
    }

    ntt_dif_openmp(ab, irt, limit);

    int invn = qpow(limit, p - 2, p);
    for (int i = 0; i < 2 * n; i++) {
        ab[i] = 1LL * ab[i] * invn % p;
    }
}

void ntt_dit_Mint(Mint *a,Mint *r,int limit){
    for (int i = 1, l = limit >> 1; i < limit; i <<= 1, l >>= 1) {
        for (int j = 0, q = 0; j < i; j++, q += l << 1) {
            Mint w = r[i + j];
            for (int k = q; k < q + l; k++) {
                Mint x = a[k] , y = a[k + l] * w;
                a[k + l] = x - y;
                a[k] = x + y;
            }
        }
    }
}
void ntt_dif_Mint(Mint *a,Mint *r,int limit){
    for (int l = 1, i = limit >> 1; i; l <<= 1, i >>= 1) {
        for (int j = 0, q = 0; j < i; j++, q += l << 1) {
            Mint w = r[i + j];
            for (int k = q; k < q + l; k++) {
                Mint x = a[k] , y = a[k + l];
                a[k + l] =(x - y) * w;
                a[k] = (x + y);
            }
        }
    }
}
void ntt_dif_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n){
    int L=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    Mint wnn = 3;
    Mint invwnn = 1 / wnn;
    rt[0] = 1;
    rt[1] = wnn.pow((p - 1) / limit / 2);
    for (int i = 2; i < limit; i++) rt[i] = rt[i-1] * rt[1];
    // for (int i = 2; i < limit; i++) printf("%d ", rt[i].get());
    irt[0] = 1;
    irt[1] = 1 / rt[1];
    for (int i = 2; i < limit; i++) irt[i] = irt[i-1] * irt[1];
    // for (int i = 2; i < limit; i++) printf("%d ", irt[i].get());
    for (int i = 0, j = 0; i < limit; i++) {
        if (i > j) {
            std::swap(rt[i], rt[j]);
            std::swap(irt[i], irt[j]);
        }
        // printf("%d ",rt[i]);
        for (int t = limit >> 1; (j ^= t) < t; t >>= 1);
    }
    ntt_dit_Mint(a,rt,limit);
    ntt_dit_Mint(b,rt,limit);
	for(int i = 0;i<limit;i++){
        // printf("%d ",b[i].get());
        ab[i]=a[i]*b[i];
        // printf("%d ",ab[i].get());
    }
	ntt_dif_Mint(ab,irt,limit);
    Mint Mlimit = limit;
    Mint invn = 1 / Mlimit;
    // std::cout<<invn<<'\n';
    for(int i = 0; i < 2 * n; i++){
        ab[i] = ab[i] * invn;
        // printf("%d\n",ab[i].get());
    } 
}
void ntt_dit_Mint_openmp(Mint *a,Mint *r,int limit){
    #pragma omp parallel num_threads(thread_count)
    {
        for (int i = 1, l = limit >> 1; i < limit; i <<= 1, l >>= 1) {
            #pragma omp for schedule(static)        
            for (int j = 0; j < i; j++){
                Mint w = r[i + j];
                for (int k = j * (l << 1); k < j * (l << 1) + l; k++){
                    Mint x = a[k] , y = a[k + l] * w;
                    a[k + l] = x - y;
                    a[k] = x + y;
                }
            }
        }
    }

}
void ntt_dif_Mint_openmp(Mint *a,Mint *r,int limit){
    #pragma omp parallel num_threads(thread_count)
    {
        for (int l = 1, i = limit >> 1; i; l <<= 1, i >>= 1) {
            #pragma omp for schedule(static)
            for (int j = 0; j < i; j++) {
                Mint w = r[i + j];
                for (int k = j * (l << 1); k < j * (l << 1) + l; k++) {
                    Mint x = a[k] , y = a[k + l];
                    a[k + l] =(x - y) * w;
                    a[k] = (x + y);
                }
            }
        }
    }
}
void ntt_dif_Mint_openmp(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n){
    int L=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    Mint wnn = 3;
    Mint invwnn = 1 / wnn;
    rt[0] = 1;
    rt[1] = wnn.pow((p - 1) / limit / 2);
    for (int i = 2; i < limit; i++) rt[i] = rt[i-1] * rt[1];
    // for (int i = 2; i < limit; i++) printf("%d ", rt[i].get());
    irt[0] = 1;
    irt[1] = 1 / rt[1];
    for (int i = 2; i < limit; i++) irt[i] = irt[i-1] * irt[1];
    // for (int i = 2; i < limit; i++) printf("%d ", irt[i].get());
    for (int i = 0, j = 0; i < limit; i++) {
        if (i > j) {
            std::swap(rt[i], rt[j]);
            std::swap(irt[i], irt[j]);
        }
        // printf("%d ",rt[i]);
        for (int t = limit >> 1; (j ^= t) < t; t >>= 1);
    }
    ntt_dit_Mint_openmp(a,rt,limit);
    ntt_dit_Mint_openmp(b,rt,limit);
    #pragma omp parallel for num_threads(thread_count) schedule(static)
	for(int i = 0;i<limit;i++){
        // printf("%d ",b[i].get());
        ab[i]=a[i]*b[i];
        // printf("%d ",ab[i].get());
    }
	ntt_dif_Mint_openmp(ab,irt,limit);
    Mint Mlimit = limit;
    Mint invn = 1 / Mlimit;
    // std::cout<<invn<<'\n';
    #pragma omp parallel for num_threads(thread_count) schedule(static)
    for(int i = 0; i < 2 * n; i++){
        ab[i] = ab[i] * invn;
        // printf("%d\n",ab[i].get());
    } 
}

void ntt_dit_x4(int *a,int *rt,int *irt,int limit,int level){
    int one = 1, imag = rt[level-2];
    int logn = __builtin_ctz(limit);
    int dw[level - 1];
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL*dw[i - 1] * irt[level - 1 - i]%p * rt[level - 3 - i]%p;//,printf("%d ",dw[i]);
    dw[level - 2] = 1LL*dw[level - 3] * irt[1]%p;
    if (logn & 1) {
        for (int i = 0; i < limit/2; i++) {
            int x = a[i], y = a[i + limit/2];
            a[i] = (1LL*x + y)% p;
            a[i + limit/2] = ((1LL*x - y)%p + p)%p;
        }
    }
    for (int e = logn & ~1; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        int w2 = 1;
        for (int i = 0; i < limit; i += m) {
            const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2%p;
            for (int j = i; j < i + m4; ++j) {
                int a0 = 1LL*a[j + m4 * 0] * one%p, a1 = 1LL*a[j + m4 * 1] * w2%p;
                int a2 = 1LL*a[j + m4 * 2] * w1%p, a3 = 1LL*a[j + m4 * 3] * w3%p;
                int t02p = (1LL*a0 + a2)%p, t13p = (1LL*a1 + a3)%p;
                int t02m = ((1LL*a0 - a2)%p + p)%p, t13m = ((1LL*a1 - a3)%p +p)%p * imag%p;
                a[j + m4 * 0] = (1LL*t02p + t13p)%p;
                a[j + m4 * 1] = ((1LL*t02p - t13p)%p+p)%p;
                a[j + m4 * 2] = (1LL*t02m + t13m)%p;
                a[j + m4 * 3] = ((1LL*t02m - t13m)%p+p)%p;
            }
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
    }
}
void ntt_dif_x4(int *a,int *rt,int *irt,int limit,int level){
    int one = 1, imag = irt[level-2];
    int logn = __builtin_ctz(limit);
    int dw[level - 1];
    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL*dw[i - 1] * rt[level - 1 - i]%p * irt[level - 3 - i]%p;
    dw[level - 2] = 1LL*dw[level - 3] * rt[1]%p;
    
    for (int e = 2; e <= logn; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        int w2 = one;
        for (int i = 0; i < limit; i += m) {
            const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2 %p;
            for (int j = i; j < i + m4; ++j) {
                int a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                int a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                int t01p = (1LL*a0 + a1)%p, t23p = (1LL*a2 + a3)%p;
                int t01m = ((1LL*a0 - a1)%p + p)%p, t23m = ((1LL*a2 - a3)%p +p)%p * imag%p;
                a[j + m4 * 0] = (1LL*t01p + t23p)%p*one%p;
                a[j + m4 * 2] = ((1LL*t01p - t23p)%p+p)%p*w1%p;
                a[j + m4 * 1] = (1LL*t01m + t23m)%p*w2%p;
                a[j + m4 * 3] = ((1LL*t01m - t23m)%p+p)%p*w3%p;
            }
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        // printf("%d ",w2);
    }
    if (logn & 1) {
        for (int i = 0; i < limit/2; i++) {
            int x = a[i], y = a[i + limit/2];
            a[i] = (1LL*x + y )% p;
            a[i + limit/2] = ((1LL*x - y)%p + p)%p;
        }
    }
}
void ntt_dif_x4(int *a,int *b,int *ab,int *rt,int *irt,int n){
    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    int wnn=3;
    wnn=qpow(wnn,(p-1)>>level,p);
    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = 1LL*rt[i - 1] * rt[i - 1] % p;
    irt[0] = qpow(wnn,p-2,p);
    for (int i = 1; i < level; i++) irt[i] = 1LL*irt[i - 1] * irt[i - 1] % p;
    ntt_dit_x4(a,rt,irt,limit,level);
    ntt_dit_x4(b,rt,irt,limit,level);
	for(int i = 0;i<limit;i++){
        ab[i]=1LL*a[i]*b[i]%p;
    }
	ntt_dif_x4(ab,rt,irt,limit,level);
    int invn=qpow(limit,p-2,p);
    // std::cout<<invn<<'\n';
    for(int i = 0; i < 2 * n; i++){
        ab[i] = (1LL * ab[i] * invn) % p;
        // printf("%d\n",ab[i]);
    } 
}
void ntt_dit_x4_openmp(int *a, int *rt, int *irt, int limit, int level) {
    int one = 1, imag = rt[level-2];
    int logn = __builtin_ctz(limit);
    int dw[level - 1];
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL*dw[i - 1] * irt[level - 1 - i] % p * rt[level - 3 - i] % p;
    dw[level - 2] = 1LL*dw[level - 3] * irt[1] % p;
    #pragma omp parallel num_threads(thread_count)
    {
        if (logn & 1) {
            #pragma omp for  
            for (int i = 0; i < limit/2; i++) {
                int x = a[i], y = a[i + limit/2];
                a[i] = (1LL * x + y) % p;
                a[i + limit/2] = ((1LL * x - y) % p + p) % p;
            }
        }

        for (int e = logn & ~1; e >= 2; e -= 2) {
            const int m = 1 << e, m4 = m >> 2;
            int w2 = 1;
            #pragma omp for
            for (int i = 0; i < limit; i += m) {
                const int w1 = 1LL * w2 * w2 % p, w3 = 1LL * w1 * w2 % p;
                for (int j = i; j < i + m4; ++j) {
                    int a0 = 1LL * a[j + m4 * 0] * one % p, a1 = 1LL * a[j + m4 * 1] * w2 % p;
                    int a2 = 1LL * a[j + m4 * 2] * w1 % p, a3 = 1LL * a[j + m4 * 3] * w3 % p;
                    int t02p = (1LL * a0 + a2) % p, t13p = (1LL * a1 + a3) % p;
                    int t02m = ((1LL * a0 - a2) % p + p) % p, t13m = ((1LL * a1 - a3) % p + p) % p * imag % p;
                    a[j + m4 * 0] = (1LL * t02p + t13p) % p;
                    a[j + m4 * 1] = ((1LL * t02p - t13p) % p + p) % p;
                    a[j + m4 * 2] = (1LL * t02m + t13m) % p;
                    a[j + m4 * 3] = ((1LL * t02m - t13m) % p + p) % p;
                }
                w2 = 1LL * w2 * dw[__builtin_ctz(~(i >> e))] % p;
            }
        }
    }
}
void ntt_dif_x4_openmp(int *a, int *rt, int *irt, int limit, int level) {
    int one = 1, imag = irt[level-2];
    int logn = __builtin_ctz(limit);
    int dw[level - 1];
    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL * dw[i - 1] * rt[level - 1 - i] % p * irt[level - 3 - i] % p;
    dw[level - 2] = 1LL * dw[level - 3] * rt[1] % p;

    #pragma omp parallel num_threads(thread_count)
    {
        for (int e = 2; e <= logn; e += 2) {
            const int m = 1 << e, m4 = m >> 2;
            int w2 = one;
            #pragma omp for
            for (int i = 0; i < limit; i += m) {
                const int w1 = 1LL * w2 * w2 % p, w3 = 1LL * w1 * w2 % p;
                for (int j = i; j < i + m4; ++j) {
                    int a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                    int a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                    int t01p = (1LL * a0 + a1) % p, t23p = (1LL * a2 + a3) % p;
                    int t01m = ((1LL * a0 - a1) % p + p) % p, t23m = ((1LL * a2 - a3) % p + p) % p * imag % p;
                    a[j + m4 * 0] = (1LL * t01p + t23p) % p * one % p;
                    a[j + m4 * 2] = ((1LL * t01p - t23p) % p + p) % p * w1 % p;
                    a[j + m4 * 1] = (1LL * t01m + t23m) % p * w2 % p;
                    a[j + m4 * 3] = ((1LL * t01m - t23m) % p + p) % p * w3 % p;
                }
                w2 = 1LL * w2 * dw[__builtin_ctz(~(i >> e))] % p;
            }
        }
        if (logn & 1) {
            #pragma omp for
            for (int i = 0; i < limit/2; i++) {
                int x = a[i], y = a[i + limit/2];
                a[i] = (1LL * x + y) % p;
                a[i + limit/2] = ((1LL * x - y) % p + p) % p;
            }
        }
    }
}
void ntt_dif_x4_openmp(int *a,int *b,int *ab,int *rt,int *irt,int n){
    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    int wnn=3;
    wnn=qpow(wnn,(p-1)>>level,p);
    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = 1LL*rt[i - 1] * rt[i - 1] % p;
    irt[0] = qpow(wnn,p-2,p);
    for (int i = 1; i < level; i++) irt[i] = 1LL*irt[i - 1] * irt[i - 1] % p;
    ntt_dit_x4_openmp(a,rt,irt,limit,level);
    ntt_dit_x4_openmp(b,rt,irt,limit,level);
    #pragma omp parallel for num_threads(thread_count) schedule(static)
	for(int i = 0;i<limit;i++){
        ab[i]=1LL*a[i]*b[i]%p;
    }
	ntt_dif_x4_openmp(ab,rt,irt,limit,level);
    int invn=qpow(limit,p-2,p);
    // std::cout<<invn<<'\n';
    #pragma omp parallel for num_threads(thread_count) schedule(static)
    for(int i = 0; i < 2 * n; i++){
        ab[i] = (1LL * ab[i] * invn) % p;
        // printf("%d\n",ab[i]);
    } 
}

void ntt_dit_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level){
    int one = 1;
    int logn = __builtin_ctz(limit);
    if (logn & 1) {
        for (int i = 0; i < limit/2; i++) {
            Mint x = a[i], y = a[i + limit/2];
            a[i] = x + y;
            a[i + limit/2] = x - y;
        }
    }
    for (int e = logn & ~1; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = 1;
        for (int i = 0; i < limit; i += m) {
            Mint w1 = w2 * w2, w3 = w1 * w2;
            for (int j = i; j < i + m4; ++j) {
                Mint a0 = a[j + m4 * 0] * one, a1 = a[j + m4 * 1] * w2;
                Mint a2 = a[j + m4 * 2] * w1, a3 = a[j + m4 * 3] * w3;
                Mint t02p = a0 + a2, t13p = a1 + a3;
                Mint t02m = a0 - a2, t13m = (a1 - a3) * imag;
                a[j + m4 * 0] = t02p + t13p;
                a[j + m4 * 1] = t02p - t13p;
                a[j + m4 * 2] = t02m + t13m;
                a[j + m4 * 3] = t02m - t13m;
            }
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))];
            // printf("%d ",__builtin_ctz(~(i >> e)));
        }
        // puts("");
    }
}
void ntt_dif_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level){
    int one = 1;
    int logn = __builtin_ctz(limit);
    
    for (int e = 2; e <= logn; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = one;
        for (int i = 0; i < limit; i += m) {
            Mint w1 = w2 * w2, w3 = w1 * w2;
            for (int j = i; j < i + m4; ++j) {
                Mint a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                Mint a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                Mint t01p = a0 + a1, t23p = a2 + a3;
                Mint t01m = a0 - a1, t23m = (a2 - a3)* imag;
                a[j + m4 * 0] = (t01p + t23p)*one;
                a[j + m4 * 2] = (t01p - t23p)*w1;
                a[j + m4 * 1] = (t01m + t23m)*w2;
                a[j + m4 * 3] = (t01m - t23m)*w3;
            }
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        // printf("%d ",w2);
    }
    if (logn & 1) {
        for (int i = 0; i < limit/2; i++) {
            Mint x = a[i], y = a[i + limit/2];
            a[i] = x + y ;
            a[i + limit/2] = x - y;
        }
    }
}
void ntt_dif_x4_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n){
    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    Mint wnn=3;
    wnn=wnn.pow((p-1)>>level);
    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = rt[i - 1] * rt[i - 1];
    irt[0] = 1/wnn;
    for (int i = 1; i < level; i++) irt[i] = irt[i - 1] * irt[i - 1];
    Mint dw[level - 1];
    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * rt[level - 1 - i] * irt[level - 3 - i];
    dw[level - 2] = dw[level - 3] * rt[1];

    ntt_dit_x4_Mint(a,dw,irt[level-2],limit,level);
    ntt_dit_x4_Mint(b,dw,irt[level-2],limit,level);

	for(int i = 0;i<limit;i++){
        ab[i]=a[i]*b[i];
    }
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * irt[level - 1 - i] * rt[level - 3 - i];//,printf("%d ",dw[i]);
    dw[level - 2] = 1LL*dw[level - 3] * irt[1];
	ntt_dif_x4_Mint(ab,dw,rt[level-2],limit,level);
    Mint Limit=limit;
    Mint invn=1/Limit;
    for(int i = 0; i < 2 * n; i++){
        ab[i] = ab[i] * invn;
    } 
}
void ntt_dit_x4_Mint_openmp(Mint *a,Mint *dw,Mint imag,int limit,int level){
    int one = 1;
    int logn = __builtin_ctz(limit);
    #pragma omp parallel num_threads(thread_count)
    {
        if (logn & 1) {
            #pragma omp for
            for (int i = 0; i < limit/2; i++) {
                Mint x = a[i], y = a[i + limit/2];
                a[i] = x + y;
                a[i + limit/2] = x - y;
            }
        }
        for (int e = logn & ~1; e >= 2; e -= 2) {
            const int m = 1 << e, m4 = m >> 2;
            Mint w2 = 1;
            #pragma omp for
            for (int i = 0; i < limit; i += m) {
                Mint w1 = w2 * w2, w3 = w1 * w2;
                for (int j = i; j < i + m4; ++j) {
                    Mint a0 = a[j + m4 * 0] * one, a1 = a[j + m4 * 1] * w2;
                    Mint a2 = a[j + m4 * 2] * w1, a3 = a[j + m4 * 3] * w3;
                    Mint t02p = a0 + a2, t13p = a1 + a3;
                    Mint t02m = a0 - a2, t13m = (a1 - a3) * imag;
                    a[j + m4 * 0] = t02p + t13p;
                    a[j + m4 * 1] = t02p - t13p;
                    a[j + m4 * 2] = t02m + t13m;
                    a[j + m4 * 3] = t02m - t13m;
                }
                w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))];
                // printf("%d ",__builtin_ctz(~(i >> e)));
            }
            // puts("");
        }
    }
}
void ntt_dif_x4_Mint_openmp(Mint *a,Mint *dw,Mint imag,int limit,int level){
    int one = 1;
    int logn = __builtin_ctz(limit);
    #pragma omp parallel num_threads(thread_count)
    {
        for (int e = 2; e <= logn; e += 2) {
            const int m = 1 << e, m4 = m >> 2;
            Mint w2 = one;
            #pragma omp for
            for (int i = 0; i < limit; i += m) {
                Mint w1 = w2 * w2, w3 = w1 * w2;
                for (int j = i; j < i + m4; ++j) {
                    Mint a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                    Mint a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                    Mint t01p = a0 + a1, t23p = a2 + a3;
                    Mint t01m = a0 - a1, t23m = (a2 - a3)* imag;
                    a[j + m4 * 0] = (t01p + t23p)*one;
                    a[j + m4 * 2] = (t01p - t23p)*w1;
                    a[j + m4 * 1] = (t01m + t23m)*w2;
                    a[j + m4 * 3] = (t01m - t23m)*w3;
                }
                w2 = w2*dw[__builtin_ctz(~(i >> e))];
            }
            // printf("%d ",w2);
        }
        if (logn & 1) {
            #pragma omp for
            for (int i = 0; i < limit/2; i++) {
                Mint x = a[i], y = a[i + limit/2];
                a[i] = x + y ;
                a[i + limit/2] = x - y;
            }
        }
    }
    

}
void ntt_dif_x4_Mint_openmp(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n){
    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    Mint wnn=3;
    wnn=wnn.pow((p-1)>>level);
    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = rt[i - 1] * rt[i - 1];
    irt[0] = 1/wnn;
    for (int i = 1; i < level; i++) irt[i] = irt[i - 1] * irt[i - 1];
    Mint dw[level - 1];
    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * rt[level - 1 - i] * irt[level - 3 - i];
    dw[level - 2] = dw[level - 3] * rt[1];

    ntt_dit_x4_Mint_openmp(a,dw,irt[level-2],limit,level);
    ntt_dit_x4_Mint_openmp(b,dw,irt[level-2],limit,level);

	for(int i = 0;i<limit;i++){
        ab[i]=a[i]*b[i];
    }
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * irt[level - 1 - i] * rt[level - 3 - i];//,printf("%d ",dw[i]);
    dw[level - 2] = 1LL*dw[level - 3] * irt[1];
	ntt_dif_x4_Mint_openmp(ab,dw,rt[level-2],limit,level);
    Mint Limit=limit;
    Mint invn=1/Limit;
    for(int i = 0; i < 2 * n; i++){
        ab[i] = ab[i] * invn;
    } 
}