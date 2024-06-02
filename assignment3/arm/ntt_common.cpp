#include <pthread.h>
#include <algorithm>
#include <iostream>
#include <semaphore.h>
#include<sys/time.h>
#include<chrono>

const int p = 998244353; // 假设一个模数p
const int numThreads = 4;

int qpow(int base, int exp, int mod) {
    int result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (1LL * result * base) % mod;
        }
        base = (1LL * base * base) % mod;
        exp /= 2;
    }
    return result;
}

struct ButterflyData {
    int *a;
    int *r;
    int wnn;
    int invwnn;
    int id;
    int limit;
    int type;
};

pthread_barrier_t barr_merge;
pthread_barrier_t barr_ntt;
int W[200005];

void* butterflyFunc(void* arg) {
    ButterflyData* data = static_cast<ButterflyData*>(arg);
    int* a = data->a;
    int* r = data->r;
    int wnn = data->wnn;
    int invwnn = data->invwnn;
    int id = data->id;
    int limit = data->limit;
    int type=data->type;

    for (int mid = 1; mid < limit; mid <<= 1) {
        if(id == 0){
            int Wn = qpow(type == 1 ? wnn : invwnn, (p - 1) / (mid << 1), p);
            W[0] = 1;
            for(int k=1;k<mid;++k){
                W[k] = 1LL*W[k-1]*Wn%p;
            }
        }
        pthread_barrier_wait(&barr_ntt);
        for (int j = 0; j < limit; j += (mid << 1)) {
            for (int k = id; k < mid; k+=numThreads) {
                int x = a[j + k], y = 1LL * W[k] * a[j + k + mid] % p;
                a[j + k] = (1LL * x + y) % p;
                a[j + k + mid] = (1LL * x - y + p) % p;
            }
        }
        pthread_barrier_wait(&barr_merge);
    }
    
    return nullptr;
}

void ntt_common(int *a, int *r, int limit, int type) {
    pthread_t threads[numThreads];
    ButterflyData threadData[numThreads];
    
    for(int i = 0; i < limit; i++) {
        if(i < r[i]){
            std::swap(a[i], a[r[i]]);
        }
    }

    int wnn = 3;
    int invwnn = qpow(wnn, p - 2, p);

    pthread_barrier_init(&barr_merge,NULL,numThreads);
    pthread_barrier_init(&barr_ntt,NULL,numThreads);

    for (int i = 0; i < numThreads; ++i) {
        threadData[i] = {a, r, wnn, invwnn, i, limit, type};
        pthread_create(&threads[i], nullptr, butterflyFunc, &threadData[i]);
    }

    for (int i = 0; i < numThreads; ++i) {
        pthread_join(threads[i], nullptr);
    }
    pthread_barrier_destroy(&barr_merge);
    pthread_barrier_destroy(&barr_ntt);
}

void ntt_common(int *a, int *b, int *ab, int *r, int n) {
    
    int L = 0, i = 0;
    int limit = 1;
    while (limit <= 2 * n - 1) {
        limit <<= 1;
        L++;
    }
    for (i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    
    ntt_common(a, r, limit, 1);
    
    ntt_common(b, r, limit, 1);

    for (i = 0; i < limit; i++) {
        if(i<10){
            std::cout<<a[i]<<" "<<b[i]<<'\n';
        }
        ab[i] = 1LL * a[i] * b[i] % p;
    }
    

    ntt_common(ab, r, limit, -1);
    int invn = qpow(limit, p - 2, p);

    for (i = 0; i < 2 * n; i++) {
        ab[i] = (1LL * ab[i] * invn) % p;
    }
}

const int maxn=2e5+50;
int n,m;
int a[maxn],b[maxn],ab[maxn],r[maxn],ir[maxn];

int nn[11]={256,512,1024,2048,4096,8192,16384,32768,65536,131072};
int qq[5]={1409,3329,7681,12289};
int main() {
    // freopen("2.in","r",stdin);
	for(int i=0;i<=9;++i){
        for(int j=0;j<=3;++j){
            long double ans=0;
            int cnt=100;
            for(int k=1;k<=100;++k){
                FILE* fp;
                // fRead(fp,a,b,nn[i],qq[j]);
                // fRead(fp,A,B,nn[i],qq[j]);
                // std::cout<<nn[i]<<" "<<qq[j]<<'\n';
                auto Start=std::chrono::high_resolution_clock::now();
                ntt_common(a,b,ab,r,nn[i]);
                // ntt_Montgomery(a,b,ab,r,nn[i]);
                // ntt_Montgomery_Mint(A,B,AB,r,nn[i]);
                // ntt_dif(a,b,ab,rt,irt,nn[i]);
                // ntt_dif_Mint(A,B,AB,RT,IRT,nn[i]);
                // ntt_dif_x4(a,b,ab,rt,irt,nn[i]);
                // ntt_dif_x4_Mint(A,B,AB,RT,IRT,nn[i]);
                auto End=std::chrono::high_resolution_clock::now();
                std::chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
                ans+=elapsed.count();
                // std::cout<<elapsed.count()<<std::endl;
            }   
            std::cout<<ans/cnt<<std::endl;
        }
    }
    int limit=1,L=0;

	for(int i=0;i<n;++i)std::cin>>a[i];
	for(int i=0;i<m;++i)std::cin>>b[i];//,std::cout<<b[i]<<'\n';
    
    ntt_common(a, b, ab, r, std::max(n,m));

    for (int i = 0; i < 5; i++) {
        std::cout << ab[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
