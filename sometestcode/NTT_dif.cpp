#include<cstdio>
#include<iostream>

const int MOD = 998244353,p=MOD;
const int G = 3;
const int MAXN = 1 << 24 | 1;

#define ulong long long 

int rt[MAXN], irt[MAXN];

inline int add(int x, int v) { return x + v >= MOD ? x + v - MOD : x + v; }
inline int dec(int x, int v) { return x - v < 0 ? x - v + MOD : x - v; }

inline int modPow(register int a, register int b) {
    register int ret = 1;
    for (; b; b >>= 1, a = (ulong)a * a % MOD)
        if (b & 1) ret = (ulong)ret * a % MOD;
    return ret;
}

inline void dit(int *a, int n) {
    for (int i = 1, mid = n >> 1; i < n; i <<= 1, mid >>= 1) {
        for (int j = 0, w, o = 0; j < i; j++, o += mid << 1) {
            w = rt[i + j];
            for (int k = o, t; k < o + mid; k++) {
                t = (ulong)a[k + mid] * w % MOD;
                a[k + mid] = dec(a[k], t);
                a[k] = add(a[k], t);
            }
        }
    }
}

inline void dif(int *a, int n) {
    for (int i = n >> 1, mid = 1; i; i >>= 1, mid <<= 1) {
        for (int j = 0, w, o = 0; j < i; j++, o += mid << 1) {
            w = irt[i + j];
            for (int k = o, t; k < o + mid; k++) {
                t = a[k + mid];
                a[k + mid] = dec(a[k], t) * (ulong)w % MOD;
                a[k] = add(a[k], t);
            }
        }
    }
}
const int maxn=2e6+50;
int n,m;
int a[maxn],b[maxn],ab[maxn],r[maxn],ir[maxn];
int main(){
	freopen("1.in","r",stdin);
	std::cin>>n>>m;
	n++,m++;
    int limit=1,mid=0;
    while(limit <= m + n-2){
        limit <<= 1, mid++;
    } 
    int wnn=3;
    int invwnn=modPow(wnn,p-2);
	rt[0] = 1;
    rt[1] = modPow(G, (MOD - 1) / limit / 2);
    for (int i = 2; i < limit; i++) rt[i] = (ulong)rt[i - 1] * rt[1] % MOD;
    irt[0] = 1;
    irt[1] = modPow(rt[1], MOD - 2);
    for (int i = 2; i < limit; i++) irt[i] = (ulong)irt[i - 1] * irt[1] % MOD;
    for (register int i = 0, j = 0; i < limit; i++) {
        if (i > j) {
            std::swap(rt[i], rt[j]);
            std::swap(irt[i], irt[j]);
        }
        for (register int t = limit >> 1; (j ^= t) < t; t >>= 1)
            ;
    }
	for(int i=0;i<n;++i)std::cin>>a[i];
	for(int i=0;i<m;++i)std::cin>>b[i];
	dit(a,limit);
	dit(b,limit);
	for(int i = 0;i<limit;i++){
        ab[i]=1LL*a[i]*b[i]%p;
    }
	dif(ab,limit);
    int invn=modPow(limit,p-2);
    for(int i = 0; i < n + m -1 ; i++){
        ab[i] = (1LL * ab[i] * invn) % p;
        printf("%d ",ab[i]);
    } 
}