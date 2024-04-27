#include<cstdio>
#include<iostream>
#include<immintrin.h>
#include<stdint.h>
#define ll long long
const ll MOD = 104857601,p=MOD;
const ll G = 3;
const ll MAXN = 1 << 24 | 1;
#define ulong long long 
#define level 22

ll rt[MAXN], irt[MAXN];

ll add(ll x, ll v) { return x + v >= MOD ? x + v - MOD : x + v; }
ll dec(ll x, ll v) { return x - v < 0 ? x - v + MOD : x - v; }

ll modPow( ll a,  ll b) {
     ll ret = 1;
    for (; b; b >>= 1, a = (ulong)a * a % MOD)
        if (b & 1) ret = (ulong)ret * a % MOD;
    return ret;
}

__m256i mod_add_epi64(__m256i a, __m256i b) {
    __m256i result = _mm256_add_epi64(a, b);
    __m256i mod = _mm256_set1_epi64x(MOD);
    __m256i cmp = _mm256_cmpgt_epi64(result, mod);
    return _mm256_sub_epi64(result, _mm256_and_si256(cmp, mod));
}

// // 使用AVX2执行向量化的模减操作
// __m256i mod_sub_epi64(__m256i a, __m256i b) {
//     __m256i result = _mm256_sub_epi64(a, b);
//     __m256i mod = _mm256_set1_epi64x(MOD);
//     __m256i cmp = _mm256_cmpgt_epi64(_mm256_setzero_si256(), result);
//     return _mm256_add_epi64(result, _mm256_and_si256(cmp, mod));
// }

// __m256i scalar_mod_mul_epi64(__m256i a, __m256i b) {
//     __m256i result = _mm256_setzero_si256();

//     for (int i = 0; i < 4; i++) {
//         long long a_i = _mm256_extract_epi64(a, i);
//         long long b_i = _mm256_extract_epi64(b, i);
//         long long r_i = a_i * b_i % MOD;
//         result = _mm256_insert_epi64(result, r_i, i);
//     }

//     return result;
// }

 void dit(ll *a, ll n) {
    ll one = 1, imag = rt[level-2];
    ll logn = __builtin_ctz(n);
    // printf("%d\n",n);
    ll dw[level - 1];
    dw[0] = rt[level - 3];
    for (ll i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * irt[level - 1 - i]%MOD * rt[level - 3 - i]%MOD;
    dw[level - 2] = dw[level - 3] * irt[1]%MOD;
    bool asd=0;
    if (logn & 1) {
        for (ll i = 0; i < n/2; i+=4) {
            // __m256i x = _mm256_loadu_si256((__m256i *)(a + i));
            // __m256i y = _mm256_loadu_si256((__m256i *)(a + i + n / 2));
            // __m256i sum = mod_add_epi64(x, y);
            // __m256i diff = mod_sub_epi64(x, y);

            // _mm256_storeu_si256((__m256i *)(a + i), sum);
            // _mm256_storeu_si256((__m256i *)(a + i + n / 2), diff);
            // if(asd==0&&a[i]==0)asd=1,printf("%lld %lld %lld %lld\n",i + 2*n,a[i+2*n],x,y);
        }
    }
    // for (ll e = logn & ~1; e >= 2; e -= 2) {
    //     const ll m = 1 << e, m4 = m >> 2;
    //     ll w2 = 1;
        
    //     for (ll i = 0; i < n; i += m) {
    //         const ll w1 = w2 * w2 %MOD, w3 = w1 * w2%MOD;
    //         // printf("%lld %lld\n",w2,w3);
    //         for (ll j = i; j < i + m4; j+=4) {
    //             __m256i a0 = scalar_mod_mul_epi64(_mm256_loadu_si256((__m256i *)(a + j + m4 * 0)), _mm256_set1_epi64x(one));
    //             __m256i a1 = scalar_mod_mul_epi64(_mm256_loadu_si256((__m256i *)(a + j + m4 * 1)), _mm256_set1_epi64x(w2));
    //             __m256i a2 = scalar_mod_mul_epi64(_mm256_loadu_si256((__m256i *)(a + j + m4 * 2)), _mm256_set1_epi64x(w1));
    //             __m256i a3 = scalar_mod_mul_epi64(_mm256_loadu_si256((__m256i *)(a + j + m4 * 3)), _mm256_set1_epi64x(w3));

    //             __m256i t02p = mod_add_epi64(a0,a2);
    //             __m256i t02m = mod_sub_epi64(a0,a2);
    //             __m256i t13p = mod_add_epi64(a1,a3);
    //             __m256i t13m = mod_sub_epi64(a1,a3);
    //             t13m = scalar_mod_mul_epi64(t13m,_mm256_set1_epi64x(imag));
    //             _mm256_storeu_si256((__m256i *)(a + j + m4*0), mod_add_epi64(t02p,t13p));
    //             _mm256_storeu_si256((__m256i *)(a + j + m4*1), mod_sub_epi64(t02p,t13p));
    //             _mm256_storeu_si256((__m256i *)(a + j + m4*2), mod_add_epi64(t02m,t13m));
    //             _mm256_storeu_si256((__m256i *)(a + j + m4*3), mod_sub_epi64(t02m,t13m));
    //             // ll a0 = a[j + m4 * 0] * one%MOD, a1 = a[j + m4 * 1] * w2%MOD;
    //             // ll a2 = a[j + m4 * 2] * w1%MOD, a3 = a[j + m4 * 3] * w3%MOD;
    //             // ll t02p = (a0 + a2)%MOD, t13p = (a1 + a3)%MOD;
    //             // ll t02m = ((a0 - a2 )%MOD+MOD)%MOD, t13m = ((a1 - a3 )%MOD+MOD)%MOD * imag%MOD;
    //             // a[j + m4 * 0] = (t02p + t13p)%MOD;
    //             // a[j + m4 * 1] = ((t02p - t13p)%MOD+MOD)%MOD;
    //             // a[j + m4 * 2] = (t02m + t13m)%MOD;
    //             // a[j + m4 * 3] = ((t02m - t13m)%MOD+MOD)%MOD;
    //         }
    //         w2 *= dw[__builtin_ctz(~(i >> e))];
    //         w2%=MOD;
    //     }
    //     // printf("%d ",w2);
    // }
    // for(ll i=0;i<2;++i)printf("%d ",a[i]);
}

 void dif(ll *a, ll n) {
    ll one = 1, imag = irt[level-2];
    ll logn = __builtin_ctz(n);
    // printf("%d\n",one);
    ll dw[level - 1];
    dw[0] = irt[level - 3];
    for (ll i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * rt[level - 1 - i]%MOD * irt[level - 3 - i]%MOD;
    dw[level - 2] = dw[level - 3] * rt[1]%MOD;
    
    for (ll e = 2; e <= logn; e += 2) {
        const ll m = 1 << e, m4 = m >> 2;
        ll w2 = one;
        for (ll i = 0; i < n; i += m) {
            const ll w1 = w2 * w2 %MOD, w3 = w1 * w2 %MOD;
            for (ll j = i; j < i + m4; ++j) {
                ll a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                ll a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                ll t01p = (a0 + a1)%MOD, t23p = (a2 + a3)%MOD;
                ll t01m = ((a0 - a1)%MOD + MOD)%MOD, t23m = ((a2 - a3)%MOD +MOD)%MOD * imag%MOD;
                a[j + m4 * 0] = (t01p + t23p)%MOD*one%MOD;
                a[j + m4 * 2] = ((t01p - t23p)%MOD+MOD)%MOD*w1%MOD;
                a[j + m4 * 1] = (t01m + t23m)%MOD*w2%MOD;
                a[j + m4 * 3] = ((t01m - t23m)%MOD+MOD)%MOD*w3%MOD;
            }
            w2 *= dw[__builtin_ctz(~(i >> e))];
            w2%=MOD;
        }
        // printf("%lld ",w2);
    }
    if (logn & 1) {
        for (ll i = 0; i < n/2; i++) {
            ll x = a[i], y = a[i + n/2];
            a[i] = (x + y )% MOD;
            a[i + n/2] = ((x - y)%MOD + MOD)%MOD;
        }
    }
}
ll n,m;
ll a[MAXN],b[MAXN],ab[MAXN],r[MAXN],ir[MAXN];
int main(){
	freopen("2.in","r",stdin);
	std::cin>>n>>m;
	n++,m++;
    ll limit=1,L=0;
    while(limit <= m + n-2){
        limit <<= 1, L++;
    }
    ll wnn=3;
    wnn=modPow(wnn,(MOD-1)>>level);
	rt[0] = wnn;
    // printf("%d\n",wnn);
    for (ll i = 1; i < level; i++) rt[i] = (ulong)rt[i - 1] * rt[i - 1] % MOD;
    irt[0] = modPow(wnn,MOD-2);
    for (ll i = 1; i < level; i++) irt[i] = (ulong)irt[i - 1] * irt[i - 1] % MOD;//,printf("%d ",irt[i]);
	for(ll i=0;i<n;++i)std::cin>>a[i];
	for(ll i=0;i<m;++i)std::cin>>b[i];
	dit(a,limit);
	dit(b,limit);
	for(ll i = 0;i<limit;i++){
        ab[i]=1LL*a[i]*b[i]%p;
    }
	dif(ab,limit);
    ll invn=modPow(limit,p-2);
    for(ll i = 0; i < n + m -1 ; i++){
        ab[i] = 1LL * ab[i] * invn % p;
        printf("%lld ",ab[i]);
    } 
}