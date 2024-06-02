#include <pthread.h>
#include <algorithm>
#include <iostream>
#include <semaphore.h>

#define Original_Q 104857601
#define Mint MontgomeryModInt32<Original_Q>

#define i32 int32_t
#define u32 uint32_t
#define i64 int64_t
#define u64 uint64_t
#define MMint u32

#define numThreads 4

const int p = Original_Q; // 假设一个模数p

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

template <std::uint32_t P> 
struct MontgomeryModInt32 {

private:
    u32 v;

    static constexpr u32 get_r() {
        u32 iv = P;

        for (u32 i = 0; i != 4; ++i)
            iv *= 2 - P * iv;
        // std::cout<<-iv<<" "<<-u64(P) % P<<'\n';
        return iv;
    }

    static constexpr u32 r = -get_r(), r2 = -u64(P) % P;

public:
    static constexpr u32 pow_mod(u32 x, u64 y) {
        if ((y %= P - 1) < 0)
            y += P - 1;

        u32 res = 1;

        for (; y != 0; y >>= 1, x = u64(x) * x % P)
            if (y & 1)
                res = u64(res) * x % P;

        return res;
    }

    static constexpr u32 get_pr() {
        u32 tmp[32] = {}, cnt = 0;
        const u64 phi = P - 1;
        u64 m = phi;

        for (u64 i = 2; i * i <= m; ++i) {
            if (m % i == 0) {
                tmp[cnt++] = i;

                while (m % i == 0)
                    m /= i;
            }
        }

        if (m > 1)
            tmp[cnt++] = m;

        for (u64 res = 2; res <= phi; ++res) {
            bool flag = true;

            for (u32 i = 0; i != cnt && flag; ++i)
                flag &= pow_mod(res, phi / tmp[i]) != 1;

            if (flag)
                return res;
        }

        return 0;
    }

    MontgomeryModInt32() = default;
    ~MontgomeryModInt32() = default;
    constexpr MontgomeryModInt32(u32 v) : v(reduce(u64(v) * r2)) {}
    constexpr MontgomeryModInt32(const MontgomeryModInt32 &rhs) : v(rhs.v) {}
    static constexpr u32 reduce(u64 x) {
        return x + (u64(u32(x) * r) * P) >> 32;
    }
    constexpr u32 getv() const {
        return v;
    }
    constexpr u32 get() const {
        u32 res = reduce(v);
        return res - (P & -(res >= P));
    }
    explicit constexpr operator u32() const {
        return get();
    }
    explicit constexpr operator i32() const {
        return i32(get());
    }
    constexpr MontgomeryModInt32 &operator=(const MontgomeryModInt32 &rhs) {
        return v = rhs.v, *this;
    }
    constexpr MontgomeryModInt32 operator-() const {
        MontgomeryModInt32 res;
        return res.v = (P << 1 & -(v != 0)) - v, res;
    }
    constexpr MontgomeryModInt32 inv() const {
        return pow(-1);
    }
    constexpr MontgomeryModInt32 &operator+=(const MontgomeryModInt32 &rhs) {
        return v += rhs.v - (P << 1), v += P << 1 & -(i32(v) < 0), *this;
    }
    constexpr MontgomeryModInt32 &operator-=(const MontgomeryModInt32 &rhs) {
        return v -= rhs.v, v += P << 1 & -(i32(v) < 0), *this;
    }
    constexpr MontgomeryModInt32 &operator*=(const MontgomeryModInt32 &rhs) {
        return v = reduce(u64(v) * rhs.v), *this;
    }
    constexpr MontgomeryModInt32 &operator/=(const MontgomeryModInt32 &rhs) {
        return this->operator*=(rhs.inv());
    }
    friend MontgomeryModInt32 operator+(const MontgomeryModInt32 &lhs,
                                        const MontgomeryModInt32 &rhs) {
        return MontgomeryModInt32(lhs) += rhs;
    }
    friend MontgomeryModInt32 operator-(const MontgomeryModInt32 &lhs,
                                        const MontgomeryModInt32 &rhs) {
        return MontgomeryModInt32(lhs) -= rhs;
    }
    friend MontgomeryModInt32 operator*(const MontgomeryModInt32 &lhs,
                                        const MontgomeryModInt32 &rhs) {
        return MontgomeryModInt32(lhs) *= rhs;
    }
    friend MontgomeryModInt32 operator/(const MontgomeryModInt32 &lhs,
                                        const MontgomeryModInt32 &rhs) {
        return MontgomeryModInt32(lhs) /= rhs;
    }
    friend std::istream &operator>>(std::istream &is, MontgomeryModInt32 &rhs) {
        return is >> rhs.v, rhs.v = reduce(u64(rhs.v) * r2), is;
    }
    friend std::ostream &operator<<(std::ostream &os, const MontgomeryModInt32 &rhs) {
        return os << rhs.get();
    }
    constexpr MontgomeryModInt32 pow(i64 y) const {
        if ((y %= P - 1) < 0)
            y += P - 1; // phi(P) = P - 1, assume P is a prime number
        // if(y==204800)std::cout<<y<<'\n';
        MontgomeryModInt32 res(1), x(*this);

        for (; y != 0; y >>= 1, x *= x)
            if (y & 1)
                res *= x;
        // if(y==204800)std::cout<<res<<'\n';
        return res;
    }
};


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
int W[100005];

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

const int maxn=2e6+50;
int n,m;
int a[maxn],b[maxn],ab[maxn],r[maxn],ir[maxn];
int main() {
    freopen("1.in","r",stdin);
	std::cin>>n>>m;
	n++,m++;
    int limit=1,L=0;

	for(int i=0;i<n;++i)std::cin>>a[i];
	for(int i=0;i<m;++i)std::cin>>b[i];

    ntt_Montgomery(a, b, ab, r, std::max(n,m));


    for (int i = 0; i < n+m-1; i++) {
        std::cout << ab[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
