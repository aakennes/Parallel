#include "ntt.h"
#include "mod_arith.h"
#include "permutation.h"
#include <cmath>
#include <map>

int qpow(int a,int x,int p){
    int ans=1;
    while(x){
        if(x&1)ans=1LL*ans*a%p;
        a=1LL*a*a%p;
        x>>=1;
    }
    return ans%p;
}

#define i32 int32_t
#define u32 uint32_t
#define i64 int64_t
#define u64 uint64_t

template <std::uint32_t P> 
struct MontgomeryModInt32 {
public:
    u32 v;
private:
    const static u32 Montr = 104857599;
    const static u32 Montr2 = 45971250; 
public:
    static  u32 pow_mod(u32 x, u64 y) {
        if ((y %= P - 1) < 0)
            y += P - 1;

        u32 res = 1;

        for (; y != 0; y >>= 1, x = u64(x) * x % P)
            if (y & 1)
                res = u64(res) * x % P;

        return res;
    }
 
    MontgomeryModInt32() = default;
    ~MontgomeryModInt32() = default;
    MontgomeryModInt32(u32 v) : v(reduce(u64(v) * Montr2)) {}
    MontgomeryModInt32(const MontgomeryModInt32 &rhs) : v(rhs.v) {}
    static  u32 reduce(u64 x) {
        return x + (u64(u32(x) * Montr) * P) >> 32;
    }
     u32 get() const {
        u32 res = reduce(v);
        return res - (P & -(res >= P));
    }
    explicit  operator u32() const {
        return get();
    }
    explicit  operator i32() const {
        return i32(get());
    }
    MontgomeryModInt32 &operator=(const MontgomeryModInt32 &rhs) {
        return v = rhs.v, *this;
    }
    MontgomeryModInt32 operator-() const {
        MontgomeryModInt32 res;
        return res.v = (P << 1 & -(v != 0)) - v, res;
    }
    MontgomeryModInt32 inv() const {
        return pow(-1);
    }
    MontgomeryModInt32 &operator+=(const MontgomeryModInt32 &rhs) {
        return v += rhs.v - (P << 1), v += P << 1 & -(i32(v) < 0), *this;
    }
    MontgomeryModInt32 &operator-=(const MontgomeryModInt32 &rhs) {
        return v -= rhs.v, v += P << 1 & -(i32(v) < 0), *this;
    }
    MontgomeryModInt32 &operator*=(const MontgomeryModInt32 &rhs) {
        return v = reduce(u64(v) * rhs.v), *this;
    }
    MontgomeryModInt32 &operator/=(const MontgomeryModInt32 &rhs) {
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
        return is >> rhs.v, rhs.v = reduce(u64(rhs.v) * Montr2), is;
    }
    friend std::ostream &operator<<(std::ostream &os, const MontgomeryModInt32 &rhs) {
        return os << rhs.get();
    }
     MontgomeryModInt32 pow(i64 y) const {
        if ((y %= P - 1) < 0)
            y += P - 1; // phi(P) = P - 1, assume P is a prime number

        MontgomeryModInt32 res(1), x(*this);

        for (; y != 0; y >>= 1, x *= x)
            if (y & 1)
                res *= x;

        return res;
    }
};

namespace hehub {

//common
u64 __pow_mod(u64 modulus, u64 base, size_t index) {
    u64 power = 1;
    size_t mask = 1;
    while (mask <= index) {
        mask <<= 1;
    }
    mask >>= 1;
    while (mask) {
        power = (u128)power * power % modulus;
        if (mask & index) {
            power = (u128)power * base % modulus;
        }
        mask >>= 1;
    }
    return power;
}

u64 __get_2nth_unity_root(u64 modulus, u64 n) {
    if ((modulus - 1) % (2 * n) != 0) {
        throw std::invalid_argument("2N doesn't divide (modulus - 1)");
    }
    u64 candidate = 2;
    while (true) {
        if (__pow_mod(modulus, candidate, (modulus - 1) / 2) == modulus - 1) {
            break;
        }
        candidate++;
    }
    auto root = __pow_mod(modulus, candidate, (modulus - 1) / (2 * n));
    return root;
}

u64 __get_nth_unity_root(u64 modulus, u64 n) {
    if ((modulus - 1) % n != 0) {
        throw std::invalid_argument("N doesn't divide (modulus - 1)");
    }
    u64 candidate = 2;
    while (true) {
        if (__pow_mod(modulus, candidate, (modulus - 1)/2) == modulus - 1) {
            break;
        }
        candidate++;
    }
    return candidate;
}

struct NTTFactors {
    NTTFactors(u64 modulus, size_t log_dimension, bool for_inverse = false) {
        const size_t log_modulus = (u64)(log2(modulus) + 0.5);
        if (log_modulus > 59) {
            throw std::invalid_argument(
                "NTT not supporting primes with bit size > 59 currently.");
        }
        size_t dimension = 1 << log_dimension;
        shuffled_indices.resize(dimension,0);
        for (size_t i = 0; i < dimension; i++) {
            shuffled_indices[i] = (shuffled_indices[i >> 1] >> 1) | ((i & 1) << (log_dimension - 1));
        }
        // printf("%lu\n",shuffled_indices[1]);
    }

    NTTFactors(const NTTFactors &copying) = default;

    NTTFactors(NTTFactors &&moving) = default;

    ~NTTFactors() {}

    std::vector<size_t> shuffled_indices;
};

std::map<std::pair<u64, u64>, NTTFactors> &ntt_factors_cache() {
    static std::map<std::pair<u64, u64>, NTTFactors> global_ntt_factors_cache;
    return global_ntt_factors_cache;
}

std::map<std::pair<u64, u64>, NTTFactors> &intt_factors_cache() {
    static std::map<std::pair<u64, u64>, NTTFactors> global_intt_factors_cache;
    return global_intt_factors_cache;
}

inline const auto &
__find_or_create_ntt_factors(const u64 modulus, const size_t log_dimension,
                             const bool for_inverse = false) {
    if (!for_inverse) {
        auto it =
            ntt_factors_cache().find(std::make_pair(modulus, log_dimension));
        if (it == ntt_factors_cache().end()) {
            ntt_factors_cache().insert(
                std::make_pair(std::make_pair(modulus, log_dimension),
                               NTTFactors(modulus, log_dimension)));
            it =
                ntt_factors_cache().find(std::make_pair(modulus, log_dimension));
        }
        return it->second;
    } else {
        auto it =
            intt_factors_cache().find(std::make_pair(modulus, log_dimension));
        if (it == intt_factors_cache().end()) {
            intt_factors_cache().insert(
                std::make_pair(std::make_pair(modulus, log_dimension),
                               NTTFactors(modulus, log_dimension, true)));
            it = intt_factors_cache().find(
                std::make_pair(modulus, log_dimension));
        }
        return it->second;
    }
}

void ntt_negacyclic_inplace_lazy(const size_t log_dimension, const u64 modulus,
                                 u64 coeffs[]) {
    const size_t dimension = 1ULL << log_dimension;
    // generate or read from cache
    const auto &ntt_factors =
        __find_or_create_ntt_factors(modulus, log_dimension);
    const u64 root_of_2nth = __get_nth_unity_root(modulus, dimension);
    auto root_of_2nth_inv =
                __pow_mod(modulus, root_of_2nth, modulus-2);
    for(int i = 0; i < dimension; i++) {
		if(i < ntt_factors.shuffled_indices[i]){
            std::swap(coeffs[i], coeffs[ntt_factors.shuffled_indices[i]]);
            
        }
    }
    for(int mid = 1; mid < dimension; mid <<= 1) {
        auto Wn = __pow_mod(modulus, root_of_2nth, (modulus - 1) / (mid << 1));
		for(int j = 0; j < dimension; j += (mid << 1)) {
			u64 w = 1;
			for(int k = 0; k < mid; k++, w = ((u128)w * Wn) % modulus) {
                auto x = coeffs[j + k];
                auto y = (u128)w * coeffs[j + k + mid] % modulus;             
                coeffs[j + k] = ((u128)x + y) % modulus,
                coeffs[j + k + mid] = ((u128)x - y + modulus) % modulus;
			}
		}
	}
}

void intt_negacyclic_inplace_lazy(const size_t log_dimension, const u64 modulus,
                                  u64 values[]) {
    const size_t dimension = 1ULL << log_dimension;
    // printf("%lu\n",modulus);
    // generate or read from cache
    const auto &intt_factors =
        __find_or_create_ntt_factors(modulus, log_dimension, true);
    const u64 root_of_2nth = __get_nth_unity_root(modulus, dimension);
    // printf("");
    auto root_of_2nth_inv =
                __pow_mod(modulus, root_of_2nth, modulus-2);
    for(int i = 0; i < dimension; i++) {
		if(i < intt_factors.shuffled_indices[i]){
            std::swap(values[i], values[intt_factors.shuffled_indices[i]]);
        }
    }
    for(int mid = 1; mid < dimension; mid <<= 1) {
        auto Wn = __pow_mod(modulus, root_of_2nth_inv, (modulus - 1) / (mid << 1));
		for(int j = 0; j < dimension; j += (mid << 1)) {
			u64 w = 1;
			for(int k = 0; k < mid; k++, w = ((u128)w * Wn) % modulus) {
                auto x = values[j + k];
                auto y = (u128)w * values[j + k + mid] % modulus;             
                values[j + k] = ((u128)x + y) % modulus,
                values[j + k + mid] = ((u128)x - y + modulus) % modulus;
			}
		}
	}
    auto invn=__pow_mod(modulus,dimension,modulus-2);
    
    // std::cout<<invn<<'\n';
    for(size_t i = 0; i < dimension; i++){
        values[i] = ((u128)values[i] * invn) % modulus;
    } 
}

void cache_ntt_factors_strict(const u64 log_dimension,
                              const std::vector<u64> &moduli) {
    for (auto modulus : moduli) {
        __find_or_create_ntt_factors(modulus, log_dimension);
        __find_or_create_ntt_factors(modulus, log_dimension, true);
    }
}

} // namespace hehub
