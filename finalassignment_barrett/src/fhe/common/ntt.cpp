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
private:
    u32 v;

    static constexpr u32 get_r() {
        u32 iv = P;
        for (u32 i = 0; i != 4; ++i)
            iv *= 2 - P * iv;
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
            y += P - 1;
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

struct NTTFactors {
    NTTFactors(u64 modulus, size_t log_dimension, bool for_inverse = false) {
        const size_t log_modulus = (u64)(log2(modulus) + 0.5);
        if (log_modulus > 59) {
            throw std::invalid_argument(
                "NTT not supporting primes with bit size > 59 currently.");
        }
        size_t dimension = 1 << log_dimension;

        const u64 root_of_2nth = __get_2nth_unity_root(modulus, dimension);
        if (!for_inverse) {
            seq.resize(dimension);
            seq_harvey.resize(dimension);
            for (size_t i = 0; i < dimension; i++) {
                seq[i] = __pow_mod(modulus, root_of_2nth,
                                   __bit_rev_naive_16(i, log_dimension));
                seq_harvey[i] = ((u128)seq[i] << 64) / modulus;
            }
        } else {
            seq.resize(dimension * 2);
            seq_harvey.resize(dimension * 2);
            auto root_of_2nth_inv =
                __pow_mod(modulus, root_of_2nth, 2 * dimension - 1);
            for (size_t l = 0; l < log_dimension; l++) {
                auto start = (1 << l) - 1;
                auto power_index_factor = 1 << (log_dimension - l);
                for (size_t i = 0; i < (1 << l); i++) {
                    auto idx = start + i;
                    seq[idx] =
                        __pow_mod(modulus, root_of_2nth_inv,
                                  __bit_rev_naive_16(i, l) * power_index_factor);
                    seq_harvey[idx] = ((u128)seq[idx] << 64) / modulus;
                }
            }
            const u64 dimension_inv = modulus - ((modulus - 1) >> log_dimension);
            const u64 dimension_inv_harvey =
                ((u128)dimension_inv << 64) / modulus;
            for (size_t i = 0; i < dimension; i++) {
                auto temp = __pow_mod(modulus, root_of_2nth_inv, i);
                auto idx = i + dimension;
                seq[idx] = mul_mod_harvey_lazy(modulus, temp, dimension_inv,
                                               dimension_inv_harvey);
                seq[idx] -= (seq[idx] >= modulus) ? modulus : 0;
                seq_harvey[idx] = ((u128)seq[idx] << 64) / modulus;
            }

            shuffled_indices.resize(dimension);
            for (size_t i = 0; i < dimension; i++) {
                shuffled_indices[i] = __bit_rev_naive_16(i, log_dimension);
            }
        }
    }

    NTTFactors(const NTTFactors &copying) = default;

    NTTFactors(NTTFactors &&moving) = default;

    ~NTTFactors() {}

    std::vector<u64> seq;

    std::vector<u64> seq_harvey;

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

    size_t level, start, data_step, h, l;
    size_t idx = 1;
    u64 temp, zeta, zeta_harvey;
    for (level = 1, data_step = dimension; level <= log_dimension;
         level++, data_step >>= 1) {
        auto gap = data_step / 2;
        for (start = 0; start < dimension; start += data_step, idx++) {
            zeta = ntt_factors.seq[idx];
            zeta_harvey = ntt_factors.seq_harvey[idx];
            for (l = start; l < start + gap; l++) {
                h = l + gap;
                temp =
                    mul_mod_harvey_lazy(modulus, coeffs[h], zeta, zeta_harvey);
                coeffs[h] = coeffs[l] + 2 * modulus - temp;
                coeffs[l] = coeffs[l] + temp;
            }
        }
    }

    const u64 log_modulus = (u64)(log2(modulus) + 0.5);
    const u64 div_fix = (modulus >= (1ULL << log_modulus)) ? 1 : 0;
    for (size_t i = 0; i < dimension; i++) {
        coeffs[i] -= ((coeffs[i] >> log_modulus) - div_fix) * modulus;
    }
}

void intt_negacyclic_inplace_lazy(const size_t log_dimension, const u64 modulus,
                                  u64 values[]) {
    const size_t dimension = 1ULL << log_dimension;
    // generate or read from cache
    const auto &intt_factors =
        __find_or_create_ntt_factors(modulus, log_dimension, true);

    u64 values_shuffled[dimension];
    const auto &shuffled_indices = intt_factors.shuffled_indices;
    for (size_t i = 0; i < dimension; i++) {
        values_shuffled[i] = values[shuffled_indices[i]];
    }

    size_t level, start, data_step, h, l;
    size_t idx = 0;
    u64 temp, zeta, zeta_harvey;
    for (level = 1, data_step = dimension; level <= log_dimension;
         level++, data_step >>= 1) {
        auto gap = data_step / 2;
        for (start = 0; start < dimension; start += data_step, idx++) {
            zeta = intt_factors.seq[idx];
            zeta_harvey = intt_factors.seq_harvey[idx];
            for (l = start; l < start + gap; l++) {
                h = l + gap;
                temp = mul_mod_harvey_lazy(modulus, values_shuffled[h], zeta,
                                           zeta_harvey);
                values_shuffled[h] = values_shuffled[l] + 2 * modulus - temp;
                values_shuffled[l] = values_shuffled[l] + temp;
            }
        }
    }

    for (size_t i = 0; i < dimension; i++) {
        values[i] = values_shuffled[shuffled_indices[i]];
    }

    const u64 log_modulus = (u64)(log2(modulus) + 0.5);
    const u64 div_fix = (modulus >= (1ULL << log_modulus)) ? 1 : 0;
    idx++;
    for (size_t i = 0; i < dimension; i++, idx++) {
        values[i] -= ((values[i] >> log_modulus) - div_fix) * modulus;
        values[i] =
            mul_mod_harvey_lazy(modulus, values[i], intt_factors.seq[idx],
                                intt_factors.seq_harvey[idx]);
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
