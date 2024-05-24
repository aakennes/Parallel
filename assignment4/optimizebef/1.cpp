#include<iostream>
#include<algorithm>
#include <immintrin.h>

using i32 = int;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;

#pragma GCC target("avx2")
#define Vec(sz, T) __attribute((vector_size(sz))) T
#define IL __inline__ __attribute__((always_inline))
#define RC(T, x) reinterpret_cast<T>(x)



constexpr u32 mod = 104857601;
constexpr u32 G = 3;
constexpr int sta_l_MB = 64;

namespace Montgo{
	struct Mont32{
		u32 Mod, Mod2, Inv, NInv, R2;
		constexpr Mont32(u32 n):Mod(n), Mod2(n << 1), Inv(1), NInv(), R2((-u64(n)) % n){
			for (int i = 0; i < 5; ++i){Inv *= 2 - n * Inv;}
			NInv = -Inv;
		}
		constexpr IL u32 reduce (u64 x)const{return (x + u64(u32(x) * NInv) * Mod) >> 32;}
		constexpr IL u32 reduce_s (u64 x)const{
			u32 r = (x >> 32) - ((u64(u32(x) * Inv) * Mod) >> 32);
			return r >> 31 ? r + Mod : r;
		}
		constexpr IL u32 mul(u32 x, u32 y)const{return reduce(u64(x) * y);}
		constexpr IL u32 mul_s(u32 x, u32 y)const{return reduce_s(u64(x) * y);}
		constexpr IL u32 In(u32 x)const{return mul(x, R2);}
		constexpr IL u32 In_s(u32 x)const{return mul_s(x, R2);}
		constexpr IL u32 Out(u32 x)const{
			u32 r = (x + (u64(x * NInv) * Mod)) >> 32;
			return __builtin_expect(r < Mod, 1) ? r : r - Mod;
		}
	};
}

namespace field_Z{
	constexpr Montgo::Mont32 Space(mod);
	constexpr u32 mod2 = Space.Mod2;
	constexpr IL u32 shrink(u32 x){return x >= mod ? x - mod : x;}
	constexpr IL u32 dilate2(u32 x){return x >> 31 ? x + mod2 : x;}
	using Z = u32;
	constexpr bool isgood(Z x){return x < mod2;}
	constexpr IL Z InZ(u32 x){return Space.In(x);}
	constexpr IL Z InZs(u32 x){return Space.In_s(x);}
	constexpr Z zero_Z(InZs(0)), one_Z(InZs(1));
	constexpr Z not_exist_Z(-1);
	constexpr IL u32 OutZ(Z x){return Space.Out(x);}
	namespace calc{
		constexpr IL Z addZ(Z x, Z y){return dilate2(x + y - mod2);}
		constexpr IL Z subZ(Z x, Z y){return dilate2(x - y);}
		constexpr IL Z mulZ(Z x, Z y){return Space.mul(x, y);}
		constexpr Z powZ(Z a, u32 b, Z r = one_Z){
			while (b){
				if (b & 1){r = mulZ(r, a);}
				a = mulZ(a, a), b >>= 1;
			}
			return r;
		}
		constexpr IL Z invZ(Z x){return powZ(x, mod - 2);}
		constexpr IL Z divZ(Z x, Z y){return powZ(y, mod - 2, x);}
		template<bool strict = true>constexpr IL Z mulZs(Z x, Z y){
			if constexpr(strict) return Space.mul_s(x, y);
			return mulZ(x, y);
		}
		constexpr Z negZ(Z x){return (!x - 1) & (mod2 - x);}
	}
	template<int fixes = 0>constexpr Z trans(Z x){
		constexpr Z o = fixes > 0 ? calc::powZ(Space.R2, fixes) : calc::powZ(1, -fixes);
		return calc::mulZs(x, o);
	}
	constexpr Z trans(Z x, int fixes){
		return calc::mulZs(x, fixes > 0 ? calc::powZ(Space.R2, fixes) : calc::powZ(1, -fixes));
	}
}
namespace SIMD{
	using i32x8 = Vec(32, i32);
	using u32x8 = Vec(32, u32);
	using i64x4 = Vec(32, i64);
	using u64x4 = Vec(32, u64);
	using I256 = __m256i;
	constexpr IL u32x8 setu32x8(u32 x){return (u32x8){x, x, x, x, x, x, x, x};}
	template<int typ>IL u32x8 shuffle(const u32x8 &x){return RC(u32x8, _mm256_shuffle_epi32(RC(I256, x), typ));}
	template<int typ>IL u32x8 blend(const u32x8 &x, const u32x8 &y){return RC(u32x8, _mm256_blend_epi32(RC(I256, x), RC(I256, y), typ));}
	IL I256 swaplohi128(const I256 &x){return _mm256_permute2x128_si256(x, x, 1);}	
	IL u32x8& x8(u32 *data){return *RC(u32x8* ,data);}
	IL const u32x8& x8(const u32 *data){return *RC(const u32x8*, data);}
	IL I256 loadu(const void* data){return _mm256_loadu_si256(RC(const __m256i_u*, data));}
	IL void storeu(const I256 &x, void* data){return _mm256_storeu_si256(RC(__m256i_u*, data), x);}
	IL u64x4 mulu32x8_fus(const u32x8 &x, const u32x8 &y){return RC(u64x4, _mm256_mul_epu32(RC(I256, x), RC(I256, y)));}
}
namespace field_Z{
	using SIMD::x8;
	using SIMD::u32x8;
	using SIMD::setu32x8;
	using Zx8 = u32x8;
	constexpr u32x8 modx8 = setu32x8(mod), mod2x8 = setu32x8(mod2), NInvx8 = setu32x8(Space.NInv);
	constexpr Zx8 one_Zx8 = setu32x8(one_Z), zerox8 = setu32x8(0u);
	IL Zx8 dilate2x8(const Zx8 &x){return x + (RC(Zx8, RC(SIMD::i32x8, x) < RC(SIMD::i32x8, zerox8)) & mod2x8);}
	IL Zx8 shrinkx8(const Zx8 &x){return x - ((x >= modx8) & modx8);}
	namespace calc{
		using namespace SIMD;
		IL Zx8 addZx8(const Zx8 &x, const Zx8 &y){return dilate2x8(x + y - mod2x8);}
		IL Zx8 subZx8(const Zx8 &x, const Zx8 &y){return dilate2x8(x - y);}
		IL Zx8 mulZx8(const Zx8 &x, const Zx8 &y){
			u32x8 z = (NInvx8 * x * y);
			return blend<0xaa>(RC(u32x8, (mulu32x8_fus(x, y) + mulu32x8_fus(z, modx8)) >> 32), RC(u32x8, (mulu32x8_fus(u32x8(u64x4(x) >> 32), u32x8(u64x4(y) >> 32)) + mulu32x8_fus(shuffle<0xf5>(z), modx8))));
		}
		IL Zx8 powZx8(Zx8 x, u32 b, Zx8 r = one_Zx8){
			while(b){
				if(b & 1){r = mulZx8(r, x);}
				x = mulZx8(x, x), b >>= 1;
			}
			return r;
		}
		IL Zx8 invZx8(const Zx8 &x){return powZx8(x, mod - 2);}
		IL Zx8 divZx8(const Zx8 &x, const Zx8 &y){return powZx8(y, mod - 2, x);}
		template<bool strict = true>IL Zx8 mulZsx8(const Zx8 &x, const Zx8 &y){
			if constexpr (strict){
				u32x8 z = (NInvx8 * x * y);
				z = blend<0xaa>(RC(u32x8, (mulu32x8_fus(x, y) + mulu32x8_fus(z, modx8)) >> 32), RC(u32x8, (mulu32x8_fus(u32x8(u64x4(x) >> 32), u32x8(u64x4(y) >> 32)) + mulu32x8_fus(shuffle<0xf5>(z), modx8)))) - modx8;
				return z + (RC(Zx8, RC(i32x8, z) < RC(i32x8, zerox8)) & modx8);
			}
			return mulZx8(x, y);
		}
		IL Zx8 negZx8(const Zx8 &x){return (x != zerox8) & (mod2x8 - x);}
	}
}
#undef Vec

namespace poly{
	//多项式基础
	namespace poly_base{
		using namespace field_Z;
		constexpr int bit_up(int x){return 1 << (32 - __builtin_clz(x));}
		constexpr int cro_32(int x){return __builtin_ctz(~x);}
		inline Z *to_align(void *mem){return RC(Z*, ((RC(u64, mem) + 31) >> 5) << 5);}
		inline bool is_align(const void* mem){return (RC(u64, mem) & 31) == 0;}
		//用于申请临时空间
        //align用于申请对齐内存
		namespace mem_helper{
			char _mem[sta_l_MB << 20];
			void *now = _mem;
			struct pre_aloc{
				void* t;
				pre_aloc(){t = now;}
				~pre_aloc(){now = t;}
			};
			void *aloc(size_t l){
				void *r = now;
				now = RC(char*, r) + l;
				return r;
			}
			Z *alocP(int l){
				Z *r = to_align(now);
				now = r + l;
				return r;
			}
		}
		using mem_helper::pre_aloc;
		using mem_helper::alocP;
	}using namespace poly_base;
	#define flx(nam, opt) void nam(Z *A,int l,const Z *B){int i=0;for(;i+7<l;i+=8){x8(A + i)=calc::opt##Zx8(x8(A+i),x8(B+i));}for(;i<l;++i){A[i]=calc::opt##Z(A[i],B[i]);}}  void nam(const Z *A,const Z *B,int l,Z *C){int i=0;for(;i+7<l;i+=8){x8(C+i)=calc::opt##Zx8(x8(A+i),x8(B+i));}for(;i<l;++i){C[i]=calc::opt##Z(A[i],B[i]);}}
		flx(dot, mul)
	#undef flx
}
namespace poly{
	//快速数论变换.基础 by QedDust413
	namespace f_n_t_t_base{
		using namespace calc;
		//A[] *= t
		template<bool strict = false>void mul_t(Z *A, int l, Z t){
			for(int j = 0; j < l; ++j){
				A[j] = mulZs<strict>(A[j], t);
			}
		}
		constexpr int mp2 = __builtin_ctz(mod - 1);
		constexpr Z _g(InZ(G));
		//2^n本原单位根表
		struct P_R_Tab{
			Z t[mp2 + 1];
			constexpr P_R_Tab(Z G) : t(){
				t[mp2] = powZ(G, (mod - 1) >> mp2);
				for(int i = mp2 - 1; ~i; --i){
					t[i] = mulZs(t[i + 1], t[i + 1]);
				}
			}
			constexpr Z operator [] (int i) const {
				return t[i];
			}
		};
		constexpr P_R_Tab rt1(_g), rt1_I(invZ(_g));
		//为传统的dif/dit提供信息
		struct ntt_info_base2{
			Z rt2[mp2 - 1], rt2_I[mp2 - 1];
			constexpr ntt_info_base2() : rt2(), rt2_I(){
				Z pr = one_Z, pr_I = one_Z;
				for(int i = 0; i < mp2 - 1; ++i){
					rt2[i] = mulZs(pr, rt1[i + 2]), rt2_I[i] = mulZs(pr_I, rt1_I[i + 2]);
					pr = mulZs(pr, rt1_I[i + 2]), pr_I = mulZs(pr_I, rt1[i + 2]);
				
				}
				
			}
		};
		//为以4为基的指令集加速dif/dit提供信息
		struct ntt_info_base4x8{
			Zx8 rt3x8[mp2 - 2], rt3x8_I[mp2 - 2], rt4ix8[mp2 - 3], rt4ix8_I[mp2 - 3], pr2, pr4, pr2_I, pr4_I;
			constexpr ntt_info_base4x8():rt3x8(), rt3x8_I(), rt4ix8(), rt4ix8_I(), pr2(), pr4(), pr2_I(), pr4_I()
			{   
				Z pr = one_Z, pr_I = one_Z;
				for(int i = 0; i < mp2 - 2; ++i){
					rt3x8[i] = setu32x8(mulZs(pr, rt1[i + 3]));
					rt3x8_I[i] = setu32x8(mulZs(pr_I, rt1_I[i + 3]));
					pr = mulZs(pr, rt1_I[i + 3]), pr_I = mulZs(pr_I, rt1[i + 3]);
				}
				pr = one_Z, pr_I = one_Z;
				for(int i = 0; i < mp2 - 3; ++i){
					{
						Z a0(one_Z), a1(mulZs(pr, rt1[i + 4])), a2(mulZs(a1, a1)), a3(mulZs(a1, a2)), 
						a4(mulZs(a1, a3)), a5(mulZs(a1, a4)), a6(mulZs(a1, a5)), a7(mulZs(a1, a6));
						rt4ix8[i] = (Zx8){a0, a1, a2, a3, a4, a5, a6, a7};
						// std::cout<<a1<<" ";
					}
					{
						Z a0(one_Z), a1(mulZs(pr_I, rt1_I[i + 4])), a2(mulZs(a1, a1)), a3(mulZs(a1, a2)), 
						a4(mulZs(a1, a3)), a5(mulZs(a1, a4)), a6(mulZs(a1, a5)), a7(mulZs(a1, a6));
						rt4ix8_I[i] = (Zx8){a0, a1, a2, a3, a4, a5, a6, a7};
					}
					pr = mulZs(pr, rt1_I[i + 4]), pr_I = mulZs(pr_I, rt1[i + 4]);
				}
				pr2 = (Zx8){one_Z, one_Z, one_Z, rt1[2], one_Z, one_Z, one_Z, rt1[2]};
				pr2_I = (Zx8){one_Z, one_Z, one_Z, rt1_I[2], one_Z, one_Z, one_Z, rt1_I[2]};
				pr4 = (Zx8){one_Z, one_Z, one_Z, one_Z, one_Z, rt1[3], rt1[2], mulZs(rt1[2], rt1[3])};
				pr4_I = (Zx8){one_Z, one_Z, one_Z, one_Z, one_Z, rt1_I[3], rt1_I[2], mulZs(rt1_I[2], rt1_I[3])};
			}
		};
		template<int typ>IL u32x8 Neg(const u32x8 &x){return blend<typ>(x, mod2x8 - x);}
	}
	//fast_number_theoretic_transform by QedDust413
    namespace f_n_t_t{
		using namespace f_n_t_t_base;
		constexpr ntt_info_base2 ib2;
		//dif-ntt 非注重性能的实现
		template<bool strict = false, int fixes = 0>void dif_base2(Z *A, int lim){
			for(int L = lim >> 1, R = lim; L; L >>= 1, R >>= 1){
				Z r = one_Z;
				for(int i = 0, k = 0; i < lim; i += R, ++k){
					for(int j = 0; j < L; ++j){
						Z x = dilate2(A[i + j] - mod2) , y = mulZ(r, A[i + j + L]);
						A[i + j] = x + y, A[i + j + L] = x - y + mod2;
					}
					r = mulZs(r, ib2.rt2[cro_32(k)]);
				}
			}
			if constexpr(fixes){
				mul_t<strict>(A, lim, trans<fixes>(one_Z));
			}
		}
		//dit-intt 非注重性能的实现
		template<bool strict = false, int fixes = 0>void dit_base2(Z *A, int lim){
			for(int L = 1, R = 2; L < lim; L <<= 1, R <<= 1){
				Z r = one_Z;
				for(int i = 0, k = 0; i < lim; i += R, ++k){
					for(int j = 0; j < L; ++j){
						Z x = A[i + j], y = A[i + j + L];
						A[i + j] = addZ(x, y), A[i + j + L] = mulZ(x - y + mod2, r);
					}
					r = mulZs(r, ib2.rt2_I[cro_32(k)]);
				}
			}
			mul_t<strict>(A, lim, trans<fixes + 1>(mod - ((mod - 1) / lim)));
		}
		constexpr Zx8 imagx8 = setu32x8(rt1[2]), imag_Ix8 = setu32x8(rt1_I[2]);
		constexpr ntt_info_base4x8 iab4;
		//以4为基的指令集加速DIF式NTT 变换长度应当至少为8 且给出的指针需对32对齐
        template<bool strict = false, int fixes = 0>void dif_base4x8(Z *A, int lim){
			int n = lim >> 3, L = n >> 1;
			Zx8 *f = RC(Zx8*, A);
			if(__builtin_ctz(n) & 1){
				for(int j = 0; j < L; ++j){
					Zx8 x = f[j], y = f[j + L];
					f[j] = x + y, f[j + L] = x - y + mod2x8;
				}
				L >>= 1;
			}
			L >>= 1;
			for(int R = L << 2; L; L >>= 2, R >>= 2){
				Zx8 r = one_Zx8, img = imagx8;
				for(int i = 0, k = 0; i < n; i += R, ++k){
					Zx8 r2 = mulZsx8(r, r), r3 = mulZsx8(r2, r);
					for(int j = 0; j < L; ++j){
						Zx8 f0 = dilate2x8(f[i + j + 0 * L] - mod2x8);
						Zx8 f1 = mulZx8(f[i + j + 1 * L], r);
						Zx8 f2 = mulZx8(f[i + j + 2 * L], r2);
						Zx8 f3 = mulZx8(f[i + j + 3 * L], r3);
						Zx8 f1f3 = mulZx8(f1 - f3 + mod2x8, img);
						Zx8 f02 = addZx8(f0, f2);
						Zx8 f13 = addZx8(f1, f3);
						Zx8 f_02 = subZx8(f0, f2);
						f[i + j + 0 * L] = f02 + f13;
						f[i + j + 1 * L] = f02 - f13 + mod2x8;
						f[i + j + 2 * L] = f_02 + f1f3;
						f[i + j + 3 * L] = f_02 - f1f3 + mod2x8;
					}
					r = mulZsx8(r, iab4.rt3x8[cro_32(k)]);
				}
			}
			
			{
				constexpr Zx8 _r = setu32x8(trans<fixes>(one_Z));
				Zx8 r = _r, pr4 = iab4.pr4, pr2 = iab4.pr2;
				
				for(int i = 0; i < n; ++i){
					Zx8& fi = f[i];
					fi = mulZx8(fi, r), fi = Neg<0xf0>(fi) + RC(Zx8, swaplohi128(RC(I256, fi)));
					fi = mulZx8(fi, pr4), fi = Neg<0xcc>(fi) + shuffle<0x4e>(fi);
					fi = mulZx8(fi, pr2), fi = addZx8(Neg<0xaa>(fi), shuffle<0xb1>(fi));
					if constexpr(strict){fi = shrinkx8(fi);}
					r = mulZsx8(r, iab4.rt4ix8[cro_32(i)]);
				}
			}
        }
		//以4为基的指令集加速DIT式INTT 变换长度应当至少为8 且给出的指针需对32对齐
		template<bool strict = false, int fixes = 0>void dit_base4x8(Z *A, int lim){
			int n = lim >> 3, L = 1;
			Zx8 *f = RC(Zx8*, A);
			{
				Zx8 r = setu32x8(trans<fixes + 1>(mod - ((mod - 1) / lim))), pr4 = iab4.pr4_I, pr2 = iab4.pr2_I;
				for(int i = 0; i < n; ++i){
					Zx8& fi = f[i];
					//
					//
					//
					//
					fi = Neg<0xaa>(fi) + shuffle<0xb1>(fi), fi = mulZx8(fi, pr2);
					fi = Neg<0xcc>(fi) + shuffle<0x4e>(fi), fi = mulZx8(fi, pr4);
					fi = Neg<0xf0>(fi) + RC(Zx8, swaplohi128(RC(I256, fi))), fi = mulZx8(fi, r);
					r = mulZsx8(r, iab4.rt4ix8_I[cro_32(i)]);
				}
			}

			for (int R = L << 2; L < (n >> 1) ; L <<= 2, R <<= 2){
				Zx8 r = one_Zx8, img = imag_Ix8;
				for(int i = 0, k = 0; i < n; i += R, ++k){
					Zx8 r2 = mulZsx8(r, r), r3 = mulZsx8(r2, r);
					for(int j = 0; j < L; ++j){
						Zx8 f0 = f[i + j + 0 * L];
						Zx8 f1 = f[i + j + 1 * L];
						Zx8 f2 = f[i + j + 2 * L];
						Zx8 f3 = f[i + j + 3 * L];
						Zx8 f2f3 = mulZx8((f2 - f3 + mod2x8), img);
						Zx8 f01 = addZx8(f0, f1);
						Zx8 f23 = addZx8(f2, f3);
						Zx8 f_01 = subZx8(f0, f1);
						f[i + j + 0 * L] = addZx8(f01, f23);
						f[i + j + 1 * L] = mulZx8(f_01 + f2f3 ,r);
						f[i + j + 2 * L] = mulZx8(f01 - f23 + mod2x8, r2);
						f[i + j + 3 * L] = mulZx8(f_01 - f2f3 + mod2x8, r3);
					}
					r = mulZsx8(r, iab4.rt3x8_I[cro_32(k)]);
				}
			}
			if(__builtin_ctz(n) & 1){
				for(int j = 0; j < L; ++j){
					Zx8 x = f[j], y = f[j + L];
					f[j] = addZx8(x, y), f[j + L] = subZx8(x, y);
				}
			}
			// if constexpr (strict){
			// 	for(int i = 0; i < n; ++i){
			// 		f[i] = shrinkx8(f[i]);
			// 	}
			// }
		}
		//1E(lim)
		template<bool strict = false, int fixes = 0>void dif(Z *A, int lim){
			lim >= 16 ? dif_base4x8<strict, fixes>(A, lim) : dif_base2<strict, fixes>(A, lim);
		}
		//1E(lim)
		template<bool strict = false, int fixes = 0>void dit(Z *A,int lim){
			lim >= 16 ? dit_base4x8<strict, fixes>(A, lim) : dit_base2<strict, fixes>(A, lim);
		}
    }
	using f_n_t_t::dif;
	using f_n_t_t::dit;
}
#undef Stati
namespace Command{
	void cut_string(){
		_Exit(0);
	}
}

using poly::alocP;
using poly::pre_aloc;
using namespace field_Z;


#include<iostream>
using namespace std;
void solve(){
	freopen("data/datain/NTT256_1409.in","r",stdin);
	int n, m;
	cin >> n >> m;
	int limit = poly::bit_up(n + m);
    // printf("%d\n",limit);
	auto F = alocP(limit), G = alocP(limit);
	for(int i = 0; i <= n; ++i){cin >> F[i];}
	std::fill(F + n + 1, F + limit, zero_Z);
	for(int i = 0; i <= m; ++i){cin >> G[i];}
	std::fill(F + m + 1, F + limit, zero_Z);
	poly::dif<false, 1>(F, limit), poly::dif<false, 1>(G, limit), poly::dot(F, limit, G), poly::dit<true, -1>(F, limit);
	for(int i = 0; i <= n + m; ++i){cout << F[i] << ' ';}
	// for(int i = 0; i < 23 - 1; ++i){
	// 	std::cout<<poly::f_n_t_t::iab4.rt3x8[i][0]<<" ";
	// 	// puts("");
	// }
}

int main(){
	std::cin.tie(nullptr) -> sync_with_stdio(false);
	solve();
	return 0;
}