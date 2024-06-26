#include<cstdio>
#include<iostream>
#include<cstring>

#include"readwrite.h"
// #include"SIMD.h"


// #define MMint u32

void ntt_common(int *a,int *b,int *ab,int *r,int n);
void ntt_common(int *a,int *r,int limit,int type);

void ntt_barrett(int *a,int *b,int *ab,int *r,int n);
void ntt_barrett(int *a,int *r,int limit,int type);

void ntt_barrett_pre(int *a,int *b,int *ab,int *r,int n);
void ntt_barrett_pre(int *a,int *r,int limit,int type);

void ntt_galois(int *a,int *b,int *ab,int *r,int n);
void ntt_galois(int *a,int *r,int limit,int type);

void ntt_Montgomery(int *a,int *b,int *ab,int *r,int n);
void ntt_Montgomery(int *a,int *r,int limit,int type);

void ntt_Montgomery_Mint(Mint *a,Mint *b,Mint *ab,int *r,int n);
void ntt_Montgomery_Mint(Mint *a,int *r,int limit,int type);

void ntt_Montgomery_MMint(MMint *a,MMint *b,MMint *ab,int *r,int n);
void ntt_Montgomery_MMint(MMint *a,int *r,int limit,int type);

void ntt_dif(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit(int *a,int *rt,int limit);
void ntt_dif(int *a,int *rt,int limit);

void ntt_dif_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);
void ntt_dit_Mint(Mint *a,Mint *rt,int limit);
void ntt_dif_Mint(Mint *a,Mint *rt,int limit);

void ntt_dif_x4(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit_x4(int *a,int *rt,int *irt,int limit,int level);
void ntt_dif_x4(int *a,int *rt,int *irt,int limit,int level);

void ntt_dif_x4_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);
void ntt_dit_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level);
void ntt_dif_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level);

// void ntt_dif_x4_avx2(MMint *a,MMint *b,MMint *ab,MMint *rt,MMint *irt,Mintx8 *RT1,Mintx8 *RT2,Mintx8 *IRT1,Mintx8 *IRT2);
// void ntt_dit_x4_avx2(MMint *a,Mintx8 *rt,Mintx8 *irt,int limit,int level);
// void ntt_dif_x4_avx2(MMint *a,Mintx8 *rt,Mintx8 *irt,int limit,int level);