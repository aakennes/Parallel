#include<cstdio>
#include<iostream>
#include<cstring>

#include"Math.h"

void ntt_common(int *a,int *b,int *ab,int *r);
void ntt_common(int *a,int *r,int limit,int type);

void ntt_Montgomery(int *a,int *b,int *ab,int *r);
void ntt_Montgomery(int *a,int *r,int limit,int type);

void ntt_Montgomery_Mint(Mint *a,Mint *b,Mint *ab,int *r);
void ntt_Montgomery_Mint(Mint *a,int *r,int limit,int type);

void ntt_dif(int *a,int *b,int *ab,int *rt,int *irt);
void ntt_dit(int *a,int *rt,int limit);
void ntt_dif(int *a,int *rt,int limit);

void ntt_dif_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt);
void ntt_dit_Mint(Mint *a,Mint *rt,int limit);
void ntt_dif_Mint(Mint *a,Mint *rt,int limit);

void ntt_dif_x4(int *a,int *b,int *ab,int *rt,int *irt);
void ntt_dit_x4(int *a,int *rt,int *irt,int limit,int level);
void ntt_dif_x4(int *a,int *rt,int *irt,int limit,int level);

void ntt_dif_x4_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt);
void ntt_dit_x4_Mint(Mint *a,Mint *rt,Mint *irt,int limit,int level);
void ntt_dif_x4_Mint(Mint *a,Mint *rt,Mint *irt,int limit,int level);

void ntt_dif_x4_avx2(u32x8 *a,u32x8 *b,u32x8 *ab,u32x8 *rt,u32x8 *irt);
void ntt_dit_x4_avx2(u32x8 *a,u32x8 *rt,u32x8 *irt,int limit,int level);
void ntt_dif_x4_avx2(u32x8 *a,u32x8 *rt,u32x8 *irt,int limit,int level);