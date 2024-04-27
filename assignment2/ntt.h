#include<cstdio>
#include<iostream>
#include<cstring>
#include<immintrin.h>

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

void ntt_dif_x4(int *a,int *b,int *ab,int *rt,int *irt);
void ntt_dit_x4(int *a,int *rt,int *irt,int limit,int level);
void ntt_dif_x4(int *a,int *rt,int *irt,int limit,int level);