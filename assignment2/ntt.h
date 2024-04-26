#include"params.h"
#include"readwrite.h"
#include"Math.h"

#include<cstdio>
#include<iostream>
#include<cstring>
#include<immintrin.h>

void ntt_common(int *a,int *b,int *ab,int *r);
void ntt_common(int *a,int *r,int limit,int type);
void ntt_Montgomery(int *a,int *b,int *ab,int *r);
void ntt_Montgomery(int *a,int *r,int limit,int type);
void ntt_dif(int *a,int *b,int *ab,int *rt,int *irt);
void ntt_dit(int *a,int *rt,int limit);
void ntt_dif(int *a,int *rt,int limit);