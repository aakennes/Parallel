#include<cstdio>
#include<iostream>
#include<cstring>
#include<omp.h>
#include<mpi.h>


#include"readwrite.h"
// #include"SIMD.h"
// #define MMint u32

void create_montgomerymodint32_type(MPI_Datatype *mpi_montgomerymodint32_type);

void ntt_common_mpi(int *a,int *b,int *ab,int *r,int n);
void ntt_common_mpi(int *a,int *r,int limit,int type, int my_rank, int processCounts);

void ntt_common_mpi_openmp(int *a,int *b,int *ab,int *r,int n);
void ntt_common_mpi_openmp(int *a,int *r,int limit,int type, int my_rank, int processCounts);

void ntt_Montgomery_Mint_mpi(Mint *a,Mint *b,Mint *ab,int *r,int n);
void ntt_Montgomery_Mint_mpi(Mint *a,int *r,int limit,int type, int my_rank, int processCounts, MPI_Datatype mpi_montgomerymodint32_type);

void ntt_Montgomery_Mint_mpi_openmp(Mint *a,Mint *b,Mint *ab,int *r,int n);
void ntt_Montgomery_Mint_mpi_openmp(Mint *a,int *r,int limit,int type, int my_rank, int processCounts, MPI_Datatype mpi_montgomerymodint32_type);

void ntt_dif_x4_mpi(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit_x4_mpi(int *a,int *dw,int imag,int limit,int level,int my_rank, int processCounts);
void ntt_dif_x4_mpi(int *a,int *dw,int imag,int limit,int level,int my_rank, int processCounts);

void ntt_dif_x4_mpi_openmp(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit_x4_mpi_openmp(int *a,int *dw,int imag,int limit,int level,int my_rank, int processCounts);
void ntt_dif_x4_mpi_openmp(int *a,int *dw,int imag,int limit,int level,int my_rank, int processCounts);

void ntt_dif_x4(int *a,int *b,int *ab,int *rt,int *irt,int n);
void ntt_dit_x4(int *a,int *rt,int *irt,int limit,int level);
void ntt_dif_x4(int *a,int *rt,int *irt,int limit,int level);

void ntt_dit_x4_Mint_mpi(Mint *a,Mint *dw,Mint imag,int limit,int level,int my_rank, int processCounts, Mint initial_w2[][4], MPI_Datatype mpi_montgomerymodint32_type);
void ntt_dif_x4_Mint_mpi(Mint *a,Mint *dw,Mint imag,int limit,int level, int my_rank, int processCounts, Mint initial_w2[][4], MPI_Datatype mpi_montgomerymodint32_type);
void ntt_dif_x4_Mint_mpi(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);

void ntt_dit_x4_Mint_mpi_openmp(Mint *a,Mint *dw,Mint imag,int limit,int level,int my_rank, int processCounts, Mint initial_w2[][4], MPI_Datatype mpi_montgomerymodint32_type);
void ntt_dif_x4_Mint_mpi_openmp(Mint *a,Mint *dw,Mint imag,int limit,int level, int my_rank, int processCounts, Mint initial_w2[][4], MPI_Datatype mpi_montgomerymodint32_type);
void ntt_dif_x4_Mint_mpi_openmp(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);

void ntt_dif_x4_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n);
void ntt_dit_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level);
void ntt_dif_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level);

// void ntt_dif_x4_avx2(MMint *a,MMint *b,MMint *ab,MMint *rt,MMint *irt,Mintx8 *RT1,Mintx8 *RT2,Mintx8 *IRT1,Mintx8 *IRT2);
// void ntt_dit_x4_avx2(MMint *a,Mintx8 *rt,Mintx8 *irt,int limit,int level);
// void ntt_dif_x4_avx2(MMint *a,Mintx8 *rt,Mintx8 *irt,int limit,int level);