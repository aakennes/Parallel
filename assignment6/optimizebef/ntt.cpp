#include"ntt.h"

#define p Original_Q

void ntt_common_mpi(int *a,int *r,int limit,int type, int my_rank, int processCounts){
    if (my_rank == 0) {
        for (int i = 0; i < limit; i++) {
            if (i < r[i]) {
                std::swap(a[i], a[r[i]]);
            }
        }
    }

    // 动态设置 SUBARRAY_SIZE
    int SUBARRAY_SIZE = limit / processCounts;
    int *local_a = new int[SUBARRAY_SIZE];  // 动态分配子数组大小

    int wnn = 3;
    int invwnn = qpow(wnn, p - 2, p);

    // 使用MPI_Scatter将数组a的子数组分发给各个进程
    MPI_Scatter(a, SUBARRAY_SIZE, MPI_INT, local_a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    for (int mid = 1; mid < SUBARRAY_SIZE; mid <<= 1) {
        int Wn = qpow(type == 1 ? wnn : invwnn, (p - 1) / (mid << 1), p);
        for (int j = 0; j < SUBARRAY_SIZE; j += (mid << 1)) {
            int w = 1;
            for (int k = 0; k < mid; k++, w = (1LL * w * Wn) % p) {
                int x = local_a[j + k], y = 1LL * w * local_a[j + k + mid] % p;
                local_a[j + k] = (1LL * x + y) % p;
                local_a[j + k + mid] = (1LL * x - y + p) % p;
            }
        }
    }
    // 使用MPI_Gather收集各个进程处理后的子数组到主进程
    MPI_Gather(local_a, SUBARRAY_SIZE, MPI_INT, a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    
    // 将小规模数据处理完毕后，在所有 0 process进行大规模数据处理后广播
    if(my_rank == 0)
    {
        for (int mid = SUBARRAY_SIZE; mid < limit; mid <<= 1) {
            int Wn = qpow(type == 1 ? wnn : invwnn, (p - 1) / (mid << 1), p);
            for (int j = 0; j < limit; j += (mid << 1)) {
                int w = 1;
                for (int k = 0; k < mid; k++, w = (1LL * w * Wn) % p) {
                    int x = a[j + k], y = 1LL * w * a[j + k + mid] % p;
                    a[j + k] = (1LL * x + y) % p;
                    a[j + k + mid] = (1LL * x - y + p) % p;
                }
            }
        }
    }

    // 此时只有 0 process 是正确的全局变量，进行广播，统一变量
    MPI_Bcast(a, limit, MPI_INT, 0, MPI_COMM_WORLD);

    delete[] local_a;
}

void ntt_common_mpi(int *a,int *b,int *ab,int *r,int n){
    int L = 0, i = 0;
    int limit = 1;
    while (limit <= 2 * n - 2) {
        limit <<= 1;
        L++;
    }

    for (i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }

    int my_rank, processCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCounts);

    // ntt_common 内部会将 process 0 初始化的 a、b 分发
    ntt_common_mpi(a, r, limit, 1, my_rank, processCounts); // 计算后，a gather in 0 process
    ntt_common_mpi(b, r, limit, 1, my_rank, processCounts); // 计算后, b gather in 0 process

    // process 0 处理点值
    if (my_rank == 0) {
        for (i = 0; i < limit; i++) {
            ab[i] = 1LL * a[i] * b[i] % p;
            // printf("%d ",ab[i]);
        }
    }
    // ntt_common会把 process 0 处理点值的结果 ab 分发
    ntt_common_mpi(ab, r, limit, -1, my_rank, processCounts); // INTT后, ab gather in 0 process

    // 对最终结果进行处理，得到正确结果
    int invn = qpow(limit, p - 2, p);
    for (i = 0; i < 2 * n; i++) {
        ab[i] = (1LL * ab[i] * invn) % p;
        // std::cout << ab[i] << " ";
    }
}

void ntt_common_mpi_openmp(int *a,int *r,int limit,int type, int my_rank, int processCounts){
    if (my_rank == 0) {
        for (int i = 0; i < limit; i++) {
            if (i < r[i]) {
                std::swap(a[i], a[r[i]]);
            }
        }
    }

    // 动态设置 SUBARRAY_SIZE
    int SUBARRAY_SIZE = limit / processCounts;
    int *local_a = new int[SUBARRAY_SIZE];  // 动态分配子数组大小

    int wnn = 3;
    int invwnn = qpow(wnn, p - 2, p);

    // 使用MPI_Scatter将数组a的子数组分发给各个进程
    MPI_Scatter(a, SUBARRAY_SIZE, MPI_INT, local_a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);

    #pragma omp parallel num_threads(thread_count)
    {
        for (int mid = 1; mid < SUBARRAY_SIZE; mid <<= 1) {
            int Wn = qpow(type == 1 ? wnn : invwnn, (p - 1) / (mid << 1), p);
            #pragma omp for schedule(static)
            for (int j = 0; j < SUBARRAY_SIZE; j += (mid << 1)) {
                int w = 1;
                for (int k = 0; k < mid; k++, w = (1LL * w * Wn) % p) {
                    int x = local_a[j + k], y = 1LL * w * local_a[j + k + mid] % p;
                    local_a[j + k] = (1LL * x + y) % p;
                    local_a[j + k + mid] = (1LL * x - y + p) % p;
                }
            }
        }
    }


    // 使用MPI_Gather收集各个进程处理后的子数组到主进程
    MPI_Gather(local_a, SUBARRAY_SIZE, MPI_INT, a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    
    // 将小规模数据处理完毕后，在所有 0 process进行大规模数据处理后广播
    if(my_rank == 0)
    {   
        #pragma omp parallel num_threads(thread_count)
        {
            for (int mid = SUBARRAY_SIZE; mid < limit; mid <<= 1) {
                int Wn = qpow(type == 1 ? wnn : invwnn, (p - 1) / (mid << 1), p);
                #pragma omp for schedule(static)
                for (int j = 0; j < limit; j += (mid << 1)) {
                    int w = 1;
                    for (int k = 0; k < mid; k++, w = (1LL * w * Wn) % p) {
                        int x = a[j + k], y = 1LL * w * a[j + k + mid] % p;
                        a[j + k] = (1LL * x + y) % p;
                        a[j + k + mid] = (1LL * x - y + p) % p;
                    }
                }
            }   
        }
    }

    // 此时只有 0 process 是正确的全局变量，进行广播，统一变量
    MPI_Bcast(a, limit, MPI_INT, 0, MPI_COMM_WORLD);

    delete[] local_a;
}

void ntt_common_mpi_openmp(int *a,int *b,int *ab,int *r,int n){
    int L = 0, i = 0;
    int limit = 1;
    while (limit <= 2 * n - 2) {
        limit <<= 1;
        L++;
    }

    for (i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }

    int my_rank, processCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCounts);

    // ntt_common 内部会将 process 0 初始化的 a、b 分发
    ntt_common_mpi_openmp(a, r, limit, 1, my_rank, processCounts); // 计算后，a gather in 0 process
    ntt_common_mpi_openmp(b, r, limit, 1, my_rank, processCounts); // 计算后, b gather in 0 process

    // process 0 处理点值
    if (my_rank == 0) {
        for (i = 0; i < limit; i++) {
            ab[i] = 1LL * a[i] * b[i] % p;
            // printf("%d ",ab[i]);
        }
    }
    // ntt_common会把 process 0 处理点值的结果 ab 分发
    ntt_common_mpi_openmp(ab, r, limit, -1, my_rank, processCounts); // INTT后, ab gather in 0 process

    // 对最终结果进行处理，得到正确结果
    int invn = qpow(limit, p - 2, p);
    for (i = 0; i < 2 * n; i++) {
        ab[i] = (1LL * ab[i] * invn) % p;
        // std::cout << ab[i] << " ";
    }
}

void create_montgomerymodint32_type(MPI_Datatype *mpi_montgomerymodint32_type) {
    MPI_Aint offsets[1];
    int block_lengths[1] = {1};
    MPI_Datatype types[1] = {MPI_UINT32_T};

    offsets[0] = offsetof(Mint, v);

    MPI_Type_create_struct(1, block_lengths, offsets, types, mpi_montgomerymodint32_type);
    MPI_Type_commit(mpi_montgomerymodint32_type);
}

void ntt_Montgomery_Mint_mpi(Mint *a,int *r,int limit,int type, int my_rank, int processCounts, MPI_Datatype mpi_montgomerymodint32_type){
    
    if (my_rank == 0) {
        for (int i = 0; i < limit; i++) {
            if (i < r[i]) {
                std::swap(a[i], a[r[i]]);
            }
        }
    }

    int SUBARRAY_SIZE = limit / processCounts;
    Mint *local_a = new Mint[SUBARRAY_SIZE];  // 动态分配子数组大小

    int wnn = 3;
    Mint Mwnn = wnn;
    Mint Minvwnn = 1/Mwnn;

    MPI_Scatter(a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);

	for(int mid = 1; mid < SUBARRAY_SIZE; mid <<= 1) {	
        Mint Moperator1 = type == 1 ? Mwnn : Minvwnn;
        int64_t Moperator2 = (p - 1) / (mid << 1);
		Mint MWn = Moperator1.pow(Moperator2);
		for(int j = 0; j < SUBARRAY_SIZE; j += (mid << 1)) {
			Mint Mw = 1;
			for(int k = 0; k < mid; k++, Mw = Mw * MWn) {                   
                Mint x = local_a[j + k];
				Mint y = Mw * local_a[j + k + mid];                    
				local_a[j + k] = (x + y),
				local_a[j + k + mid] = (x - y);
			}
		}
	}
    // 使用MPI_Gather收集各个进程处理后的子数组到主进程
    MPI_Gather(local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);
    
    // 将小规模数据处理完毕后，在所有 0 process进行大规模数据处理后广播
    if(my_rank == 0)
    {
        for (int mid = SUBARRAY_SIZE; mid < limit; mid <<= 1) {
            Mint Moperator1 = type == 1 ? Mwnn : Minvwnn;
            int64_t Moperator2 = (p - 1) / (mid << 1);
            Mint MWn = Moperator1.pow(Moperator2);
            for(int j = 0; j < limit; j += (mid << 1)) {
                Mint Mw = 1;
                for(int k = 0; k < mid; k++, Mw = Mw * MWn) {                   
                    Mint x = a[j + k];
                    Mint y = Mw * a[j + k + mid];                    
                    a[j + k] = (x + y),
                    a[j + k + mid] = (x - y);
                }
            }
        }
    }

    MPI_Bcast(a, limit, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);

    delete[] local_a;
}

void ntt_Montgomery_Mint_mpi(Mint *a,Mint *b,Mint *ab,int *r,int n){
    MPI_Datatype mpi_montgomerymodint32_type;
    create_montgomerymodint32_type(&mpi_montgomerymodint32_type);
    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
	for(i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }

    int my_rank, processCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCounts);
    ntt_Montgomery_Mint_mpi(a,r,limit, 1, my_rank, processCounts, mpi_montgomerymodint32_type);
    ntt_Montgomery_Mint_mpi(b,r,limit, 1, my_rank, processCounts, mpi_montgomerymodint32_type);
    // printf("%u\n",a[0].getv());
	for(i = 0; i < limit; i++){
        ab[i]= a[i] * b[i];
        // std::cout<<ab[i]<<' ';
    } 
	ntt_Montgomery_Mint_mpi(ab,r,limit, -1, my_rank, processCounts, mpi_montgomerymodint32_type);
    Mint Mlimit = limit;
    Mint invn = (1 / Mlimit);
    // printf("%u %u\n",Mlimit.getv(),invn.getv());
    for(i = 0; i < 2*n; i++){
        ab[i] = (ab[i] * invn);
        // std::cout<<ab[i]<<' ';
    } 
    
}

void ntt_Montgomery_Mint_mpi_openmp(Mint *a,int *r,int limit,int type, int my_rank, int processCounts, MPI_Datatype mpi_montgomerymodint32_type){
    
    if (my_rank == 0) {
        for (int i = 0; i < limit; i++) {
            if (i < r[i]) {
                std::swap(a[i], a[r[i]]);
            }
        }
    }

    int SUBARRAY_SIZE = limit / processCounts;
    Mint *local_a = new Mint[SUBARRAY_SIZE];  // 动态分配子数组大小

    int wnn = 3;
    Mint Mwnn = wnn;
    Mint Minvwnn = 1/Mwnn;

    MPI_Scatter(a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);


    #pragma omp parallel num_threads(thread_count)
    {

        for(int mid = 1; mid < SUBARRAY_SIZE; mid <<= 1) {	
            Mint Moperator1 = type == 1 ? Mwnn : Minvwnn;
            int64_t Moperator2 = (p - 1) / (mid << 1);
            Mint MWn = Moperator1.pow(Moperator2);
            #pragma omp for schedule(static)
            for(int j = 0; j < SUBARRAY_SIZE; j += (mid << 1)) {
                Mint Mw = 1;
                for(int k = 0; k < mid; k++, Mw = Mw * MWn) {                   
                    Mint x = local_a[j + k];
                    Mint y = Mw * local_a[j + k + mid];                    
                    local_a[j + k] = (x + y),
                    local_a[j + k + mid] = (x - y);
                }
            }
        }

    }

    // 使用MPI_Gather收集各个进程处理后的子数组到主进程
    MPI_Gather(local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);
    
    // 将小规模数据处理完毕后，在所有 0 process进行大规模数据处理后广播
    if(my_rank == 0)
    {
        #pragma omp parallel num_threads(thread_count)
        {

            for (int mid = SUBARRAY_SIZE; mid < limit; mid <<= 1) {
                Mint Moperator1 = type == 1 ? Mwnn : Minvwnn;
                int64_t Moperator2 = (p - 1) / (mid << 1);
                Mint MWn = Moperator1.pow(Moperator2);
                #pragma omp for schedule(static)
                for(int j = 0; j < limit; j += (mid << 1)) {
                    Mint Mw = 1;
                    for(int k = 0; k < mid; k++, Mw = Mw * MWn) {                   
                        Mint x = a[j + k];
                        Mint y = Mw * a[j + k + mid];                    
                        a[j + k] = (x + y),
                        a[j + k + mid] = (x - y);
                    }
                }
            }

        }

    }

    MPI_Bcast(a, limit, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);

    delete[] local_a;
}

void ntt_Montgomery_Mint_mpi_openmp(Mint *a,Mint *b,Mint *ab,int *r,int n){
    MPI_Datatype mpi_montgomerymodint32_type;
    create_montgomerymodint32_type(&mpi_montgomerymodint32_type);
    int L=0,i=0;
    int limit=1;
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
	for(i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }

    int my_rank, processCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCounts);
    ntt_Montgomery_Mint_mpi_openmp(a,r,limit, 1, my_rank, processCounts, mpi_montgomerymodint32_type);
    ntt_Montgomery_Mint_mpi_openmp(b,r,limit, 1, my_rank, processCounts, mpi_montgomerymodint32_type);
    // printf("%u\n",a[0].getv());
	for(i = 0; i < limit; i++){
        ab[i]= a[i] * b[i];
        // std::cout<<ab[i]<<' ';
    } 
	ntt_Montgomery_Mint_mpi_openmp(ab,r,limit, -1, my_rank, processCounts, mpi_montgomerymodint32_type);
    Mint Mlimit = limit;
    Mint invn = (1 / Mlimit);
    // printf("%u %u\n",Mlimit.getv(),invn.getv());
    for(i = 0; i < 2*n; i++){
        ab[i] = (ab[i] * invn);
        // std::cout<<ab[i]<<' ';
    } 
    
}

void ntt_dit_x4_mpi(int *a,int *dw,int imag,int limit,int level,int my_rank, int threadCounts, int initial_w2[][4]){
    int one = 1;
    int logn = __builtin_ctz(limit);
    
    if (my_rank==0&&(logn&1)) {
        for (int i = 0; i < limit/2; i++) {
            int x = a[i], y = a[i + limit/2];
            a[i] = (1LL*x + y)% p;
            a[i + limit/2] = ((1LL*x - y)%p + p)%p;
        }
    }
    int log_max=(logn & ~1)-(__builtin_ctz(threadCounts));
    int SUBARRAY_SIZE = limit / threadCounts;
    int *local_a = new int[SUBARRAY_SIZE];
    if (my_rank==0) {
        for (int e = logn & ~1; e > log_max; e -= 2) {
            const int m = 1 << e, m4 = m >> 2;
            int w2 = 1;
            for (int i = 0; i < limit; i += m) {
                const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2%p;
                for (int j = i; j < i + m4; ++j) {
                    int a0 = 1LL*a[j + m4 * 0] * one%p, a1 = 1LL*a[j + m4 * 1] * w2%p;
                    int a2 = 1LL*a[j + m4 * 2] * w1%p, a3 = 1LL*a[j + m4 * 3] * w3%p;
                    int t02p = (1LL*a0 + a2)%p, t13p = (1LL*a1 + a3)%p;
                    int t02m = ((1LL*a0 - a2)%p + p)%p, t13m = ((1LL*a1 - a3)%p +p)%p * imag%p;
                    a[j + m4 * 0] = (1LL*t02p + t13p)%p;
                    a[j + m4 * 1] = ((1LL*t02p - t13p)%p+p)%p;
                    a[j + m4 * 2] = (1LL*t02m + t13m)%p;
                    a[j + m4 * 3] = ((1LL*t02m - t13m)%p+p)%p;
                }
                w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
            }
        }
    }

    
    MPI_Scatter(a, SUBARRAY_SIZE, MPI_INT, local_a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    for (int e = log_max; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        int w2 = initial_w2[e][my_rank];
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2%p;
            for (int j = i; j < i + m4; ++j) {
                int a0 = 1LL*local_a[j + m4 * 0] * one%p, a1 = 1LL*local_a[j + m4 * 1] * w2%p;
                int a2 = 1LL*local_a[j + m4 * 2] * w1%p, a3 = 1LL*local_a[j + m4 * 3] * w3%p;
                int t02p = (1LL*a0 + a2)%p, t13p = (1LL*a1 + a3)%p;
                int t02m = ((1LL*a0 - a2)%p + p)%p, t13m = ((1LL*a1 - a3)%p +p)%p * imag%p;
                local_a[j + m4 * 0] = (1LL*t02p + t13p)%p;
                local_a[j + m4 * 1] = ((1LL*t02p - t13p)%p+p)%p;
                local_a[j + m4 * 2] = (1LL*t02m + t13m)%p;
                local_a[j + m4 * 3] = ((1LL*t02m - t13m)%p+p)%p;
            }
            w2 = 1LL*w2*dw[__builtin_ctz(~((i+my_rank*SUBARRAY_SIZE) >> e))]%p;
        }
    }
    MPI_Gather(local_a, SUBARRAY_SIZE, MPI_INT, a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(a, limit, MPI_INT, 0, MPI_COMM_WORLD);

    delete[] local_a;
}
void ntt_dif_x4_mpi(int *a,int *dw,int imag,int limit,int level, int my_rank, int threadCounts, int initial_w2[][4]){
    int one = 1;
    int logn = __builtin_ctz(limit);

    int logmax = (logn & ~1);
    int SUBARRAY_SIZE = limit / threadCounts;
    int *local_a = new int[SUBARRAY_SIZE];

    MPI_Scatter(a, SUBARRAY_SIZE, MPI_INT, local_a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    for (int e = 2; e < logmax; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        int w2 = initial_w2[e][my_rank];
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2 %p;
            for (int j = i; j < i + m4; ++j) {
                int a0 = local_a[j + m4 * 0], a1 = local_a[j + m4 * 1];
                int a2 = local_a[j + m4 * 2], a3 = local_a[j + m4 * 3];
                int t01p = (1LL*a0 + a1)%p, t23p = (1LL*a2 + a3)%p;
                int t01m = ((1LL*a0 - a1)%p + p)%p, t23m = ((1LL*a2 - a3)%p +p)%p * imag%p;
                local_a[j + m4 * 0] = (1LL*t01p + t23p)%p*one%p;
                local_a[j + m4 * 2] = ((1LL*t01p - t23p)%p+p)%p*w1%p;
                local_a[j + m4 * 1] = (1LL*t01m + t23m)%p*w2%p;
                local_a[j + m4 * 3] = ((1LL*t01m - t23m)%p+p)%p*w3%p;
            }
            w2 = 1LL*w2*dw[__builtin_ctz(~((i+my_rank*SUBARRAY_SIZE) >> e))]%p;
        }
    }
    MPI_Gather(local_a, SUBARRAY_SIZE, MPI_INT, a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    if(my_rank==0){
        for (int e = logmax; e <= (logn & ~1); e += 2) {
            const int m = 1 << e, m4 = m >> 2;
            int w2 = one;
            for (int i = 0; i < limit; i += m) {
                const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2 %p;
                for (int j = i; j < i + m4; ++j) {
                    int a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                    int a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                    int t01p = (1LL*a0 + a1)%p, t23p = (1LL*a2 + a3)%p;
                    int t01m = ((1LL*a0 - a1)%p + p)%p, t23m = ((1LL*a2 - a3)%p +p)%p * imag%p;
                    a[j + m4 * 0] = (1LL*t01p + t23p)%p*one%p;
                    a[j + m4 * 2] = ((1LL*t01p - t23p)%p+p)%p*w1%p;
                    a[j + m4 * 1] = (1LL*t01m + t23m)%p*w2%p;
                    a[j + m4 * 3] = ((1LL*t01m - t23m)%p+p)%p*w3%p;
                }
                w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
            }
        }
        if (logn&1) {
            for (int i = 0; i < limit/2; i++) {
                int x = a[i], y = a[i + limit/2];
                a[i] = (1LL*x + y )% p;
                a[i + limit/2] = ((1LL*x - y)%p + p)%p;
            }
        }
    }
    MPI_Bcast(a, limit, MPI_INT, 0, MPI_COMM_WORLD);

    delete[] local_a;
}
void ntt_dif_x4_mpi(int *a,int *b,int *ab,int *rt,int *irt,int n){
    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    
    while(limit <= 2 * n-1){
        limit <<= 1, L++;
    } 
    int logn = __builtin_ctz(limit);
    int wnn=3;
    wnn=qpow(wnn,(p-1)>>level,p);
    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = 1LL*rt[i - 1] * rt[i - 1] % p;
    irt[0] = qpow(wnn,p-2,p);
    for (int i = 1; i < level; i++) irt[i] = 1LL*irt[i - 1] * irt[i - 1] % p;

    int dw[level - 1];
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL*dw[i - 1] * irt[level - 1 - i]%p * rt[level - 3 - i]%p;//,printf("%d ",dw[i]);
    dw[level - 2] = 1LL*dw[level - 3] * irt[1]%p;
    int imag = rt[level-2];

    int my_rank, threadCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCounts);

    int log_max=(logn & ~1)-(__builtin_ctz(threadCounts));
    int SUBARRAY_SIZE = limit / threadCounts;
    int initial_w2[level][4];
    
    for (int e = log_max; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        int w2 = 1;
        initial_w2[e][0] = 1;
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][1] = w2;
        for (int i = SUBARRAY_SIZE; i < SUBARRAY_SIZE*2; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][2] = w2;
        for (int i = SUBARRAY_SIZE*2; i < SUBARRAY_SIZE*3; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][3] = w2;
    }

    ntt_dit_x4_mpi(a,dw,imag,limit,level,my_rank,threadCounts,initial_w2);
    ntt_dit_x4_mpi(b,dw,imag,limit,level,my_rank,threadCounts,initial_w2);

    if (my_rank == 0) {
        for(int i = 0;i<limit;i++){
            ab[i]=1LL*a[i]*b[i]%p;   
        }
    }

    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL*dw[i - 1] * rt[level - 1 - i]%p * irt[level - 3 - i]%p;
    dw[level - 2] = 1LL*dw[level - 3] * rt[1]%p;
    imag = irt[level-2];

    log_max=(logn & ~1);
    for (int e = 2; e <= log_max; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        int w2 = 1;
        initial_w2[e][0] = 1;
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][1] = w2;
        for (int i = SUBARRAY_SIZE; i < SUBARRAY_SIZE*2; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][2] = w2;
        for (int i = SUBARRAY_SIZE*2; i < SUBARRAY_SIZE*3; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][3] = w2;
    }
    

	ntt_dif_x4_mpi(ab,dw,imag,limit,level,my_rank,threadCounts,initial_w2);
    int invn=qpow(limit,p-2,p);
    for(int i = 0; i < 2 * n; i++){
        ab[i] = (1LL * ab[i] * invn) % p;
        // std::cout<<ab[i] << " ";
    }
}

void ntt_dit_x4_mpi_openmp(int *a,int *dw,int imag,int limit,int level,int my_rank, int threadCounts, int initial_w2[][4]){
    int one = 1;
    int logn = __builtin_ctz(limit);
    
    if (my_rank==0&&(logn&1)) {
        for (int i = 0; i < limit/2; i++) {
            int x = a[i], y = a[i + limit/2];
            a[i] = (1LL*x + y)% p;
            a[i + limit/2] = ((1LL*x - y)%p + p)%p;
        }
    }
    int log_max=(logn & ~1)-(__builtin_ctz(threadCounts));
    int SUBARRAY_SIZE = limit / threadCounts;
    int *local_a = new int[SUBARRAY_SIZE];
    if (my_rank==0) {
        #pragma omp parallel num_threads(thread_count)
        {
            for (int e = logn & ~1; e > log_max; e -= 2) {
                const int m = 1 << e, m4 = m >> 2;
                int w2 = 1;
                #pragma omp for schedule(static)
                for (int i = 0; i < limit; i += m) {
                    const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2%p;
                    for (int j = i; j < i + m4; ++j) {
                        int a0 = 1LL*a[j + m4 * 0] * one%p, a1 = 1LL*a[j + m4 * 1] * w2%p;
                        int a2 = 1LL*a[j + m4 * 2] * w1%p, a3 = 1LL*a[j + m4 * 3] * w3%p;
                        int t02p = (1LL*a0 + a2)%p, t13p = (1LL*a1 + a3)%p;
                        int t02m = ((1LL*a0 - a2)%p + p)%p, t13m = ((1LL*a1 - a3)%p +p)%p * imag%p;
                        a[j + m4 * 0] = (1LL*t02p + t13p)%p;
                        a[j + m4 * 1] = ((1LL*t02p - t13p)%p+p)%p;
                        a[j + m4 * 2] = (1LL*t02m + t13m)%p;
                        a[j + m4 * 3] = ((1LL*t02m - t13m)%p+p)%p;
                    }
                    w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
                }
            }
        }
    }

    
    MPI_Scatter(a, SUBARRAY_SIZE, MPI_INT, local_a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    #pragma omp parallel num_threads(thread_count)
    {   

        for (int e = log_max; e >= 2; e -= 2) {
            const int m = 1 << e, m4 = m >> 2;
            int w2 = initial_w2[e][my_rank];
            #pragma omp for schedule(static)
            for (int i = 0; i < SUBARRAY_SIZE; i += m) {
                const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2%p;
                for (int j = i; j < i + m4; ++j) {
                    int a0 = 1LL*local_a[j + m4 * 0] * one%p, a1 = 1LL*local_a[j + m4 * 1] * w2%p;
                    int a2 = 1LL*local_a[j + m4 * 2] * w1%p, a3 = 1LL*local_a[j + m4 * 3] * w3%p;
                    int t02p = (1LL*a0 + a2)%p, t13p = (1LL*a1 + a3)%p;
                    int t02m = ((1LL*a0 - a2)%p + p)%p, t13m = ((1LL*a1 - a3)%p +p)%p * imag%p;
                    local_a[j + m4 * 0] = (1LL*t02p + t13p)%p;
                    local_a[j + m4 * 1] = ((1LL*t02p - t13p)%p+p)%p;
                    local_a[j + m4 * 2] = (1LL*t02m + t13m)%p;
                    local_a[j + m4 * 3] = ((1LL*t02m - t13m)%p+p)%p;
                }
                w2 = 1LL*w2*dw[__builtin_ctz(~((i+my_rank*SUBARRAY_SIZE) >> e))]%p;
            }
        }
    }

    MPI_Gather(local_a, SUBARRAY_SIZE, MPI_INT, a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(a, limit, MPI_INT, 0, MPI_COMM_WORLD);

    delete[] local_a;
}
void ntt_dif_x4_mpi_openmp(int *a,int *dw,int imag,int limit,int level, int my_rank, int threadCounts, int initial_w2[][4]){
    int one = 1;
    int logn = __builtin_ctz(limit);

    int logmax = (logn & ~1);
    int SUBARRAY_SIZE = limit / threadCounts;
    int *local_a = new int[SUBARRAY_SIZE];

    MPI_Scatter(a, SUBARRAY_SIZE, MPI_INT, local_a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    #pragma omp parallel num_threads(thread_count)
    {   
        for (int e = 2; e < logmax; e += 2) {
            const int m = 1 << e, m4 = m >> 2;
            int w2 = initial_w2[e][my_rank];
            #pragma omp for schedule(static)
            for (int i = 0; i < SUBARRAY_SIZE; i += m) {
                const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2 %p;
                for (int j = i; j < i + m4; ++j) {
                    int a0 = local_a[j + m4 * 0], a1 = local_a[j + m4 * 1];
                    int a2 = local_a[j + m4 * 2], a3 = local_a[j + m4 * 3];
                    int t01p = (1LL*a0 + a1)%p, t23p = (1LL*a2 + a3)%p;
                    int t01m = ((1LL*a0 - a1)%p + p)%p, t23m = ((1LL*a2 - a3)%p +p)%p * imag%p;
                    local_a[j + m4 * 0] = (1LL*t01p + t23p)%p*one%p;
                    local_a[j + m4 * 2] = ((1LL*t01p - t23p)%p+p)%p*w1%p;
                    local_a[j + m4 * 1] = (1LL*t01m + t23m)%p*w2%p;
                    local_a[j + m4 * 3] = ((1LL*t01m - t23m)%p+p)%p*w3%p;
                }
                w2 = 1LL*w2*dw[__builtin_ctz(~((i+my_rank*SUBARRAY_SIZE) >> e))]%p;
            }
        }
    }

    MPI_Gather(local_a, SUBARRAY_SIZE, MPI_INT, a, SUBARRAY_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    if(my_rank==0){
        #pragma omp parallel num_threads(thread_count)
        {
            for (int e = logmax; e <= (logn & ~1); e += 2) {
                const int m = 1 << e, m4 = m >> 2;
                int w2 = one;
                #pragma omp for schedule(static)
                for (int i = 0; i < limit; i += m) {
                    const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2 %p;
                    for (int j = i; j < i + m4; ++j) {
                        int a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                        int a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                        int t01p = (1LL*a0 + a1)%p, t23p = (1LL*a2 + a3)%p;
                        int t01m = ((1LL*a0 - a1)%p + p)%p, t23m = ((1LL*a2 - a3)%p +p)%p * imag%p;
                        a[j + m4 * 0] = (1LL*t01p + t23p)%p*one%p;
                        a[j + m4 * 2] = ((1LL*t01p - t23p)%p+p)%p*w1%p;
                        a[j + m4 * 1] = (1LL*t01m + t23m)%p*w2%p;
                        a[j + m4 * 3] = ((1LL*t01m - t23m)%p+p)%p*w3%p;
                    }
                    w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
                }
            }
        }

        if (logn&1) {
            for (int i = 0; i < limit/2; i++) {
                int x = a[i], y = a[i + limit/2];
                a[i] = (1LL*x + y )% p;
                a[i + limit/2] = ((1LL*x - y)%p + p)%p;
            }
        }
    }
    MPI_Bcast(a, limit, MPI_INT, 0, MPI_COMM_WORLD);

    delete[] local_a;
}
void ntt_dif_x4_mpi_openmp(int *a,int *b,int *ab,int *rt,int *irt,int n){
    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    
    while(limit <= 2 * n-1){
        limit <<= 1, L++;
    } 
    int logn = __builtin_ctz(limit);
    int wnn=3;
    wnn=qpow(wnn,(p-1)>>level,p);
    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = 1LL*rt[i - 1] * rt[i - 1] % p;
    irt[0] = qpow(wnn,p-2,p);
    for (int i = 1; i < level; i++) irt[i] = 1LL*irt[i - 1] * irt[i - 1] % p;

    int dw[level - 1];
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL*dw[i - 1] * irt[level - 1 - i]%p * rt[level - 3 - i]%p;//,printf("%d ",dw[i]);
    dw[level - 2] = 1LL*dw[level - 3] * irt[1]%p;
    int imag = rt[level-2];

    int my_rank, threadCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCounts);

    int log_max=(logn & ~1)-(__builtin_ctz(threadCounts));
    int SUBARRAY_SIZE = limit / threadCounts;
    int initial_w2[level][4];
    
    for (int e = log_max; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        int w2 = 1;
        initial_w2[e][0] = 1;
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][1] = w2;
        for (int i = SUBARRAY_SIZE; i < SUBARRAY_SIZE*2; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][2] = w2;
        for (int i = SUBARRAY_SIZE*2; i < SUBARRAY_SIZE*3; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][3] = w2;
    }

    ntt_dit_x4_mpi_openmp(a,dw,imag,limit,level,my_rank,threadCounts,initial_w2);
    ntt_dit_x4_mpi_openmp(b,dw,imag,limit,level,my_rank,threadCounts,initial_w2);

    if (my_rank == 0) {
        for(int i = 0;i<limit;i++){
            ab[i]=1LL*a[i]*b[i]%p;   
        }
    }

    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL*dw[i - 1] * rt[level - 1 - i]%p * irt[level - 3 - i]%p;
    dw[level - 2] = 1LL*dw[level - 3] * rt[1]%p;
    imag = irt[level-2];

    log_max=(logn & ~1);
    for (int e = 2; e <= log_max; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        int w2 = 1;
        initial_w2[e][0] = 1;
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][1] = w2;
        for (int i = SUBARRAY_SIZE; i < SUBARRAY_SIZE*2; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][2] = w2;
        for (int i = SUBARRAY_SIZE*2; i < SUBARRAY_SIZE*3; i += m) {
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        initial_w2[e][3] = w2;
    }
    

	ntt_dif_x4_mpi_openmp(ab,dw,imag,limit,level,my_rank,threadCounts,initial_w2);
    int invn=qpow(limit,p-2,p);
    for(int i = 0; i < 2 * n; i++){
        ab[i] = (1LL * ab[i] * invn) % p;
        // std::cout<<ab[i] << " ";
    }
}

void ntt_dit_x4_Mint_mpi(Mint *a,Mint *dw,Mint imag,int limit,int level,int my_rank, int processCounts, Mint initial_w2[][4], MPI_Datatype mpi_montgomerymodint32_type){
    Mint one = 1;
    int logn = __builtin_ctz(limit);
    
    if (my_rank==0&&(logn&1)) {
        for (int i = 0; i < limit/2; i++) {
            Mint x = a[i], y = a[i + limit/2];
            a[i] = x + y;
            a[i + limit/2] = x - y;
        }
    }
    int log_max=(logn & ~1)-(__builtin_ctz(processCounts));
    int SUBARRAY_SIZE = limit / processCounts;
    Mint *local_a = new Mint[SUBARRAY_SIZE];
    if (my_rank==0) {
        for (int e = logn & ~1; e > log_max; e -= 2) {
            const int m = 1 << e, m4 = m >> 2;
            Mint w2 = one;
            for (int i = 0; i < limit; i += m) {
                const Mint w1 = w2 * w2, w3 = w1 * w2;
                for (int j = i; j < i + m4; ++j) {
                    Mint a0 = a[j + m4 * 0] * one, a1 = a[j + m4 * 1] * w2;
                    Mint a2 = a[j + m4 * 2] * w1, a3 = a[j + m4 * 3] * w3;
                    Mint t02p = a0 + a2, t13p = a1 + a3;
                    Mint t02m = a0 - a2, t13m = (a1 - a3) * imag;
                    a[j + m4 * 0] = t02p + t13p;
                    a[j + m4 * 1] = t02p - t13p;
                    a[j + m4 * 2] = t02m + t13m;
                    a[j + m4 * 3] = t02m - t13m;
                }
                w2 = w2*dw[__builtin_ctz(~(i >> e))];
            }
        }
    }

    
    MPI_Scatter(a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);
    for (int e = log_max; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = initial_w2[e][my_rank];
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            const Mint w1 = w2 * w2, w3 = w1 * w2;
            for (int j = i; j < i + m4; ++j) {
                Mint a0 = local_a[j + m4 * 0] * one, a1 = local_a[j + m4 * 1] * w2;
                Mint a2 = local_a[j + m4 * 2] * w1, a3 = local_a[j + m4 * 3] * w3;
                Mint t02p = a0 + a2, t13p = a1 + a3;
                Mint t02m = a0 - a2, t13m = (a1 - a3) * imag;
                local_a[j + m4 * 0] = t02p + t13p;
                local_a[j + m4 * 1] = t02p - t13p;
                local_a[j + m4 * 2] = t02m + t13m;
                local_a[j + m4 * 3] = t02m - t13m;
            }
            w2 = w2*dw[__builtin_ctz(~((i+my_rank*SUBARRAY_SIZE) >> e))];
        }
    }
    MPI_Gather(local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);

    MPI_Bcast(a, limit, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);

    delete[] local_a;
}
void ntt_dif_x4_Mint_mpi(Mint *a,Mint *dw,Mint imag,int limit,int level, int my_rank, int processCounts, Mint initial_w2[][4], MPI_Datatype mpi_montgomerymodint32_type){
    Mint one = 1;
    int logn = __builtin_ctz(limit);

    int logmax = (logn & ~1);
    int SUBARRAY_SIZE = limit / processCounts;
    Mint *local_a = new Mint[SUBARRAY_SIZE];

    MPI_Scatter(a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);
    for (int e = 2; e < logmax; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = initial_w2[e][my_rank];
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            const Mint w1 = w2 * w2, w3 = w1 * w2;
            for (int j = i; j < i + m4; ++j) {
                Mint a0 = local_a[j + m4 * 0], a1 = local_a[j + m4 * 1];
                Mint a2 = local_a[j + m4 * 2], a3 = local_a[j + m4 * 3];
                Mint t01p = a0 + a1, t23p = a2 + a3;
                Mint t01m = a0 - a1, t23m = (a2 - a3) * imag;
                local_a[j + m4 * 0] = (t01p + t23p)*one;
                local_a[j + m4 * 2] = (t01p - t23p)*w1;
                local_a[j + m4 * 1] = (t01m + t23m)*w2;
                local_a[j + m4 * 3] = (t01m - t23m)*w3;
            }
            w2 = w2*dw[__builtin_ctz(~((i+my_rank*SUBARRAY_SIZE) >> e))];
        }
    }
    MPI_Gather(local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);
    if(my_rank==0){
        for (int e = logmax; e <= (logn & ~1); e += 2) {
            const int m = 1 << e, m4 = m >> 2;
            Mint w2 = one;
            for (int i = 0; i < limit; i += m) {
                const Mint w1 = w2 * w2, w3 = w1 * w2;
                for (int j = i; j < i + m4; ++j) {
                    Mint a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                    Mint a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                    Mint t01p = a0 + a1, t23p = a2 + a3;
                    Mint t01m = a0 - a1, t23m = (a2 - a3) * imag;
                    a[j + m4 * 0] = (t01p + t23p)*one;
                    a[j + m4 * 2] = (t01p - t23p)*w1;
                    a[j + m4 * 1] = (t01m + t23m)*w2;
                    a[j + m4 * 3] = (t01m - t23m)*w3;
                }
                w2 = w2*dw[__builtin_ctz(~(i >> e))];
            }
        }
        if (logn&1) {
            
            for (int i = 0; i < limit/2; i++) {
                Mint x = a[i], y = a[i + limit/2];
                a[i] = x + y;
                a[i + limit/2] = x - y;
            }
           
        }
    }
    MPI_Bcast(a, limit, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);

    delete[] local_a;
}
void ntt_dif_x4_Mint_mpi(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n){
    MPI_Datatype mpi_montgomerymodint32_type;
    create_montgomerymodint32_type(&mpi_montgomerymodint32_type);

    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    
    while(limit <= 2 * n-1){
        limit <<= 1, L++;
    } 
    int logn = __builtin_ctz(limit);

    Mint wnn=3;
    wnn=wnn.pow((p-1)>>level);

    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = rt[i - 1] * rt[i - 1];
    irt[0] = 1/wnn;
    for (int i = 1; i < level; i++) irt[i] = irt[i - 1] * irt[i - 1];

    Mint dw[level - 1];
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * irt[level - 1 - i] * rt[level - 3 - i];//,printf("%d ",dw[i]);
    dw[level - 2] = dw[level - 3] * irt[1];
    Mint imag = rt[level-2];

    int my_rank, processCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCounts);

    int log_max=(logn & ~1)-(__builtin_ctz(processCounts));
    int SUBARRAY_SIZE = limit / processCounts;
    Mint initial_w2[level][4];
    Mint one=1;
    for (int e = log_max; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = one;
        initial_w2[e][0] = one;
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][1] = w2;
        for (int i = SUBARRAY_SIZE; i < SUBARRAY_SIZE*2; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][2] = w2;
        for (int i = SUBARRAY_SIZE*2; i < SUBARRAY_SIZE*3; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][3] = w2;
    }

    ntt_dit_x4_Mint_mpi(a,dw,imag,limit,level,my_rank,processCounts,initial_w2,mpi_montgomerymodint32_type);
    ntt_dit_x4_Mint_mpi(b,dw,imag,limit,level,my_rank,processCounts,initial_w2,mpi_montgomerymodint32_type);

    if (my_rank == 0) {
        for(int i = 0;i<limit;i++){
            ab[i]=a[i]*b[i];   
        }
    }

    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * rt[level - 1 - i] * irt[level - 3 - i];
    dw[level - 2] = dw[level - 3] * rt[1];
    imag = irt[level-2];

    log_max=(logn & ~1);
    for (int e = 2; e <= log_max; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = one;
        initial_w2[e][0] = one;
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][1] = w2;
        for (int i = SUBARRAY_SIZE; i < SUBARRAY_SIZE*2; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][2] = w2;
        for (int i = SUBARRAY_SIZE*2; i < SUBARRAY_SIZE*3; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][3] = w2;
    }
    

	ntt_dif_x4_Mint_mpi(ab,dw,imag,limit,level,my_rank,processCounts,initial_w2,mpi_montgomerymodint32_type);

    Mint Limit=limit;
    Mint invn=1/Limit;
    for(int i = 0; i < 2 * n; i++){
        ab[i] = ab[i] * invn;
        // std::cout << ab[i] << " ";
    }
    
}
void ntt_dit_x4_Mint_mpi_openmp(Mint *a,Mint *dw,Mint imag,int limit,int level,int my_rank, int processCounts, Mint initial_w2[][4], MPI_Datatype mpi_montgomerymodint32_type){
    Mint one = 1;
    int logn = __builtin_ctz(limit);
    
    if (my_rank==0&&(logn&1)) {
        for (int i = 0; i < limit/2; i++) {
            Mint x = a[i], y = a[i + limit/2];
            a[i] = x + y;
            a[i + limit/2] = x - y;
        }
    }
    int log_max=(logn & ~1)-(__builtin_ctz(processCounts));
    int SUBARRAY_SIZE = limit / processCounts;
    Mint *local_a = new Mint[SUBARRAY_SIZE];
    if (my_rank==0) {
        #pragma omp parallel num_threads(thread_count)
        {
            for (int e = logn & ~1; e > log_max; e -= 2) {
                const int m = 1 << e, m4 = m >> 2;
                Mint w2 = one;
                #pragma omp for schedule(static)
                for (int i = 0; i < limit; i += m) {
                    const Mint w1 = w2 * w2, w3 = w1 * w2;
                    for (int j = i; j < i + m4; ++j) {
                        Mint a0 = a[j + m4 * 0] * one, a1 = a[j + m4 * 1] * w2;
                        Mint a2 = a[j + m4 * 2] * w1, a3 = a[j + m4 * 3] * w3;
                        Mint t02p = a0 + a2, t13p = a1 + a3;
                        Mint t02m = a0 - a2, t13m = (a1 - a3) * imag;
                        a[j + m4 * 0] = t02p + t13p;
                        a[j + m4 * 1] = t02p - t13p;
                        a[j + m4 * 2] = t02m + t13m;
                        a[j + m4 * 3] = t02m - t13m;
                    }
                    w2 = w2*dw[__builtin_ctz(~(i >> e))];
                }
            }
        }

    }

    
    MPI_Scatter(a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);
    #pragma omp parallel num_threads(thread_count)
    {
        for (int e = log_max; e >= 2; e -= 2) {
            const int m = 1 << e, m4 = m >> 2;
            Mint w2 = initial_w2[e][my_rank];
            #pragma omp for schedule(static)
            for (int i = 0; i < SUBARRAY_SIZE; i += m) {
                const Mint w1 = w2 * w2, w3 = w1 * w2;
                for (int j = i; j < i + m4; ++j) {
                    Mint a0 = local_a[j + m4 * 0] * one, a1 = local_a[j + m4 * 1] * w2;
                    Mint a2 = local_a[j + m4 * 2] * w1, a3 = local_a[j + m4 * 3] * w3;
                    Mint t02p = a0 + a2, t13p = a1 + a3;
                    Mint t02m = a0 - a2, t13m = (a1 - a3) * imag;
                    local_a[j + m4 * 0] = t02p + t13p;
                    local_a[j + m4 * 1] = t02p - t13p;
                    local_a[j + m4 * 2] = t02m + t13m;
                    local_a[j + m4 * 3] = t02m - t13m;
                }
                w2 = w2*dw[__builtin_ctz(~((i+my_rank*SUBARRAY_SIZE) >> e))];
            }
        }
    }

    MPI_Gather(local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);

    MPI_Bcast(a, limit, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);

    delete[] local_a;
}
void ntt_dif_x4_Mint_mpi_openmp(Mint *a,Mint *dw,Mint imag,int limit,int level, int my_rank, int processCounts, Mint initial_w2[][4], MPI_Datatype mpi_montgomerymodint32_type){
    Mint one = 1;
    int logn = __builtin_ctz(limit);

    int logmax = (logn & ~1);
    int SUBARRAY_SIZE = limit / processCounts;
    Mint *local_a = new Mint[SUBARRAY_SIZE];

    MPI_Scatter(a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);
    #pragma omp parallel num_threads(thread_count)
    {
        for (int e = 2; e < logmax; e += 2) {
            const int m = 1 << e, m4 = m >> 2;
            Mint w2 = initial_w2[e][my_rank];
            for (int i = 0; i < SUBARRAY_SIZE; i += m) {
                const Mint w1 = w2 * w2, w3 = w1 * w2;
                for (int j = i; j < i + m4; ++j) {
                    Mint a0 = local_a[j + m4 * 0], a1 = local_a[j + m4 * 1];
                    Mint a2 = local_a[j + m4 * 2], a3 = local_a[j + m4 * 3];
                    Mint t01p = a0 + a1, t23p = a2 + a3;
                    Mint t01m = a0 - a1, t23m = (a2 - a3) * imag;
                    local_a[j + m4 * 0] = (t01p + t23p)*one;
                    local_a[j + m4 * 2] = (t01p - t23p)*w1;
                    local_a[j + m4 * 1] = (t01m + t23m)*w2;
                    local_a[j + m4 * 3] = (t01m - t23m)*w3;
                }
                w2 = w2*dw[__builtin_ctz(~((i+my_rank*SUBARRAY_SIZE) >> e))];
            }
        }
    }

    MPI_Gather(local_a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, a, SUBARRAY_SIZE, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);
    if(my_rank==0){
        #pragma omp parallel num_threads(thread_count)
        {

            for (int e = logmax; e <= (logn & ~1); e += 2) {
                const int m = 1 << e, m4 = m >> 2;
                Mint w2 = one;
                #pragma omp for schedule(static)
                for (int i = 0; i < limit; i += m) {
                    const Mint w1 = w2 * w2, w3 = w1 * w2;
                    for (int j = i; j < i + m4; ++j) {
                        Mint a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                        Mint a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                        Mint t01p = a0 + a1, t23p = a2 + a3;
                        Mint t01m = a0 - a1, t23m = (a2 - a3) * imag;
                        a[j + m4 * 0] = (t01p + t23p)*one;
                        a[j + m4 * 2] = (t01p - t23p)*w1;
                        a[j + m4 * 1] = (t01m + t23m)*w2;
                        a[j + m4 * 3] = (t01m - t23m)*w3;
                    }
                    w2 = w2*dw[__builtin_ctz(~(i >> e))];
                }
            }

        }


        if (logn&1) {
            
            for (int i = 0; i < limit/2; i++) {
                Mint x = a[i], y = a[i + limit/2];
                a[i] = x + y;
                a[i + limit/2] = x - y;
            }
           
        }
    }
    MPI_Bcast(a, limit, mpi_montgomerymodint32_type, 0, MPI_COMM_WORLD);

    delete[] local_a;
}
void ntt_dif_x4_Mint_mpi_openmp(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n){
    MPI_Datatype mpi_montgomerymodint32_type;
    create_montgomerymodint32_type(&mpi_montgomerymodint32_type);

    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    
    while(limit <= 2 * n-1){
        limit <<= 1, L++;
    } 
    int logn = __builtin_ctz(limit);

    Mint wnn=3;
    wnn=wnn.pow((p-1)>>level);

    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = rt[i - 1] * rt[i - 1];
    irt[0] = 1/wnn;
    for (int i = 1; i < level; i++) irt[i] = irt[i - 1] * irt[i - 1];

    Mint dw[level - 1];
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * irt[level - 1 - i] * rt[level - 3 - i];//,printf("%d ",dw[i]);
    dw[level - 2] = dw[level - 3] * irt[1];
    Mint imag = rt[level-2];

    int my_rank, processCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCounts);

    int log_max=(logn & ~1)-(__builtin_ctz(processCounts));
    int SUBARRAY_SIZE = limit / processCounts;
    Mint initial_w2[level][4];
    Mint one=1;
    for (int e = log_max; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = one;
        initial_w2[e][0] = one;
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][1] = w2;
        for (int i = SUBARRAY_SIZE; i < SUBARRAY_SIZE*2; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][2] = w2;
        for (int i = SUBARRAY_SIZE*2; i < SUBARRAY_SIZE*3; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][3] = w2;
    }

    ntt_dit_x4_Mint_mpi_openmp(a,dw,imag,limit,level,my_rank,processCounts,initial_w2,mpi_montgomerymodint32_type);
    ntt_dit_x4_Mint_mpi_openmp(b,dw,imag,limit,level,my_rank,processCounts,initial_w2,mpi_montgomerymodint32_type);

    if (my_rank == 0) {
        for(int i = 0;i<limit;i++){
            ab[i]=a[i]*b[i];   
        }
    }

    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * rt[level - 1 - i] * irt[level - 3 - i];
    dw[level - 2] = dw[level - 3] * rt[1];
    imag = irt[level-2];

    log_max=(logn & ~1);
    for (int e = 2; e <= log_max; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = one;
        initial_w2[e][0] = one;
        for (int i = 0; i < SUBARRAY_SIZE; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][1] = w2;
        for (int i = SUBARRAY_SIZE; i < SUBARRAY_SIZE*2; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][2] = w2;
        for (int i = SUBARRAY_SIZE*2; i < SUBARRAY_SIZE*3; i += m) {
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        initial_w2[e][3] = w2;
    }
    

	ntt_dif_x4_Mint_mpi_openmp(ab,dw,imag,limit,level,my_rank,processCounts,initial_w2,mpi_montgomerymodint32_type);

    Mint Limit=limit;
    Mint invn=1/Limit;
    for(int i = 0; i < 2 * n; i++){
        ab[i] = ab[i] * invn;
    }
    
}


void ntt_dit_x4(int *a,int *rt,int *irt,int limit,int level){
    int one = 1, imag = rt[level-2];
    int logn = __builtin_ctz(limit);
    int dw[level - 1];
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL*dw[i - 1] * irt[level - 1 - i]%p * rt[level - 3 - i]%p;//,printf("%d ",dw[i]);
    dw[level - 2] = 1LL*dw[level - 3] * irt[1]%p;
    if (logn & 1) {
        for (int i = 0; i < limit/2; i++) {
            int x = a[i], y = a[i + limit/2];
            a[i] = (1LL*x + y)% p;
            a[i + limit/2] = ((1LL*x - y)%p + p)%p;
        }
    }
    for (int e = logn & ~1; e >= 2; e -= 2) {
        // std::cout<<a[15]<<'\n';
        const int m = 1 << e, m4 = m >> 2;
        int w2 = 1;
        for (int i = 0; i < limit; i += m) {
            const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2%p;
            for (int j = i; j < i + m4; ++j) {
                int a0 = 1LL*a[j + m4 * 0] * one%p, a1 = 1LL*a[j + m4 * 1] * w2%p;
                int a2 = 1LL*a[j + m4 * 2] * w1%p, a3 = 1LL*a[j + m4 * 3] * w3%p;
                int t02p = (1LL*a0 + a2)%p, t13p = (1LL*a1 + a3)%p;
                int t02m = ((1LL*a0 - a2)%p + p)%p, t13m = ((1LL*a1 - a3)%p +p)%p * imag%p;
                a[j + m4 * 0] = (1LL*t02p + t13p)%p;
                a[j + m4 * 1] = ((1LL*t02p - t13p)%p+p)%p;
                a[j + m4 * 2] = (1LL*t02m + t13m)%p;
                a[j + m4 * 3] = ((1LL*t02m - t13m)%p+p)%p;
            }
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
    }
}
void ntt_dif_x4(int *a,int *rt,int *irt,int limit,int level){
    int one = 1, imag = irt[level-2];
    int logn = __builtin_ctz(limit);
    int dw[level - 1];
    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = 1LL*dw[i - 1] * rt[level - 1 - i]%p * irt[level - 3 - i]%p;
    dw[level - 2] = 1LL*dw[level - 3] * rt[1]%p;
    
    for (int e = 2; e <= logn; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        int w2 = one;
        for (int i = 0; i < limit; i += m) {
            const int w1 = 1LL*w2 * w2 %p, w3 = 1LL*w1 * w2 %p;
            for (int j = i; j < i + m4; ++j) {
                int a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                int a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                int t01p = (1LL*a0 + a1)%p, t23p = (1LL*a2 + a3)%p;
                int t01m = ((1LL*a0 - a1)%p + p)%p, t23m = ((1LL*a2 - a3)%p +p)%p * imag%p;
                a[j + m4 * 0] = (1LL*t01p + t23p)%p*one%p;
                a[j + m4 * 2] = ((1LL*t01p - t23p)%p+p)%p*w1%p;
                a[j + m4 * 1] = (1LL*t01m + t23m)%p*w2%p;
                a[j + m4 * 3] = ((1LL*t01m - t23m)%p+p)%p*w3%p;
            }
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))]%p;
        }
        // printf("%d ",w2);
    }
    if (logn & 1) {
        for (int i = 0; i < limit/2; i++) {
            int x = a[i], y = a[i + limit/2];
            a[i] = (1LL*x + y )% p;
            a[i + limit/2] = ((1LL*x - y)%p + p)%p;
        }
    }
}
void ntt_dif_x4(int *a,int *b,int *ab,int *rt,int *irt,int n){
    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    int wnn=3;
    wnn=qpow(wnn,(p-1)>>level,p);
    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = 1LL*rt[i - 1] * rt[i - 1] % p;
    irt[0] = qpow(wnn,p-2,p);
    for (int i = 1; i < level; i++) irt[i] = 1LL*irt[i - 1] * irt[i - 1] % p;
    ntt_dit_x4(a,rt,irt,limit,level);
    ntt_dit_x4(b,rt,irt,limit,level);
	for(int i = 0;i<limit;i++){
        ab[i]=1LL*a[i]*b[i]%p;
    }
    // printf("%d %d %d\n",a[0],b[0],ab[0]);
	ntt_dif_x4(ab,rt,irt,limit,level);
    int invn=qpow(limit,p-2,p);
    for(int i = 0; i < 2 * n; i++){
        ab[i] = (1LL * ab[i] * invn) % p;
        // printf("%d\n",ab[i]);
    } 
}

void ntt_dit_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level){
    int one = 1;
    int logn = __builtin_ctz(limit);
    if (logn & 1) {
        for (int i = 0; i < limit/2; i++) {
            Mint x = a[i], y = a[i + limit/2];
            a[i] = x + y;
            a[i + limit/2] = x - y;
        }
    }
    for (int e = logn & ~1; e >= 2; e -= 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = 1;
        for (int i = 0; i < limit; i += m) {
            Mint w1 = w2 * w2, w3 = w1 * w2;
            for (int j = i; j < i + m4; ++j) {
                Mint a0 = a[j + m4 * 0] * one, a1 = a[j + m4 * 1] * w2;
                Mint a2 = a[j + m4 * 2] * w1, a3 = a[j + m4 * 3] * w3;
                Mint t02p = a0 + a2, t13p = a1 + a3;
                Mint t02m = a0 - a2, t13m = (a1 - a3) * imag;
                a[j + m4 * 0] = t02p + t13p;
                a[j + m4 * 1] = t02p - t13p;
                a[j + m4 * 2] = t02m + t13m;
                a[j + m4 * 3] = t02m - t13m;
            }
            w2 = 1LL*w2*dw[__builtin_ctz(~(i >> e))];
            // printf("%d ",__builtin_ctz(~(i >> e)));
        }
        // puts("");
    }
}
void ntt_dif_x4_Mint(Mint *a,Mint *dw,Mint imag,int limit,int level){
    int one = 1;
    int logn = __builtin_ctz(limit);
    
    for (int e = 2; e <= logn; e += 2) {
        const int m = 1 << e, m4 = m >> 2;
        Mint w2 = one;
        for (int i = 0; i < limit; i += m) {
            Mint w1 = w2 * w2, w3 = w1 * w2;
            for (int j = i; j < i + m4; ++j) {
                Mint a0 = a[j + m4 * 0], a1 = a[j + m4 * 1];
                Mint a2 = a[j + m4 * 2], a3 = a[j + m4 * 3];
                Mint t01p = a0 + a1, t23p = a2 + a3;
                Mint t01m = a0 - a1, t23m = (a2 - a3)* imag;
                a[j + m4 * 0] = (t01p + t23p)*one;
                a[j + m4 * 2] = (t01p - t23p)*w1;
                a[j + m4 * 1] = (t01m + t23m)*w2;
                a[j + m4 * 3] = (t01m - t23m)*w3;
            }
            w2 = w2*dw[__builtin_ctz(~(i >> e))];
        }
        // printf("%d ",w2);
    }
    if (logn & 1) {
        for (int i = 0; i < limit/2; i++) {
            Mint x = a[i], y = a[i + limit/2];
            a[i] = x + y ;
            a[i + limit/2] = x - y;
        }
    }
}
void ntt_dif_x4_Mint(Mint *a,Mint *b,Mint *ab,Mint *rt,Mint *irt,int n){
    int L=0;
    int limit=1;
    int level=__builtin_ctzll(p - 1);
    while(limit <= 2 * n-2){
        limit <<= 1, L++;
    } 
    Mint wnn=3;
    wnn=wnn.pow((p-1)>>level);
    rt[0] = wnn;
    for (int i = 1; i < level; i++) rt[i] = rt[i - 1] * rt[i - 1];
    irt[0] = 1/wnn;
    for (int i = 1; i < level; i++) irt[i] = irt[i - 1] * irt[i - 1];
    Mint dw[level - 1];
    dw[0] = irt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * rt[level - 1 - i] * irt[level - 3 - i];
    dw[level - 2] = dw[level - 3] * rt[1];

    ntt_dit_x4_Mint(a,dw,irt[level-2],limit,level);
    ntt_dit_x4_Mint(b,dw,irt[level-2],limit,level);

	for(int i = 0;i<limit;i++){
        ab[i]=a[i]*b[i];
    }
    dw[0] = rt[level - 3];
    for (int i = 1; i < level - 2; i++)
        dw[i] = dw[i - 1] * irt[level - 1 - i] * rt[level - 3 - i];//,printf("%d ",dw[i]);
    dw[level - 2] = 1LL*dw[level - 3] * irt[1];
	ntt_dif_x4_Mint(ab,dw,rt[level-2],limit,level);
    Mint Limit=limit;
    Mint invn=1/Limit;
    for(int i = 0; i < 2 * n; i++){
        ab[i] = ab[i] * invn;
        // std::cout << ab[i] << " ";
    } 
}

