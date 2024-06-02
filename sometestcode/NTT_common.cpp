#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <vector>

const int p = 998244353; // 假设一个模数p
// 快速幂计算 (base^exp) % mod
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

// NTT 核心函数
void ntt_common(int *a, int *r, int limit, int type, int my_rank, int threadCounts) {
    // thread 0 负责位反转，然后分发给各个线程
    if (my_rank == 0) {
        for (int i = 0; i < limit; i++) {
            if (i < r[i]) {
                std::swap(a[i], a[r[i]]);
            }
        }
    }

    // 动态设置 SUBARRAY_SIZE
    int SUBARRAY_SIZE = limit / threadCounts;
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
    
    // 将小规模数据处理完毕后，在所有 0 thread进行大规模数据处理后广播
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

    // 此时只有 0 thread 是正确的全局变量，进行广播，统一变量
    MPI_Bcast(a, limit, MPI_INT, 0, MPI_COMM_WORLD);

    delete[] local_a;
}

// NTT 入口函数，调用两次 NTT 进行卷积
void ntt_common(int *a, int *b, int *ab, int *r, int n) {
    // 全局变量
    int L = 0, i = 0;
    int limit = 1;
    while (limit <= 2 * n - 2) {
        limit <<= 1;
        L++;
    }

    for (i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }

    int my_rank, threadCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCounts);

    // ntt_common 内部会将 thread 0 初始化的 a、b 分发
    ntt_common(a, r, limit, 1, my_rank, threadCounts); // 计算后，a gather in 0 thread
    ntt_common(b, r, limit, 1, my_rank, threadCounts); // 计算后, b gather in 0 thread

    // thread 0 处理点值
    if (my_rank == 0) {
        for (i = 0; i < limit; i++) {
            ab[i] = 1LL * a[i] * b[i] % p;
        }
    }
    // ntt_common会把 thread 0 处理点值的结果 ab 分发
    ntt_common(ab, r, limit, -1, my_rank, threadCounts); // INTT后, ab gather in 0 thread

    // 对最终结果进行处理，得到正确结果
    int invn = qpow(limit, p - 2, p);
    for (i = 0; i < 2 * n; i++) {
        ab[i] = (1LL * ab[i] * invn) % p;
    }

}

const int maxn = 3e5 + 50;
int n, m;
int a[maxn], b[maxn], ab[maxn], r[maxn];

int main(int argc, char **argv) {
    freopen("1.in","r",stdin);
    std::cin >> n >> m;
    n++, m++;
    for (int i = 0; i < n; ++i) std::cin >> a[i];
    for (int i = 0; i < m; ++i) std::cin >> b[i];

    double elapsed;

    MPI_Init(NULL, NULL);
    auto start_time = MPI_Wtime();
    ntt_common(a, b, ab, r, n);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    auto end_time = MPI_Wtime();
    double loc_elapsed = end_time - start_time;
    MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(my_rank == 0){
        // for (int i = 0; i < n + m; i++) {
        //     std::cout << ab[i] << " ";
        // }
        // std::cout << std::endl;
        printf("Elapsed time = %lf ms\n", elapsed * 1000);
    }
    MPI_Finalize();
    return 0;
}
