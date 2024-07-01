#include "fhe/ckks/ckks.h"
#include <cmath>
#include <iostream>
#include <numeric>
#include <immintrin.h>
#include <sys/time.h>
#include <chrono>
#include <pthread.h>  // Don't forget to include this for pthreads
#include <unistd.h>

using namespace hehub;
#define NUM_THREADS 4

struct ThreadArgs {
    int start;
    int end;
    CkksCt* ct_squared;
    const CkksSk* sk;
    CkksCt* local_sum;
};
pthread_barrier_t barr_merge;
pthread_barrier_t barr_merge1;
pthread_barrier_t barr_merge2;
void* thread_func(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    int start = args->start;
    int end = args->end;
    const CkksCt* ct_square = args->ct_squared;
    CkksCt* local_sum = args->local_sum;
    const CkksSk* sk = args->sk;
    for (int i = start; i <= end; i++) {
        // pthread_barrier_wait(&barr_merge);
        if (i == start) {
            local_sum[i] = ct_square[i];
        } else {
            pthread_barrier_wait(&barr_merge);
            local_sum[i] = ckks::add(local_sum[i-1], ct_square[i]);
        }
    }
    // pthread_barrier_wait(&barr_merge);
    return nullptr;
}

void Calc(int n){
    pthread_barrier_init(&barr_merge,NULL,NUM_THREADS);

    int precision_bits = 30;
    auto params = ckks::create_params(n, precision_bits);
   
    CkksSk sk(params);
    RlweKsk relin_key;
    CkksCt ct_squared[102];
    CkksCt local_sum[102];
    relin_key = get_relin_key(sk, params.additional_mod);
    for (int i = 1; i <= 100; ++i) {
        auto pt = ckks::encode(1.0 / i, params);
        auto ct = ckks::encrypt(pt, sk);
        ct_squared[i] = ckks::mult(ct, ct, relin_key);
    }
    pthread_t threads[NUM_THREADS];
    ThreadArgs thread_args[NUM_THREADS];
    double results[NUM_THREADS];
    int chunk_size = 100 / NUM_THREADS;
    int start = 1;
    for (int t = 0; t < NUM_THREADS; ++t) {
        int end = start + chunk_size - 1;
        // printf("%d %d\n",start,end);
        thread_args[t].start = start;
        thread_args[t].end = end > 100 ? 100 : end;
        thread_args[t].ct_squared = ct_squared;
        thread_args[t].sk = &sk;
        thread_args[t].local_sum = local_sum;
        if (pthread_create(&threads[t], nullptr, thread_func, &thread_args[t]) != 0) {
            std::cerr << "Error creating thread " << t << std::endl;
            exit(1);
        }
        start = end + 1;
    }
    double sum=0;
    CkksCt ct_sum;
    for (int t = 0; t < NUM_THREADS; ++t) {
        // puts("HERER----------------");
        if (pthread_join(threads[t], NULL) != 0) {
            std::cerr << "Error joining thread " << t << std::endl;
            exit(1);
        }
        sum += ckks::decode(ckks::decrypt(local_sum[(t+1)*25], sk));
    }
    pthread_barrier_destroy(&barr_merge);
    // std::cout << sum << std::endl;
}

int nn[11] = {4096, 8192, 16384, 32768};

int main() {
    for(int i = 2; i <= 2; ++i){
        long double ans = 0;
        int cnt = 1;
        for(int j = 1; j <= cnt; ++j){
            auto Start = std::chrono::high_resolution_clock::now();
            Calc(nn[i]);
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::ratio<1,1000>> elapsed = End - Start;
            ans += elapsed.count();
        }
        std::cout << ans / cnt << std::endl;
    }
    return 0;
}
