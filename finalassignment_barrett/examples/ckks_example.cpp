#include "fhe/ckks/ckks.h"
#include <cmath>
#include <iostream>
#include <numeric>
#include<immintrin.h>
#include<sys/time.h>
#include<chrono>

using namespace hehub;
int nn[11]={4096,8192,16384,32768};
int main() {
    freopen("1.out","w",stdout);
    for(int i=0;i<=4;++i){
        long double ans=0;
        int cnt=5;
        for(int j=1;j<=5;++j){
            auto Start=std::chrono::high_resolution_clock::now();
            int precision_bits = 30;
            auto params = ckks::create_params(nn[i], precision_bits);
            CkksSk sk(params);
            auto relin_key = get_relin_key(sk, params.additional_mod);
            
            CkksCt ct_sum;
            for (int i = 1; i <= 100; i++) {
                auto pt = ckks::encode(1.0 / i, params);
                auto ct = ckks::encrypt(pt, sk);
                auto ct_squared = ckks::mult(ct, ct, relin_key);

                if (i == 1) {
                    ct_sum = ct_squared;
                } else {
                    ct_sum = ckks::add(ct_sum, ct_squared);
                }
            }

            double sum = ckks::decode(ckks::decrypt(ct_sum, sk));
            auto End=std::chrono::high_resolution_clock::now();
            std::chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
            ans+=elapsed.count();
            
        }
        std::cout<<ans/cnt<<std::endl;
    }
    
    // std::cout << "(" << sum << ", " << M_PI * M_PI / 6 << ")" << std::endl;
}
