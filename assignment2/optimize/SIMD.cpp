#include"SIMD.h"
#include<cstring>
int nn[11]={256,512,1024,2048,4096,8192,16384,32768,65536,131072};
int qq[5]={1409,3329,7681,12289};
int f[400005],g[400005];
int main(){
	// freopen("../data/datain/NTT256_1409.in","r",stdin);
	freopen("1.out","w",stdout);
	for(int i=0;i<=9;++i){
        for(int j=0;j<=3;++j){
			long double ans=0;
            int cnt=100;
			
			std::string str1="../data/datain/NTT";
			std::string str2=std::to_string(nn[i]);
			std::string str3=std::to_string(qq[i]);
			std::string strin=str1+str2+"_"+str3+".in";
			char charArrayin[strin.size() + 1];
			std::copy(strin.begin(), strin.end(), charArrayin);
			charArrayin[strin.size()] = '\0';
			freopen(charArrayin,"r",stdin);
			int n=nn[i],m=n;
			int limit = poly::bit_up(n + m - 2);
			auto F = alocP(limit), G = alocP(limit);
			
			for(int i = 0; i < n; ++i){cin >> f[i];}
			for(int i = 0; i < m; ++i){cin >> f[i];}
			
			for(int k=1;k<=100;++k){
				memcpy(F,f,sizeof(f));
				memcpy(G,g,sizeof(g));
				auto Start=std::chrono::high_resolution_clock::now();
				fill(F + n , F + limit, zero_Z);
				fill(G + m , G + limit, zero_Z);
				poly::dif<false, 1>(F, limit), poly::dif<false, 1>(G, limit), poly::dot(F, limit, G), poly::dit<true, -1>(F, limit);
				auto End=std::chrono::high_resolution_clock::now();
				std::chrono::duration<double,std::ratio<1,1000>>elapsed=End-Start;
				ans+=elapsed.count();
			}
            std::cout<<ans/cnt<<std::endl;
        }
    }
	
	// for(int i = 0; i < n + m -1; ++i){cout << F[i] << ' ';}
	return 0;
}