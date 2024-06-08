#include"SIMD.h"
#include<cstring>
int nn[11]={256,512,1024,2048,4096,8192,16384,32768,65536,131072};
int qq[5]={1409,3329,7681,12289};
int f[400005],g[400005];
int main(){
	 
	FILE* nmread=fopen("result.in","r");
    int xx,yy;
    // std::cin>>xx>>yy;
    fscanf(nmread,"%d%d",&xx,&yy);
    fclose(nmread);
    // std::cout<<xx<<" "<<yy<<'\n';
    FILE* fp;
	long double ans=0;
	
	std::string str1="../data/datain/NTT";
	std::string str2=std::to_string(nn[xx]);
	std::string str3=std::to_string(qq[yy]);
	std::string strin=str1+str2+"_"+str3+".in";
	char charArrayin[strin.size() + 1];
	std::copy(strin.begin(), strin.end(), charArrayin);
	charArrayin[strin.size()] = '\0';
	FILE* fp=fopen(charArrayin, "r+");

	int n=nn[xx],m=n;
	int limit = poly::bit_up(n + m - 2);
	auto F = alocP(limit), G = alocP(limit);
	
	for(int i = 0; i < n; ++i){fscanf(fp,"%d",&f[i]);}
	for(int i = 0; i < m; ++i){fscanf(fp,"%d",&g[i]);}
	
	memcpy(F,f,sizeof(f));
	memcpy(G,g,sizeof(g));
	fill(F + n , F + limit, zero_Z);
	fill(G + m , G + limit, zero_Z);

	double elapsed;
    MPI_Init(NULL, NULL);
    auto start_time = MPI_Wtime();

	int my_rank, threadCounts;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &threadCounts);

    int log_max=(logn & ~1)-(__builtin_ctz(threadCounts));
    int SUBARRAY_SIZE = limit / threadCounts;

	poly::dif<false, 1>(F, limit), poly::dif<false, 1>(G, limit);

	poly::dot(F, limit, G);
	
	poly::dit<true, -1>(F, limit);


	int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    auto end_time = MPI_Wtime();
    double loc_elapsed = end_time - start_time;
    MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(my_rank == 0){
        // for(int i=0;i<nn[0]*2;++i)std::cout<<ab[i]<<" ";
        FILE* fpwrite=fopen("test.out","a+");
        // std::cout<<loc_elapsed*1000<<'\n';
        fprintf(fpwrite,"%lf\n",loc_elapsed*1000);
    }
    MPI_Finalize();

	return 0;
}