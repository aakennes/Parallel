#include"../params.h"
// #include"readwrite.h"
#include"ntt.h"

#include<cstdio>
#include<iostream>
#include<cstring>

const int maxn=1e6+5;

int p=Original_Q;
int n=Original_N;

int countw,wn[maxn];
int r[maxn];
int rt[maxn],irt[maxn];
Mint RT[maxn],IRT[maxn];

int a[maxn],b[maxn],ab[maxn];
Mint A[maxn],B[maxn],AB[maxn];
MMint aa[maxn],bb[maxn],aabb[maxn];

void poly_mul(){
    int aa[maxn],bb[maxn],aabb[maxn];
    memcpy(aa,a,sizeof(a));
    memcpy(bb,b,sizeof(b));
    memcpy(aabb,ab,sizeof(ab));
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            aabb[i+j]+=1LL*aa[i]*bb[j]%p;
            aabb[i+j]%=p;
        }
    } 
    for(int i=0;i<2*n-1;++i){
        std::cout<<aabb[i]<<" ";
    }
}
int nn[11]={256,512,1024,2048,4096,8192,16384,32768,65536,131072};
int qq[5]={1409,3329,7681,12289};
#include<sys/time.h>
void time_test(){
    freopen("1.out","w",stdout);
    for(int i=0;i<=0;++i){
        for(int j=0;j<=0;++j){
            long double ans=0;
            FILE* fp;
            fRead(fp,A,B,nn[i],qq[j]);
            // ntt_common(a,b,ab,r,nn[i]);
            // ntt_dif(a,b,ab,rt,irt,nn[i]);
            // ntt_dif_x4(a,b,ab,rt,irt,nn[i]);
            // ntt_Montgomery(a,b,ab,r,nn[i]);
            // ntt_Montgomery_Mint(A,B,AB,r,nn[i]);
            // ntt_dif_Mint(A,B,AB,RT,IRT,nn[i]);
            // ntt_dif_x4_Mint(A,B,AB,RT,IRT,nn[i]);
            ntt_common_mpi(a,b,ab,r,nn[i]);
            // std::cout<<elapsed.count()<<std::endl;
            std::cout<<ans<<std::endl;
        }
    }
    
}

int main(){
    FILE* nmread=fopen("result.in","r");
    // freopen("result.in","r",stdin);
    int xx,yy;
    // std::cin>>xx>>yy;
    fscanf(nmread,"%d%d",&xx,&yy);
    fclose(nmread);
    // std::cout<<xx<<" "<<yy<<'\n';
    FILE* fp;
    fRead(fp,a,b,nn[xx],qq[yy]);
    // fRead(fp,A,B,nn[xx],qq[yy]);
    // freopen("1.in","r",stdin);
    // std::cin>>n>>n;
    // for(int i=0;i<n;++i)std::cin>>a[i];
    // for(int i=0;i<n;++i)std::cin>>b[i];
    double elapsed;
    MPI_Init(NULL, NULL);
    auto start_time = MPI_Wtime();
    ntt_common_mpi(a,b,ab,r,nn[xx]);
    // ntt_Montgomery_Mint_mpi(A,B,AB,r,nn[xx]);
    // ntt_dif_x4_mpi(a,b,ab,rt,irt,nn[xx]);
    // ntt_dif_x4_Mint_mpi(A,B,AB,RT,IRT,nn[xx]);

    
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
    // ntt_dif_x4(a,b,ab,rt,irt,nn[xx]);
    
    // poly_mul();
    // ntt_common(a,b,ab,r);
    // ntt_dif(a,b,ab,rt,irt);
    
    // for(int i=0;i<nn[0];++i)std::cout<<a[i]<<" ";
    // ntt_Montgomery(a,b,ab,r,nn[0]);
    // for(int i=0;i<nn[0];++i)aa[i]=In((u32)a[i]);
    // for(int i=0;i<nn[0];++i)bb[i]=In((u32)b[i]);
    // std::cout<<aa[0]<<"\n";
    // ntt_Montgomery_MMint(aa,bb,aabb,r,nn[0]);
    // std::cout<<get(aabb[1])<<'\n';
    // for(int i=0;i<nn[0]*2-1;++i)std::cout<<get(aabb[i])<<" ";
    // fWrite(fp,ab,n);
    // puts("a--------------");
    // fRead(fp,A,B,nn[0],qq[0]);
    // std::cout<<A[0].getv()<<'\n';
    // ntt_Montgomery_Mint(A,B,AB,r,nn[0]);
    // std::cout<<AB[0]<<'\n';
    // ntt_dif_Mint(A,B,AB,RT,IRT);
    // ntt_dif_x4_Mint(A,B,AB,RT,IRT);
    // fWrite(fp,AB,n);

    // fclose(fp);
    return 0;
}