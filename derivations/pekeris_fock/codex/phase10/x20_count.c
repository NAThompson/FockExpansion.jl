/* Small serial point counts for the auxiliary e2 K3. Tables bounded at N<=7000. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
static int p,n,N,rel[3];
static uint16_t *adds,*muls;
static int add0(int a,int b){int r=0,k=1;for(int j=0;j<n;j++){r+=((a%p+b%p)%p)*k;a/=p;b/=p;k*=p;}return r;}
static int mul0(int a,int b){int aa[3]={0},bb[3]={0},c[5]={0};for(int j=0;j<n;j++){aa[j]=a%p;bb[j]=b%p;a/=p;b/=p;}for(int j=0;j<n;j++)for(int k=0;k<n;k++)c[j+k]=(c[j+k]+aa[j]*bb[k])%p;for(int j=2*n-2;j>=n;j--)for(int k=0;k<n;k++)c[j-n+k]=(c[j-n+k]+c[j]*rel[k])%p;int r=0,k=1;for(int j=0;j<n;j++){r+=c[j]*k;k*=p;}return r;}
#define AD(a,b) adds[(size_t)(a)*N+(b)]
#define MU(a,b) muls[(size_t)(a)*N+(b)]
static int neg(int a){return MU(p-1,a);}
static int pow0(int a,int k){int r=1;while(k){if(k&1)r=MU(r,a);a=MU(a,a);k>>=1;}return r;}
int main(int argc,char **argv){if(argc!=8)return 2;p=atoi(argv[1]);n=atoi(argv[2]);rel[0]=atoi(argv[3]);rel[1]=atoi(argv[4]);rel[2]=atoi(argv[5]);int u=atoi(argv[6]),v=atoi(argv[7]);N=1;for(int j=0;j<n;j++)N*=p;if(n>3||N>7000)return 3;size_t bytes=(size_t)N*N*sizeof(uint16_t);adds=malloc(bytes);muls=malloc(bytes);int8_t *chi=calloc(N,sizeof(int8_t));if(!adds||!muls||!chi)return 4;
for(int a=0;a<N;a++)for(int b=0;b<N;b++){AD(a,b)=add0(a,b);MU(a,b)=mul0(a,b);}for(int a=1;a<N;a++)chi[MU(a,a)]=1;for(int a=1;a<N;a++)if(!chi[a])chi[a]=-1;
int A=AD(1,neg(MU(v,v)));long long total=0;int croots=0,kroots=0;
for(int qq=0;qq<=N;qq++){int pp1,pp2,D;if(qq==N){pp1=1;pp2=0;D=1;}else{int q2=MU(qq,qq),q4=MU(q2,q2);pp1=AD(q4,p-1);pp2=MU(2,q2);D=AD(1,q4);}int px=AD(MU(u,pp1),MU(v,pp2));int B=AD(pp2,neg(MU(v,px)));int C=AD(MU(D,D),neg(MU(px,px)));int a2=MU(2,B),a4=MU(A,C);if(!C)croots++;if(!pp1)kroots++;
for(int xx=0;xx<N;xx++){int rhs=MU(xx,AD(MU(xx,AD(xx,a2)),a4));total+=chi[rhs];}}
printf("{\"p\":%d,\"n\":%d,\"N\":%d,\"u\":%d,\"v\":%d,\"char_sum\":%lld,\"C_roots\":%d,\"cross_roots\":%d,\"points\":%lld,\"table_bytes\":%zu}\n",p,n,N,u,v,total,croots,kroots,(long long)(N+1)*(N+1)+total+(long long)N*(croots+kroots),2*bytes);free(adds);free(muls);free(chi);return 0;}
