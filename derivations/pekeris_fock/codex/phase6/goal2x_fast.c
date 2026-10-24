#include <stdint.h>
#define P 1004535809ULL
uint32_t modpow(uint32_t a,uint32_t n){uint64_t r=1;while(n){if(n&1)r=r*a%P;a=(uint64_t)a*a%P;n>>=1;}return r;}
void convolution(int n,const uint32_t*x,const uint32_t*y,uint32_t*out,const int*off,const int*aa,const int*bb){
for(int k=0;k<n;k++){uint64_t sum=0;int count=0;for(int j=off[k];j<off[k+1];j++){sum+=(uint64_t)x[aa[j]]*y[bb[j]];if(++count==8){sum%=P;count=0;}}out[k]=sum%P;}}
void inverse(int n,const uint32_t*x,uint32_t*out,const int*off,const int*aa,const int*bb){
out[0]=modpow(x[0],P-2);for(int k=1;k<n;k++){uint64_t sum=0;int count=0;for(int j=off[k];j<off[k+1];j++){if(aa[j]==0)continue;sum+=(uint64_t)x[aa[j]]*out[bb[j]];if(++count==8){sum%=P;count=0;}}out[k]=(P-sum%P)*out[0]%P;}}
