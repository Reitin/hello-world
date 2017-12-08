#include <string.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc, char *argv[])

    /**********************************************************************/
    /* Two-dim FFT/
    /**********************************************************************/
{

FILE *bin;
fftwf_complex *Yin,*Yout;
fftwf_plan px,py;

int i,j,k,l,m,n,Nch,NN,MM,A1,A2,A3,A4,A5,M;
int adr1,adr2,adr3,adr4,bdr1,bdr2,bdr3,bdr4;
float *spIN,*spWorkR,*spWorkI,spR,spI,*T,scaleM,scaleN;
float fw,df,fend,cf,sf,sum,rNN,*CF,*SF;
char name1[128],name2[128];
//float *spectr1,*spectr2,*spectr3,*spectr4;
//float *spectr5,*spectr6,*spectr7,*spectr8;
strcpy(name1,argv[1]);
strcpy(name2,argv[2]);
printf("Nch=");
scanf("%d",&Nch);
printf("%d\n",Nch);

printf("mp=");
scanf("%d",&m);
printf("%d\n",m);
NN=Nch*m;
MM=m*m;
rNN=(float) NN;
spIN=new float[NN];
spWorkR=new float[NN];
spWorkI=new float[NN];
CF=new float[MM]; SF=new float[MM];
T=new float[m];
scaleM=1.0/(float) m;
scaleN=1.0/(float) Nch;
Yin=(fftwf_complex *) fftwf_malloc(NN* sizeof(fftwf_complex));
Yout=(fftwf_complex *) fftwf_malloc(NN * sizeof(fftwf_complex));

px = fftwf_plan_dft_1d(m,Yin,Yout,FFTW_FORWARD,FFTW_ESTIMATE);
py = fftwf_plan_dft_1d(Nch,Yin,Yout,FFTW_FORWARD,FFTW_ESTIMATE);
df=2.0*3.141592/4768.1695; //rags04ap
//df=2.0*3.141592/5685.0;
//df=2.0*3.141592/12150.0;
fend=-df*(float) m/2;

bin=fopen(name1,"rb");
fread(spIN,4,NN,bin);
fclose(bin);
sum=0.0;
for(i=0;i<NN;i++)
    sum+=spIN[i];
    sum=sum/rNN;
    A1=192; A2=208; A3=592; A4=881; A5=976;
for(j=0;j<m;j++)
    {
     T[j]=(spIN[j*Nch]-0.04167467)*86400.0;
     //T[j]=(spIN[j*Nch]-0.906264)*86400.0;
      //spIN[j*Nch+A1]=spIN[j*Nch+A1+10];
      //spIN[j*Nch+A2]=spIN[j*Nch+A2+10];
      //spIN[j*Nch+A3]=spIN[j*Nch+A3+10];
      //spIN[j*Nch+A4]=spIN[j*Nch+A4+10];
      //spIN[j*Nch+A5]=spIN[j*Nch+A5+10];
printf("%d	%f\n",j,T[j]);

     spIN[j*Nch]=0.0;
     for(i=0;i<Nch;i++)
        {
        Yin[i][0]=spIN[j*Nch+i]-sum; Yin[i][1]=0.0;
        }
fftwf_execute(py);

for(i=0;i<Nch;i++)
     {
      spWorkR[i+j*Nch]=Yout[i][0];
      spWorkI[i+j*Nch]=Yout[i][1];
     }
   }
/**SIN-COS-matrix*********/
for(k=0;k<m/2;k++)
    {
    fw=df*(float) k;
    for(i=0;i<m;i++)
      {
       M=i+m*k;
       CF[M]=cosf(T[i]*fw); SF[M]=sinf(T[i]*fw);
       }
       //printf("k=%d M=%d\n",k,M);
     }  
for(k=m/2;k<m;k++)
      {
        fw=fend+df*(float) (k-m/2);
	for(i=0;i<m;i++)
	{
	M=i+m*k;
        CF[M]=cosf(T[i]*fw); SF[M]=sinf(T[i]*fw);
        }
	//printf("k1=%d M1=%d\n",k,M);
      }
/****END of SIN-COS matrix******/

for(i=0;i<NN;i++)
   spIN[i]=0.0;
for(j=0;j<Nch;j++)
   {
printf("j1=%d\n",j);
    for(i=0;i<m;i++)
       {
        Yin[i][0]=spWorkR[i*Nch+j];
        Yin[i][1]=spWorkI[i*Nch+j];
       }

/****direct spectrum******************/
   for(k=0;k<m/2;k++)
      {
      // fw=df*(float) k;
       spR=0.0; spI=0.0;
//printf("%d	%e\n",k,fw);
       for(i=0;i<m;i++)
          {
	   M=i+k*m;
           cf=CF[M]; sf=SF[M];
           spR+=scaleM*(Yin[i][0]*cf+Yin[i][1]*sf);
           spI+=scaleM*(Yin[i][1]*cf-Yin[i][0]*sf);
          }
           spIN[k*Nch+j]=spR*spR+spI*spI;         
     }
  for(k=m/2;k<m;k++)
      {
       //fw=fend+df*(float) (k-m/2);
       spR=0.0; spI=0.0;
       for(i=0;i<m;i++)
          {
	  M=i+k*m;
           cf=CF[M]; sf=SF[M];
            spR+=scaleM*(Yin[i][0]*cf+Yin[i][1]*sf);
           spI+=scaleM*(Yin[i][1]*cf-Yin[i][0]*sf);
          }
          
           spIN[k*Nch+j]=spR*spR+spI*spI;         
   }
//fftwf_execute(px);
//printf("j2=%d\n",j);
    //  for(i=0;i<m;i++)
      //   spIN[i*Nch+j]=log(0.1+Yout[i][0]*Yout[i][0]+Yout[i][1]*Yout[i][1]);

      }
spIN[0]=spIN[1];
for(j=0;j<m/2;j++)
    {
     for(i=0;i<Nch/2;i++)
         {
          adr1=i+j*Nch;
          bdr1=Nch*m/2+Nch/2+adr1;
          spWorkR[bdr1]=log10f(0.1+spIN[adr1]);
          adr2=i+j*Nch+Nch/2;
          bdr2=Nch*m/2+adr1;
          spWorkR[bdr2]=log10f(0.1+spIN[adr2]); 
          adr3=bdr1;
          bdr3=adr1;
          spWorkR[bdr3]=log10f(0.1+spIN[adr3]);
          adr4=bdr2;
          bdr4=adr2;
          spWorkR[bdr4]=log10f(0.1+spIN[adr4]);
         }
   }
//for(i=401400;i<401500;i++)
//printf("i=%d	Y=%e\n",i,spIN[i]);
bin=fopen(name2,"wb");
fwrite(spWorkR,4,NN,bin);
//fwrite(spIN,4,NN,bin);
fclose(bin);
}

