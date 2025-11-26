#include<stdio.h>
#include<math.h>
#include<complex.h>
#include <string>
#include <fstream>

double pi=acos(-1);
double kx,ky;
double w,w1,w2,w3,w4;

double S;
double h;
double J;
double K;
double G;
double J3;
double phi;
double v0,v1,v2Re,v2Im;

FILE *gap;
FILE *gap0;
FILE *gaph;
FILE *gapphi;
FILE *spectral110;
FILE *spectral12Re0;
FILE *spectral12Im0;
FILE *spectral13Re0;
FILE *spectral13Im0;
FILE *spectral14Re0;
FILE *spectral14Im0;
FILE *G11Re0;
FILE *G11Im0;
FILE *G12Re0;
FILE *G12Im0;
FILE *G13Re0;
FILE *G13Im0;
FILE *G14Re0;
FILE *G14Im0;
FILE *test;
FILE *bstphi;
FILE *bsth;
FILE *bst;
FILE *hcritical;
FILE *bstexistence;
FILE *bstJ;
FILE *bstKG;
FILE *gapJ;
FILE *gapKG;


double Ree1(double phi,double kx,double ky);
double Ime1(double phi,double kx,double ky);
double Ree2(double phi,double kx,double ky);
double Ime2(double phi,double kx,double ky);
double e0(double h);
double spectrum(double h,double phi,double kx,double ky);
double spectrumterm1(double h,double phi,double kx,double ky);
double spectrumterm2(double h,double phi,double kx,double ky);

double A(double h,double phi,double w, double kx, double ky);
double B(double h,double phi,double w, double kx, double ky);
double C11(double h,double phi,double w, double kx, double ky);
double D11(double h,double phi,double w, double kx, double ky);
double C12Re(double h,double phi,double w, double kx, double ky);
double C12Im(double h,double phi,double w, double kx, double ky);
double D12Re(double h,double phi,double w, double kx, double ky);
double D12Im(double h,double phi,double w, double kx, double ky);
double C13Re(double h,double phi,double w, double kx, double ky);
double C13Im(double h,double phi,double w, double kx, double ky);
double D13Re(double h,double phi,double w, double kx, double ky);
double D13Im(double h,double phi,double w, double kx, double ky);
double C14Re(double h,double phi,double w, double kx, double ky);
double C14Im(double h,double phi,double w, double kx, double ky);
double D14Re(double h,double phi,double w, double kx, double ky);
double D14Im(double h,double phi,double w, double kx, double ky);

double spectral11k(double h,double phi,double w,double kx,double ky);
double spectral12kRe(double h,double phi,double w,double kx,double ky);
double spectral12kIm(double h,double phi,double w,double kx,double ky);
double spectral13kRe(double h,double phi,double w,double kx,double ky);
double spectral13kIm(double h,double phi,double w,double kx,double ky);
double spectral14kRe(double h,double phi,double w,double kx,double ky);
double spectral14kIm(double h,double phi,double w,double kx,double ky);
//double complex determinant(double complex G11,double complex G12,double complex G13,double complex G14,double complex G21,double complex G22,double complex G23,double complex G24,double complex G31,double complex G32,double complex G33,double complex G34,double complex G41,double complex G42,double complex G43,double complex G44);

int main()
{
	S=1.;
	J=-0.5;
	K=-5.0;
	G=2.5;
	J3=0.5;
	
//	S=1.;
//	J=1.;
//	K=0;
//	G=0;
//	J3=0;
	
	double J0,K0,G0,a;
	
	double stepskx,stepsky,stepsw;
	double dkx,dky,dw;
	double wmin,wmax;
	int i,j;
	int j1,j2,j3,j4;
	int realcond,bstexists;
	double diff;
	double add11,add33,add13Re,add13Im,add31Re,add31Im;
	double Gap,GapBZx,GapBZy,energy,boundstate;
	double phimin,phimax,stepsphi,dphi;
	double hcrit,hmax,stepsh,dh;
//	double complex det;

	
	stepskx=400.0;
	stepsky=400.0;
	stepsw=2000.0;
	dkx=2*pi/(sqrt(3.0)*stepskx);
	dky=4*pi/(3*stepsky);
	wmin=-10.0;
	wmax=10.0;
	dw=(wmax-wmin)/stepsw;
	stepsphi=12.0;
	phimin=0.0;
	phimax=2*pi/3.0;
	dphi=(phimax-phimin)/stepsphi;
	hmax=15;
	dh=0.5;
//	stepsh=25;
//	dh=hmax/stepsh;
	
	
	double spectral11[(int)stepsw];
	double spectral12Re[(int)stepsw];
	double spectral12Im[(int)stepsw];
	double spectral13Re[(int)stepsw];
	double spectral13Im[(int)stepsw];
	double spectral14Re[(int)stepsw];
	double spectral14Im[(int)stepsw];
	double G11Re[(int)stepsw];
	double G11Im[(int)stepsw];
	double G12Re[(int)stepsw];
	double G12Im[(int)stepsw];
	double G13Re[(int)stepsw];
	double G13Im[(int)stepsw];
	double G14Re[(int)stepsw];
	double G14Im[(int)stepsw];
	double T11[(int)stepsw];
	double T13[(int)stepsw];
	double det1Re[(int)stepsw];
	double det2Re[(int)stepsw];
	double det1Im[(int)stepsw];
	double det2Im[(int)stepsw];
	double PolesRe[(int)stepsw];
	double PolesIm[(int)stepsw];
	double discriminant[(int)stepsw];
	double solution1[(int)stepsw];
	double solution2[(int)stepsw];
	
	std::ofstream discr;
	std::ofstream sol1;
	std::ofstream sol2;
	std::ofstream polesRe;
	std::ofstream polesIm;
	
/*	double complex G11,G12,G13,G14,G21,G22,G23,G24,G31,G32,G33,G34,G41,G42,G43,G44;*/
	hcrit=S*(2*J+K-0.5*G+6*J3+sqrt(K*K-G*K+9*G*G/4.0));
	phi=pi/2.+pi/3.;
	h=2*hcrit;

	
/*	bstexistence=fopen("bstexistence.txt","w");
	
	for(K=-8;K<-5;K+=0.5)
	{
	for(G=4;G<10;G+=0.5)
	{
	printf("K=%lf \t G=%lf \n",K,G);*/
	
	test=fopen("test.txt","w");
	
	//hcritical=fopen("hcritical.txt","w");
	
	//bst=fopen("boundstates.txt","w");
	bsth=fopen("Tboundstates.dat","w");
	//bstphi=fopen("boundstatesphi.txt","w");
	
	//gap=fopen("gap.txt","w");
	gaph=fopen("gaph.dat","w");
	//gapphi=fopen("gapphi.txt","w");
	
		
	//for(phi=phimin;phi<phimax;phi+=dphi)
	//{
		printf("\nphi=%lf \n",phi);	
		
	//	h=hmax;
	//	realcond=1;	
	//	bstexists=0;
		/*for(kx=0.5*dkx;kx<2*pi/sqrt(3.0);kx+=dkx)
			{
				for(ky=(-kx/sqrt(3.0)+0.5*dky);ky<(-kx/sqrt(3.0)+4*pi/3.0);ky+=dky)
				{
					diff=spectrumterm1(h,phi,kx,ky)-spectrumterm2(h,phi,kx,ky);
					
					if(diff<0)
					{
						realcond=0;
						break;
					}
					
				}
				if(realcond==0)
				{
					break;
				}
			}*/
		
	//	if(realcond==1)
	//	realcond=1;	
	//	while(realcond==1)
		for(h=4.;h<4.1;h+=0.2) 
		{
			printf("h=%lf \n",h);
			
			Gap=wmax;
			for(kx=0.5*dkx;kx<2*pi/sqrt(3.0);kx+=dkx)
			{
				for(ky=(-kx/sqrt(3.0)+0.5*dky);ky<(-kx/sqrt(3.0)+4*pi/3.0);ky+=dky)
				{
					energy=spectrum(h,phi,kx,ky);
					if(energy<Gap)
					{
						Gap=energy;
						GapBZx=kx;
						GapBZy=ky;
					}
					fprintf(test,"%lf \t %lf \t %lf \t %lf \n",kx,ky,energy,Gap);	
				}
			}
			fprintf(test,"\n%lf \t %lf \t %lf \n",GapBZx,GapBZy,Gap);
			
			printf("gap=%lf \n",Gap);
			//fprintf(gap,"%lf \t %lf \t %lf\n",phi,h,Gap);
			fprintf(gaph,"%lf \t %lf\n",h,Gap);
			//fprintf(gapphi,"%lf \t %lf\n",phi,Gap);*/
			
			gap0=fopen("gap0.txt","w");
			fprintf(gap0,"%lf \t %d",Gap,0);
			fclose(gap0);
			
			
			for(i=0;i<stepsw;i++)
			{
				spectral11[i]=0;
				spectral12Re[i]=0;
				spectral12Im[i]=0;
				spectral13Re[i]=0;
				spectral13Im[i]=0;
				spectral14Re[i]=0;
				spectral14Im[i]=0;
				//PolesRe[i]=0;
				//PolesIm[i]=0;
			}
			for(i=0;i<stepsw;++i)
			{
				G11Re[i]=0;
				G11Im[i]=0;
				G12Re[i]=0;
				G12Im[i]=0;
				G13Re[i]=0;
				G13Im[i]=0;
				G14Re[i]=0;
				G14Im[i]=0;
			}
			
			//SPECTRAL FUNCTION
			
			for(kx=0.5*dkx;kx<2*pi/sqrt(3.0);kx+=dkx)
			{
				for(ky=(-kx/sqrt(3.0)+0.5*dky);ky<(-kx/sqrt(3.0)+4*pi/3.0);ky+=dky)
				{
					//printf("kx=%lf \t ky=%lf\n",kx,ky);
					w1=sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
					+sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
					+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)-Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
					+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0));
					w2=sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
					-sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
					+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)-Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
					+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0));
					w3=-sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
					+sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
					+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)-Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
					+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0));
					w4=-sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
					-sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
					+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)-Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
					+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0));
					//printf("%lf \t %lf \n",e1me2(kx,ky),e1pe2(kx,ky));
					//printf("w1=%lf \t w2=%lf \t w3=%lf \t w4=%lf \n",w1,w2,w3,w4);
					
					j1=(int)((w1-wmin)/dw);
					j2=(int)((w2-wmin)/dw);
					j3=(int)((w3-wmin)/dw);
					j4=(int)((w4-wmin)/dw);
					//printf("%d \t %d \t %d \t %d \n\n",j1,j2,j3,j4);
					
					if(j1<stepsw)
					{
						spectral11[j1]+=2*pi/(dw*stepskx*stepsky)*spectral11k(h,phi,w1,kx,ky);
						spectral12Re[j1]+=2*pi/(dw*stepskx*stepsky)*spectral12kRe(h,phi,w1,kx,ky);
						spectral12Im[j1]+=2*pi/(dw*stepskx*stepsky)*spectral12kIm(h,phi,w1,kx,ky);
						spectral13Re[j1]+=2*pi/(dw*stepskx*stepsky)*spectral13kRe(h,phi,w1,kx,ky);
						spectral13Im[j1]+=2*pi/(dw*stepskx*stepsky)*spectral13kIm(h,phi,w1,kx,ky);
						spectral14Re[j1]+=2*pi/(dw*stepskx*stepsky)*spectral14kRe(h,phi,w1,kx,ky);
						spectral14Im[j1]+=2*pi/(dw*stepskx*stepsky)*spectral14kIm(h,phi,w1,kx,ky);
					}
					if(j2<stepsw)
					{
						spectral11[j2]+=2*pi/(dw*stepskx*stepsky)*spectral11k(h,phi,w2,kx,ky);
						spectral12Re[j2]+=2*pi/(dw*stepskx*stepsky)*spectral12kRe(h,phi,w2,kx,ky);
						spectral12Im[j2]+=2*pi/(dw*stepskx*stepsky)*spectral12kIm(h,phi,w2,kx,ky);
						spectral13Re[j2]+=2*pi/(dw*stepskx*stepsky)*spectral13kRe(h,phi,w2,kx,ky);
						spectral13Im[j2]+=2*pi/(dw*stepskx*stepsky)*spectral13kIm(h,phi,w2,kx,ky);
						spectral14Re[j2]+=2*pi/(dw*stepskx*stepsky)*spectral14kRe(h,phi,w2,kx,ky);
						spectral14Im[j2]+=2*pi/(dw*stepskx*stepsky)*spectral14kIm(h,phi,w2,kx,ky);
					}
					if(j3<stepsw)
					{
						spectral11[j3]+=2*pi/(dw*stepskx*stepsky)*spectral11k(h,phi,w3,kx,ky);
						spectral12Re[j3]+=2*pi/(dw*stepskx*stepsky)*spectral12kRe(h,phi,w3,kx,ky);
						spectral12Im[j3]+=2*pi/(dw*stepskx*stepsky)*spectral12kIm(h,phi,w3,kx,ky);
						spectral13Re[j3]+=2*pi/(dw*stepskx*stepsky)*spectral13kRe(h,phi,w3,kx,ky);
						spectral13Im[j3]+=2*pi/(dw*stepskx*stepsky)*spectral13kIm(h,phi,w3,kx,ky);
						spectral14Re[j3]+=2*pi/(dw*stepskx*stepsky)*spectral14kRe(h,phi,w3,kx,ky);
						spectral14Im[j3]+=2*pi/(dw*stepskx*stepsky)*spectral14kIm(h,phi,w3,kx,ky);
					}
					if(j4<stepsw)
					{
						spectral11[j4]+=2*pi/(dw*stepskx*stepsky)*spectral11k(h,phi,w4,kx,ky);
						spectral12Re[j4]+=2*pi/(dw*stepskx*stepsky)*spectral12kRe(h,phi,w4,kx,ky);
						spectral12Im[j4]+=2*pi/(dw*stepskx*stepsky)*spectral12kIm(h,phi,w4,kx,ky);
						spectral13Re[j4]+=2*pi/(dw*stepskx*stepsky)*spectral13kRe(h,phi,w4,kx,ky);
						spectral13Im[j4]+=2*pi/(dw*stepskx*stepsky)*spectral13kIm(h,phi,w4,kx,ky);
						spectral14Re[j4]+=2*pi/(dw*stepskx*stepsky)*spectral14kRe(h,phi,w4,kx,ky);
						spectral14Im[j4]+=2*pi/(dw*stepskx*stepsky)*spectral14kIm(h,phi,w4,kx,ky);
					}
				}
			}
			
			spectral110=fopen("spectral110.txt","w");
			spectral12Re0=fopen("spectral12Re0.txt","w");
			spectral12Im0=fopen("spectral12Im0.txt","w");
			spectral13Re0=fopen("spectral13Re0.txt","w");
			spectral13Im0=fopen("spectral13Im0.txt","w");
			spectral14Re0=fopen("spectral14Re0.txt","w");
			spectral14Im0=fopen("spectral14Im0.txt","w");
			for(i=0;i<stepsw;++i)
			{
				w=wmin+i*dw;
				fprintf(spectral110,"%lf \t %lf\n",w,spectral11[i]);
				fprintf(spectral12Re0,"%lf \t %lf\n",w,spectral12Re[i]);
				fprintf(spectral12Im0,"%lf \t %lf\n",w,spectral12Im[i]);
				fprintf(spectral13Re0,"%lf \t %lf\n",w,spectral13Re[i]);
				fprintf(spectral13Im0,"%lf \t %lf\n",w,spectral13Im[i]);
				fprintf(spectral14Re0,"%lf \t %lf\n",w,spectral14Re[i]);
				fprintf(spectral14Im0,"%lf \t %lf\n",w,spectral14Im[i]);
			}
			fclose(spectral110);
			fclose(spectral12Re0);
			fclose(spectral12Im0);
			fclose(spectral13Re0);
			fclose(spectral13Im0);
			fclose(spectral14Re0);
			fclose(spectral14Im0);
			
			//GREEN'S FUNCTION
			
			for(j=0;j<stepsw;++j)
			{
				for(i=0;i<j;++i)
				{
					G11Re[j]+=(dw/(2*pi))*spectral11[i]/(dw*(j-i));
					G12Re[j]+=(dw/(2*pi))*spectral12Re[i]/(dw*(j-i));
					G12Im[j]+=(dw/(2*pi))*spectral12Im[i]/(dw*(j-i));
					G13Re[j]+=(dw/(2*pi))*spectral13Re[i]/(dw*(j-i));
					G13Im[j]+=(dw/(2*pi))*spectral13Im[i]/(dw*(j-i));
					G14Re[j]+=(dw/(2*pi))*spectral14Re[i]/(dw*(j-i));
					G14Im[j]+=(dw/(2*pi))*spectral14Im[i]/(dw*(j-i));
					//printf("%lf \t %lf \n",ReG11[j],ReG13[j]);
					//printf("%d,%d \t %lf \t %lf \n",j,i,(dw/pi)*ImG11[i]/(dw*(j-i)),(dw/pi)*ImG13[i]/(dw*(j-i)));
				}
				G11Re[j]+=spectral11[j]*(dw/(2*pi))*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G12Re[j]+=spectral12Re[j]*(dw/(2*pi))*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G12Im[j]+=spectral12Im[j]*(dw/(2*pi))*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G13Re[j]+=spectral13Re[j]*(dw/(2*pi))*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G13Im[j]+=spectral13Im[j]*(dw/(2*pi))*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G14Re[j]+=spectral14Re[j]*(dw/(2*pi))*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G14Im[j]+=spectral14Im[j]*(dw/(2*pi))*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				//printf("%lf \t %lf \n",ReG11[j],ReG13[j]);
				//printf("%d,%d \t %lf \t %lf \n",j,j,ImG11[j]*(dw/pi)*log((j*dw+dw*0.5)/(j*dw-dw*0.5)),ImG13[j]*(dw/pi)*log((j*dw+dw*0.5)/(j*dw-dw*0.5)));
				for(i=j+1;i<stepsw;i++)
				{
					G11Re[j]+=(dw/(2*pi))*spectral11[i]/(dw*(j-i));
					G12Re[j]+=(dw/(2*pi))*spectral12Re[i]/(dw*(j-i));
					G12Im[j]+=(dw/(2*pi))*spectral12Im[i]/(dw*(j-i));
					G13Re[j]+=(dw/(2*pi))*spectral13Re[i]/(dw*(j-i));
					G13Im[j]+=(dw/(2*pi))*spectral13Im[i]/(dw*(j-i));
					G14Re[j]+=(dw/(2*pi))*spectral14Re[i]/(dw*(j-i));
					G14Im[j]+=(dw/(2*pi))*spectral14Im[i]/(dw*(j-i));
					//printf("%lf \t %lf \n",ReG11[j],ReG13[j]);
					//printf("%d,%d \t %lf \t %lf \n",j,i,(dw/pi)*ImG11[i]/(dw*(j-i)),(dw/pi)*ImG13[i]/(dw*(j-i)));
				}
				
				G11Im[j]-=spectral11[j]/2.0;
				G12Re[j]+=spectral12Im[j]/2.0;
				G12Im[j]-=spectral12Re[j]/2.0;
				G13Re[j]+=spectral13Im[j]/2.0;
				G13Im[j]-=spectral13Re[j]/2.0;
				G14Re[j]+=spectral14Im[j]/2.0;
				G14Im[j]-=spectral14Re[j]/2.0;
			}
		
			G11Re0=fopen("G11Re0.txt","w");
			G11Im0=fopen("G11Im0.txt","w");
			G12Re0=fopen("G12Re0.txt","w");
			G12Im0=fopen("G12Im0.txt","w");
			G13Re0=fopen("G13Re0.txt","w");
			G13Im0=fopen("G13Im0.txt","w");
			G14Re0=fopen("G14Re0.txt","w");
			G14Im0=fopen("G14Im0.txt","w");
			
			for(i=0;i<stepsw;++i)
			{
				w=wmin+i*dw;
				fprintf(G11Re0,"%lf \t %lf\n",w,G11Re[i]);
				fprintf(G11Im0,"%lf \t %lf\n",w,G11Im[i]);
				fprintf(G12Re0,"%lf \t %lf\n",w,G12Re[i]);
				fprintf(G12Im0,"%lf \t %lf\n",w,G12Im[i]);
				fprintf(G13Re0,"%lf \t %lf\n",w,G13Re[i]);
				fprintf(G13Im0,"%lf \t %lf\n",w,G13Im[i]);
				fprintf(G14Re0,"%lf \t %lf\n",w,G14Re[i]);
				fprintf(G14Im0,"%lf \t %lf\n",w,G14Im[i]);
			}
			fclose(G11Re0);	
			fclose(G11Im0);				
			fclose(G12Re0);
			fclose(G12Im0);
			fclose(G13Re0);		
			fclose(G13Im0);
			fclose(G14Re0);
			fclose(G14Im0);
			
			/*BOUNDSTATES*/
			
			K0=K;
			G0=G;
			J0=J+2.;
			bstJ=fopen("Tboundstates.txt","w");
			gapJ=fopen("gapJ.txt","w");
			//bstKG=fopen("bstKG.txt","w");
			//gapKG=fopen("gapKG.txt","w");
			
			//ONLY VALID FOR DJ IMPURITY
			for(i=0;i<stepsw;++i){
				j=stepsw-i;
					
				discriminant[i]=(G11Re[i]-G12Re[i]+G11Re[j]-G12Re[j])*(G11Re[i]-G12Re[i]+G11Re[j]-G12Re[j])
				-4.*(G11Re[i]-G12Re[i])*(G11Re[j]-G12Re[j])-4.*(G11Im[i]-G12Im[i])*(G11Im[j]-G12Im[j])
				+4.*(G13Re[i]-G14Re[i])*(G13Re[j]-G14Re[j])+4.*(G13Im[i]-G14Im[i])*(G13Im[j]-G14Im[j]);
				
				solution1[i]=-(G11Re[i]-G12Re[i]+G11Re[j]-G12Re[j])+sqrt(discriminant[i]);
				solution2[i]=-(G11Re[i]-G12Re[i]+G11Re[j]-G12Re[j])-sqrt(discriminant[i]);
			}
			
			discr.open("discriminant_hb60=" + std::to_string(h) + ".dat");
			sol1.open("solution1_hb60=" + std::to_string(h) + ".dat");
			sol2.open("solution2_hb60=" + std::to_string(h) + ".dat");
			for(i=0;i<stepsw;++i)
			{
				w=wmin+i*dw;
				discr<< w << "\t" << discriminant[i] << std::endl;
				sol1<< w << "\t" << solution1[i] << std::endl;
				sol2<< w << "\t" << solution2[i] << std::endl;
			}
			discr.close();
			sol1.close();
			sol2.close();			
		
			
			//for(a=-5.0;a<5.0;a+=0.5)
			for(J0=J;J0<J+4.;J0+=1.)
			{	
				//K0=a*K;
				//G0=a*G;
				printf("J0=%lf \t K0=%lf \t G0=%lf \n",J0,K0,G0);
				
				v0=-0.5*S*(J0-J+2*(K0-K)*cos(phi)*cos(phi)/3.0-(G0-G)*(1-2*cos(2*phi))/3.0);
				v1=0.5*S*(J0-J+(K0-K)*(1+2*sin(phi)*sin(phi))/6.0+(G0-G)*(1-2*cos(2*phi))/6.0);
				v2Re=S*((K0-K)*cos(2*phi)+(G0-G)*(3+2*cos(2*phi)))/12.0;
				v2Im=S*(K0-K-G0+G)*sin(phi)*sqrt(2.0)/3.0;
	
				for(i=0;i<stepsw;++i)
				{
					/*j1=(int)((wmax-wmin)/(2.0*dw));
					j=2*j1-i;*/
					j=stepsw-i;
					
					det1Re[i]=-1+(v0+v1)*(G11Re[i]+G12Re[i]+G11Re[j]+G12Re[j])+v2Re*(G13Re[i]+G14Re[i]+G13Re[j]+G14Re[j])-v2Im*(G13Im[i]+G14Im[i]+G13Im[j]+G14Im[j])
					-((v0+v1)*(v0+v1)-v2Re*v2Re-v2Im*v2Im)*((G11Re[i]+G12Re[i])*(G11Re[j]+G12Re[j])+(G11Im[i]+G12Im[i])*(G11Im[j]+G12Im[j])-(G13Re[i]+G14Re[i])*(G13Re[j]+G14Re[j])-(G13Im[i]+G14Im[i])*(G13Im[j]+G14Im[j]));
					
					det1Im[i]=(v0+v1)*(G11Im[i]+G12Im[i]-G11Im[j]-G12Im[j])+v2Im*(G13Re[i]+G14Re[i]-G13Re[j]-G14Re[j])+v2Re*(G13Im[i]+G14Im[i]-G13Im[j]-G14Im[j])
					-((v0+v1)*(v0+v1)-v2Re*v2Re-v2Im*v2Im)*((G11Im[i]+G12Im[i])*(G11Re[j]+G12Re[j])-(G11Re[i]+G12Re[i])*(G11Im[j]+G12Im[j])-(G13Im[i]+G14Im[i])*(G13Re[j]+G14Re[j])+(G13Re[i]+G14Re[i])*(G13Im[j]+G14Im[j]));
					
					det2Re[i]=-1+(v0-v1)*(G11Re[i]-G12Re[i]+G11Re[j]-G12Re[j])-v2Re*(G13Re[i]-G14Re[i]+G13Re[j]-G14Re[j])+v2Im*(G13Im[i]-G14Im[i]+G13Im[j]-G14Im[j])
					-((v0-v1)*(v0-v1)-v2Re*v2Re-v2Im*v2Im)*((G11Re[i]-G12Re[i])*(G11Re[j]-G12Re[j])+(G11Im[i]-G12Im[i])*(G11Im[j]-G12Im[j])-(G13Re[i]-G14Re[i])*(G13Re[j]-G14Re[j])-(G13Im[i]-G14Im[i])*(G13Im[j]-G14Im[j]));
					
					det2Im[i]=(v0-v1)*(G11Im[i]-G12Im[i]-G11Im[j]+G12Im[j])-v2Im*(G13Re[i]-G14Re[i]-G13Re[j]+G14Re[j])-v2Re*(G13Im[i]-G14Im[i]-G13Im[j]+G14Im[j])
					-((v0-v1)*(v0-v1)-v2Re*v2Re-v2Im*v2Im)*((G11Im[i]-G12Im[i])*(G11Re[j]-G12Re[j])-(G11Re[i]-G12Re[i])*(G11Im[j]-G12Im[j])-(G13Im[i]-G14Im[i])*(G13Re[j]-G14Re[j])+(G13Re[i]-G14Re[i])*(G13Im[j]-G14Im[j]));
					
					PolesRe[i]=det1Re[i]*det2Re[i]-det1Im[i]*det2Im[i];
					PolesIm[i]=det1Re[i]*det2Im[i]+det1Im[i]*det2Re[i];
				}
				
				polesRe.open("polesRe_hb60=" + std::to_string(h) + "_dJ=" + std::to_string(J0-J) + ".dat");
				polesIm.open("polesIm_hb60=" + std::to_string(h) + "_dJ=" + std::to_string(J0-J) + ".dat");

				for(i=0;i<stepsw;++i)
				{
					w=wmin+i*dw;
					polesRe<< w<< "\t"<< PolesRe[i]<< std::endl;
					polesIm<< w<< "\t"<< PolesIm[i]<< std::endl;
				}
				polesRe.close();
				polesIm.close();				
				
				
				j1=(int)(-wmin/dw)+1;
				j2=(int)((Gap-wmin)/dw);
				printf("\n %d \t %d \n %lf \t %lf \n",j1,j2,PolesRe[j1],PolesRe[j2]);
				if(PolesRe[j1]*PolesRe[j2]<0)
				{
					while(abs(j1-j2)>1)
					{
						i=(int)((j1+j2)/2.0);
						printf("%d \t %d \t %d \n",j1,j2,i);
						if(PolesRe[j1]*PolesRe[i]<0)
						{
							j2=i;
						}
						else
						{
							j1=i;
						}
					}
					i=(int)((j1+j2)/2.0);
					boundstate=wmin+i*dw;
					bstexists=1;
				}
				else
				{
					boundstate=0;
				}
				
				printf("boundstate=%lf \n",boundstate);
				fprintf(bstJ,"%lf \t %lf \n",J0-J,boundstate);
			//	fprintf(gapJ,"%lf \t %lf \n",J0-J,Gap);
				//fprintf(bstKG,"%lf \t %lf \n",a,boundstate);
				//fprintf(gapKG,"%lf \t %lf \n",a,Gap);
				//printf("existence=%d \n",bstexists);
				//fprintf(bst,"%lf \t %lf \t %lf \n",phi,h,boundstate);
//				fprintf(bsth,"%lf \t %lf \n",h,boundstate);
				//fprintf(bstphi,"%lf \t %lf \n",phi,boundstate);
			
			}
//			fclose(bstJ);
//			fclose(gapJ);
			//fclose(bstKG);
			//fclose(gapKG);
			
			/*if(bstexists==1)
			{
				break;
			}
			
			h-=dh;
			
			for(kx=0.5*dkx;kx<2*pi/sqrt(3.0);kx+=dkx)
			{
				for(ky=(-kx/sqrt(3.0)+0.5*dky);ky<(-kx/sqrt(3.0)+4*pi/3.0);ky+=dky)
				{
					diff=spectrumterm1(h,phi,kx,ky)-spectrumterm2(h,phi,kx,ky);
					
					if(diff<0)
					{
						realcond=0;
						break;
					}
					
				}
				if(realcond==0)
				{
					hcrit=h+dh;
				//	fprintf(hcritical,"%lf \t %lf \n",phi,hcrit);
					break;
				}
			}
		}*/
	//	fprintf(bstexistence,"%lf \t %lf \t %d \n",K,G,bstexists);
	}
	fclose(test);
	
	//fclose(gap);
	fclose(gaph);
	//fclose(gapphi);
	
	//fclose(bst);
	fclose(bsth);
	//fclose(bstphi);
	
/*	}
	}
	
	fclose(bstexistence);*/
	
	return 0;
}

double Ree1(double phi,double kx,double ky)
{
	double value;
	value=(S*(-((2*G+K)*sin((sqrt(3)*kx)/2.)*sin((3*ky)/2.)*sin(2*phi)*(2*sqrt(3) + sqrt(3)*cos(2*phi) + 3*sin(2*phi))) + 6*J3*cos(3*ky)*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) + cos((sqrt(3)*kx)/2.)*cos((3*ky)/2.)*(2*(G + 6*J + 2*K) + (2*G + K)*cos(2*phi))*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) - (-G - 6*J - 2*K - 12*J3*cos(sqrt(3)*kx) + (2*G + K)*cos(2*phi))*(2 + cos(2*phi) + sqrt(3)*sin(2*phi))))/(4.*pow(3*cos(phi) + sqrt(3)*sin(phi),2));
	return value;
}

double Ime1(double phi,double kx,double ky)
{
	double value;
	value=(S*((2*G + K)*cos((3*ky)/2.)*sin((sqrt(3)*kx)/2.)*sin(2*phi)*(2*sqrt(3) + sqrt(3)*cos(2*phi) + 3*sin(2*phi)) + cos((sqrt(3)*kx)/2.)*(2*(G + 6*J + 2*K) + (2*G + K)*cos(2*phi))*sin((3*ky)/2.)*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) + 6*J3*sin(3*ky)*(2 + cos(2*phi) + sqrt(3)*sin(2*phi))))/(4.*pow(3*cos(phi) + sqrt(3)*sin(phi),2));
	return value;
}


double Ree2(double phi,double kx,double ky)
{
	double value;
	value=(S*(14*G + K + cos((sqrt(3)*kx - 3*ky)/2.)*(14*G + K - 2*(2*G + K)*cos(4*phi)) + 2*((4*G - K)*cos((sqrt(3)*kx + 3*ky)/2.) + (7*G + 2*K + 2*(G - K)*cos((sqrt(3)*kx)/2.)*cos((3*ky)/2.))*cos(2*phi) + (2*G + K)*pow(cos((sqrt(3)*kx + 3*ky)/4.),2)*cos(4*phi) + sqrt(6)*(-G + K)*cos(3*phi)*sin((sqrt(3)*kx - 3*ky)/2.) - sqrt(6)*(G - K)*cos(phi)*(2*sin((sqrt(3)*kx - 3*ky)/2.) + 3*sin((sqrt(3)*kx + 3*ky)/2.)) + 3*sqrt(2)*(-G + K)*sin((sqrt(3)*kx + 3*ky)/2.)*sin(phi) + sqrt(3)*(3*G + 6*G*cos((sqrt(3)*kx)/2.)*cos((3*ky)/2.) + 2*(2*G + K)*sin((sqrt(3)*kx)/2.)*sin((3*ky)/2.))*sin(2*phi) - sqrt(2)*(G - K)*(sin((sqrt(3)*kx - 3*ky)/2.) + 2*sin((sqrt(3)*kx + 3*ky)/2.))*sin(3*phi) + sqrt(3)*(2*G + K)*pow(sin((sqrt(3)*kx + 3*ky)/4.),2)*sin(4*phi))))/(8.*pow(3*cos(phi) + sqrt(3)*sin(phi),2));
	return value;
}

double Ime2(double phi,double kx,double ky)
{
	double value;
	value=-(S*(2*sqrt(2)*(G - K)*sin(phi)*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) + cos((sqrt(3)*kx)/2.)*((-6*G + (2*G + K)*cos(2*phi))*sin((3*ky)/2.) + 2*sqrt(2)*(-G + K)*cos((3*ky)/2.)*sin(phi))*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) + sin((sqrt(3)*kx)/2.)*(2*sqrt(3) + sqrt(3)*cos(2*phi) + 3*sin(2*phi))*(2*sqrt(2)*(G - K)*cos(phi)*sin((3*ky)/2.) + (2*G + K)*cos((3*ky)/2.)*sin(2*phi))))/(4.*pow(3*cos(phi) + sqrt(3)*sin(phi),2));
	return value;
}

double e0(double h)
{
	double value;
	value=S*(h/S-3*J-K+G-3*J3)/2.0;
	return value;
}

double A(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=pow(e0(h),4) + pow(Ime1(phi,kx,ky),4) + pow(Ime2(phi,-kx,-ky),2)*pow(Ime2(phi,kx,ky),2) - 2*Ime2(phi,-kx,-ky)*Ime2(phi,kx,ky)*pow(Ree1(phi,kx,ky),2) + pow(Ree1(phi,kx,ky),4) + pow(Ime2(phi,kx,ky),2)*pow(Ree2(phi,-kx,-ky),2) - 2*pow(Ree1(phi,kx,ky),2)*Ree2(phi,-kx,-ky)*Ree2(phi,kx,ky) + pow(Ime2(phi,-kx,-ky),2)*pow(Ree2(phi,kx,ky),2) + pow(Ree2(phi,-kx,-ky),2)*pow(Ree2(phi,kx,ky),2) + 4*Ime1(phi,kx,ky)*Ree1(phi,kx,ky)*(-(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)) + Ime2(phi,-kx,-ky)*Ree2(phi,kx,ky)) + (pow(Ime2(phi,-kx,-ky),2) + pow(Ime2(phi,kx,ky),2) - 2*pow(Ree1(phi,kx,ky),2) + pow(Ree2(phi,-kx,-ky),2) + pow(Ree2(phi,kx,ky),2))*pow(w,2) + pow(w,4) + 2*pow(Ime1(phi,kx,ky),2)*(Ime2(phi,-kx,-ky)*Ime2(phi,kx,ky) + pow(Ree1(phi,kx,ky),2) + Ree2(phi,-kx,-ky)*Ree2(phi,kx,ky) - pow(w,2)) - pow(e0(h),2)*(2*pow(Ime1(phi,kx,ky),2) + pow(Ime2(phi,-kx,-ky),2) + pow(Ime2(phi,kx,ky),2) + 2*pow(Ree1(phi,kx,ky),2) + pow(Ree2(phi,-kx,-ky),2) + pow(Ree2(phi,kx,ky),2) + 2*pow(w,2));
	return value;
}

double B(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=2*w*(-2*pow(e0(h),2) - 2*pow(Ime1(phi,kx,ky),2) + pow(Ime2(phi,-kx,-ky),2) + pow(Ime2(phi,kx,ky),2) - 2*pow(Ree1(phi,kx,ky),2) + pow(Ree2(phi,-kx,-ky),2) + pow(Ree2(phi,kx,ky),2) + 2*pow(w,2));
	return value;	
}

double C11(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-pow(e0(h),3) - pow(e0(h),2)*w + w*(-pow(Ime1(phi,kx,ky),2) + pow(Ime2(phi,-kx,-ky),2) - pow(Ree1(phi,kx,ky),2) + pow(Ree2(phi,-kx,-ky),2) + pow(w,2)) + e0(h)*(pow(Ime1(phi,kx,ky),2) + pow(Ime2(phi,-kx,-ky),2) + pow(Ree1(phi,kx,ky),2) + pow(Ree2(phi,-kx,-ky),2) + pow(w,2));
	return value;
}

double D11(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-pow(Ime1(phi,kx,ky),2) + pow(Ime2(phi,-kx,-ky),2) - pow(Ree1(phi,kx,ky),2) + pow(Ree2(phi,-kx,-ky),2) - (e0(h) - 3*w)*(e0(h) + w);
	return value;
}

double C12Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=pow(e0(h),2)*Ree1(phi,kx,ky) - pow(Ime1(phi,kx,ky),2)*Ree1(phi,kx,ky) + Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky) - Ime2(phi,-kx,-ky)*Ree2(phi,kx,ky)) + 2*e0(h)*Ree1(phi,kx,ky)*w + Ree1(phi,kx,ky)*(Ime2(phi,-kx,-ky)*Ime2(phi,kx,ky) - pow(Ree1(phi,kx,ky),2) + Ree2(phi,-kx,-ky)*Ree2(phi,kx,ky) + pow(w,2));
	return value;
}

double C12Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-(pow(e0(h),2)*Ime1(phi,kx,ky)) + pow(Ime1(phi,kx,ky),3) - Ime2(phi,kx,ky)*Ree1(phi,kx,ky)*Ree2(phi,-kx,-ky) + Ime2(phi,-kx,-ky)*Ree1(phi,kx,ky)*Ree2(phi,kx,ky) - 2*e0(h)*Ime1(phi,kx,ky)*w + Ime1(phi,kx,ky)*(Ime2(phi,-kx,-ky)*Ime2(phi,kx,ky) + pow(Ree1(phi,kx,ky),2) + Ree2(phi,-kx,-ky)*Ree2(phi,kx,ky) - pow(w,2));
	return value;
}

double D12Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=2*Ree1(phi,kx,ky)*(e0(h) + w);
	return value;
}

double D12Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-2*Ime1(phi,kx,ky)*(e0(h) + w);
	return value;
}

double C13Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=e0(h)*Ime1(phi,kx,ky)*(Ime2(phi,-kx,-ky) - Ime2(phi,kx,ky)) - e0(h)*Ree1(phi,kx,ky)*(Ree2(phi,-kx,-ky) + Ree2(phi,kx,ky)) + Ime1(phi,kx,ky)*(Ime2(phi,-kx,-ky) + Ime2(phi,kx,ky))*w + Ree1(phi,kx,ky)*(-Ree2(phi,-kx,-ky) + Ree2(phi,kx,ky))*w;
	return value;
}

double C13Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=e0(h)*(Ime2(phi,-kx,-ky) + Ime2(phi,kx,ky))*Ree1(phi,kx,ky) + e0(h)*Ime1(phi,kx,ky)*(Ree2(phi,-kx,-ky) - Ree2(phi,kx,ky)) + Ime2(phi,-kx,-ky)*Ree1(phi,kx,ky)*w - Ime2(phi,kx,ky)*Ree1(phi,kx,ky)*w + Ime1(phi,kx,ky)*(Ree2(phi,-kx,-ky) + Ree2(phi,kx,ky))*w;
	return value;
}

double D13Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=Ime1(phi,kx,ky)*(Ime2(phi,-kx,-ky) + Ime2(phi,kx,ky)) + Ree1(phi,kx,ky)*(-Ree2(phi,-kx,-ky) + Ree2(phi,kx,ky));
	return value;
}

double D13Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=Ime2(phi,-kx,-ky)*Ree1(phi,kx,ky) - Ime2(phi,kx,ky)*Ree1(phi,kx,ky) + Ime1(phi,kx,ky)*(Ree2(phi,-kx,-ky) + Ree2(phi,kx,ky));
	return value;
}

double C14Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-2*Ime1(phi,kx,ky)*Ime2(phi,-kx,-ky)*Ree1(phi,kx,ky) - pow(Ime1(phi,kx,ky),2)*Ree2(phi,-kx,-ky) + pow(Ree1(phi,kx,ky),2)*Ree2(phi,-kx,-ky) + pow(e0(h),2)*Ree2(phi,kx,ky) - Ree2(phi,kx,ky)*(pow(Ime2(phi,-kx,-ky),2) + pow(Ree2(phi,-kx,-ky),2) + pow(w,2));
	return value;
}

double C14Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=pow(Ime1(phi,kx,ky),2)*Ime2(phi,-kx,-ky) - pow(e0(h),2)*Ime2(phi,kx,ky) + pow(Ime2(phi,-kx,-ky),2)*Ime2(phi,kx,ky) - Ime2(phi,-kx,-ky)*pow(Ree1(phi,kx,ky),2) - 2*Ime1(phi,kx,ky)*Ree1(phi,kx,ky)*Ree2(phi,-kx,-ky) + Ime2(phi,kx,ky)*(pow(Ree2(phi,-kx,-ky),2) + pow(w,2));
	return value;
}

double D14Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-2*Ree2(phi,kx,ky)*w;
	return value;
}

double D14Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=2*Ime2(phi,kx,ky)*w;
	return value;
}

double spectral11k(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C11(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D11(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}

double spectral12kRe(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C12Re(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D12Re(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}

double spectral12kIm(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C12Im(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D12Im(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}

double spectral13kRe(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C13Re(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D13Re(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}

double spectral13kIm(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C13Im(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D13Im(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}

double spectral14kRe(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C14Re(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D14Re(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}

double spectral14kIm(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C14Im(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D14Im(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}


double spectrum(double h,double phi,double kx,double ky)
{
	double value;
	value=sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
	-sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
	+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)-Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
	+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0));
	return value;
}

double spectrumterm1(double h,double phi,double kx,double ky)
{
	double value;
	value=e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0;
	return value;
}

double spectrumterm2(double h,double phi,double kx,double ky)
{
	double value;
	value=sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
	+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)-Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
	+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0);
	return value;
}


/*double complex determinant(double complex G11,double complex G12,double complex G13,double complex G14,double complex G21,double complex G22,double complex G23,double complex G24,double complex G31,double complex G32,double complex G33,double complex G34,double complex G41,double complex G42,double complex G43,double complex G44)
{
	double complex value;
	value=-((-(G43*v0) - G44*v1 - G42*conj(v2))*(((-1 + G11*v0 + G12*v1 + G14*v2)*(-1 + G22*v0 + G21*v1 + G23*v2) - (G12*v0 + G11*v1 + G13*v2)*(G21*v0 + G22*v1 + G24*v2))*(-(G34*v0) - G33*v1 - G31*conj(v2)) + (-(G31*v0) - G32*v1 - G34*v2)*(-((-1 + G22*v0 + G21*v1 + G23*v2)*(G14*v0 + G13*v1 + G11*conj(v2))) + (G12*v0 + G11*v1 + G13*v2)*(G24*v0 + G23*v1 + G21*conj(v2))) - (-(G32*v0) - G31*v1 - G33*v2)*(-((G21*v0 + G22*v1 + G24*v2)*(G14*v0 + G13*v1 + G11*conj(v2))) + (-1 + G11*v0 + G12*v1 + G14*v2)*(G24*v0 + G23*v1 + G21*conj(v2))))) + (1 - G44*v0 - G43*v1 - G41*conj(v2))*(((-1 + G11*v0 + G12*v1 + G14*v2)*(-1 + G22*v0 + G21*v1 + G23*v2) - (G12*v0 + G11*v1 + G13*v2)*(G21*v0 + G22*v1 + G24*v2))*(1 - G33*v0 - G34*v1 - G32*conj(v2)) + (-(G31*v0) - G32*v1 - G34*v2)*(-((-1 + G22*v0 + G21*v1 + G23*v2)*(G13*v0 + G14*v1 + G12*conj(v2))) + (G12*v0 + G11*v1 + G13*v2)*(G23*v0 + G24*v1 + G22*conj(v2))) - (-(G32*v0) - G31*v1 - G33*v2)*(-((G21*v0 + G22*v1 + G24*v2)*(G13*v0 + G14*v1 + G12*conj(v2))) + (-1 + G11*v0 + G12*v1 + G14*v2)*(G23*v0 + G24*v1 + G22*conj(v2)))) - (-(G41*v0) - G42*v1 - G44*v2)*(-((1 - G33*v0 - G34*v1 - G32*conj(v2))*(-((-1 + G22*v0 + G21*v1 + G23*v2)*(G14*v0 + G13*v1 + G11*conj(v2))) + (G12*v0 + G11*v1 + G13*v2)*(G24*v0 + G23*v1 + G21*conj(v2)))) + (-(G34*v0) - G33*v1 - G31*conj(v2))*(-((-1 + G22*v0 + G21*v1 + G23*v2)*(G13*v0 + G14*v1 + G12*conj(v2))) + (G12*v0 + G11*v1 + G13*v2)*(G23*v0 + G24*v1 + G22*conj(v2))) + (-(G32*v0) - G31*v1 - G33*v2)*((G13*v0 + G14*v1 + G12*conj(v2))*(G24*v0 + G23*v1 + G21*conj(v2)) - (G14*v0 + G13*v1 + G11*conj(v2))*(G23*v0 + G24*v1 + G22*conj(v2)))) + (-(G42*v0) - G41*v1 - G43*v2)*(-((1 - G33*v0 - G34*v1 - G32*conj(v2))*(-((G21*v0 + G22*v1 + G24*v2)*(G14*v0 + G13*v1 + G11*conj(v2))) + (-1 + G11*v0 + G12*v1 + G14*v2)*(G24*v0 + G23*v1 + G21*conj(v2)))) + (-(G34*v0) - G33*v1 - G31*conj(v2))*(-((G21*v0 + G22*v1 + G24*v2)*(G13*v0 + G14*v1 + G12*conj(v2))) + (-1 + G11*v0 + G12*v1 + G14*v2)*(G23*v0 + G24*v1 + G22*conj(v2))) + (-(G31*v0) - G32*v1 - G34*v2)*((G13*v0 + G14*v1 + G12*conj(v2))*(G24*v0 + G23*v1 + G21*conj(v2)) - (G14*v0 + G13*v1 + G11*conj(v2))*(G23*v0 + G24*v1 + G22*conj(v2))));
	return value;
}*/
