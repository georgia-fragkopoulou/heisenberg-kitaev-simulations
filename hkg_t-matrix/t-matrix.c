#include<stdio.h>
#include<math.h>

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

FILE *gap;
FILE *gap0;
FILE *gaph;
FILE *gapphi;
FILE *spectral110;
FILE *spectral330;
FILE *spectral13Re0;
FILE *spectral13Im0;
FILE *spectral31Re0;
FILE *spectral31Im0;
FILE *G110;
FILE *G330;
FILE *G13Re0;
FILE *G13Im0;
FILE *G31Re0;
FILE *G31Im0;
FILE *polesRe;
FILE *polesIm;
FILE *test;
FILE *bstphi;
FILE *bsth;
FILE *bst;
FILE *hcritical;
FILE *bstexistence;

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
double C33(double h,double phi,double w, double kx, double ky);
double D11(double h,double phi,double w, double kx, double ky);
double D33(double h,double phi,double w, double kx, double ky);
double C13Re(double h,double phi,double w, double kx, double ky);
double C13Im(double h,double phi,double w, double kx, double ky);
double D13Re(double h,double phi,double w, double kx, double ky);
double D13Im(double h,double phi,double w, double kx, double ky);
double C31Re(double h,double phi,double w, double kx, double ky);
double C31Im(double h,double phi,double w, double kx, double ky);
double D31Re(double h,double phi,double w, double kx, double ky);
double D31Im(double h,double phi,double w, double kx, double ky);

double spectral11k(double h,double phi,double w,double kx,double ky);
double spectral33k(double h,double phi,double w,double kx,double ky);
double spectral13kRe(double h,double phi,double w,double kx,double ky);
double spectral13kIm(double h,double phi,double w,double kx,double ky);
double spectral31kRe(double h,double phi,double w,double kx,double ky);
double spectral31kIm(double h,double phi,double w,double kx,double ky);


main()
{
	S=0.5;
	J=-0.5;
	K=-5.0;
	G=2.5;
	J3=0.5;
	
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

	
	stepskx=300.0;
	stepsky=300.0;
	stepsw=600.0;
	dkx=2*pi/(sqrt(3.0)*stepskx);
	dky=4*pi/(3*stepsky);
	wmin=-15.0;
	wmax=15.0;
	dw=(wmax-wmin)/stepsw;
	stepsphi=12.0;
	phimin=0.0;
	phimax=2*pi/3.0;
	dphi=(phimax-phimin)/stepsphi;
	hmax=15;
	dh=0.5;
//	stepsh=25;
//	dh=hmax/stepsh;
	
	
	double spectral11[int(stepsw)];
	double spectral33[int(stepsw)];
	double spectral13Re[int(stepsw)];
	double spectral13Im[int(stepsw)];
	double spectral31Re[int(stepsw)];
	double spectral31Im[int(stepsw)];
	double G11[int(stepsw)];
	double G33[int(stepsw)];
	double G13Re[int(stepsw)];
	double G13Im[int(stepsw)];
	double G31Re[int(stepsw)];
	double G31Im[int(stepsw)];
	double T11[int(stepsw)];
	double T13[int(stepsw)];
	double PolesRe[int(stepsw)];
	double PolesIm[int(stepsw)];
	
	
	phi=pi/2;
	h=6.0;
	
/*	bstexistence=fopen("bstexistence.txt","w");
	
	for(K=-8;K<-5;K+=0.5)
	{
	for(G=4;G<10;G+=0.5)
	{
	printf("K=%lf \t G=%lf \n",K,G);*/
	
	test=fopen("test.txt","w");
	
	//hcritical=fopen("hcritical.txt","w");
	
	//bst=fopen("boundstates.txt","w");
	//bsth=fopen("boundstatesh.txt","w");
	//bstphi=fopen("boundstatesphi.txt","w");
	
	//gap=fopen("gap.txt","w");
	//gaph=fopen("gaph.txt","w");
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
					fprintf(test,"%lf \t %lf \t %lf \n",kx,ky,energy);	
				}
			}
			fprintf(test,"\n%lf \t %lf \t %lf \n",GapBZx,GapBZy,Gap);
			
			printf("gap=%lf \n",Gap);
			//fprintf(gap,"%lf \t %lf \t %lf\n",phi,h,Gap);
			//fprintf(gaph,"%lf \t %lf\n",h,Gap);
			//fprintf(gapphi,"%lf \t %lf\n",phi,Gap);
			
			gap0=fopen("gap0.txt","w");
			fprintf(gap0,"%lf \t %d",Gap,0);
			fclose(gap0);
			
			for(i=0;i<stepsw;i++)
			{
				spectral11[i]=0;
				spectral33[i]=0;
				spectral13Re[i]=0;
				spectral13Im[i]=0;
				spectral31Re[i]=0;
				spectral31Im[i]=0;
				PolesRe[i]=0;
				PolesIm[i]=0;
				//printf("%lf \n",PolesRe[i]);
			}
			for(i=0;i<stepsw;++i)
			{
				G11[i]=0;
				G33[i]=0;
				G13Re[i]=0;
				G13Im[i]=0;
				G31Re[i]=0;
				G31Im[i]=0;
			}
			
			for(kx=0.5*dkx;kx<2*pi/sqrt(3.0);kx+=dkx)
			{
				for(ky=(-kx/sqrt(3.0)+0.5*dky);ky<(-kx/sqrt(3.0)+4*pi/3.0);ky+=dky)
				{
					//printf("kx=%lf \t ky=%lf\n",kx,ky);
					w1=sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
					+sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
					+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
					+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0));
					w2=sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
					-sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
					+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
					+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0));
					w3=-sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
					+sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
					+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
					+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0));
					w4=-sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
					-sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
					+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
					+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0));
					//printf("%lf \t %lf \n",e1me2(kx,ky),e1pe2(kx,ky));
					//printf("w1=%lf \t w2=%lf \t w3=%lf \t w4=%lf \n",w1,w2,w3,w4);
					
					j1=int((w1-wmin)/dw);
					j2=int((w2-wmin)/dw);
					j3=int((w3-wmin)/dw);
					j4=int((w4-wmin)/dw);
					//printf("%d \t %d \t %d \t %d \n\n",j1,j2,j3,j4);
					
					if(j1<stepsw)
					{
						spectral11[j1]+=pi/(dw*stepskx*stepsky)*spectral11k(h,phi,w1,kx,ky);
						spectral33[j1]+=pi/(dw*stepskx*stepsky)*spectral33k(h,phi,w1,kx,ky);
						spectral13Re[j1]+=pi/(dw*stepskx*stepsky)*spectral13kRe(h,phi,w1,kx,ky);
						spectral13Im[j1]+=pi/(dw*stepskx*stepsky)*spectral13kIm(h,phi,w1,kx,ky);
						spectral31Re[j1]+=pi/(dw*stepskx*stepsky)*spectral31kRe(h,phi,w1,kx,ky);
						spectral31Im[j1]+=pi/(dw*stepskx*stepsky)*spectral31kIm(h,phi,w1,kx,ky);
					}
					if(j2<stepsw)
					{
						spectral11[j2]+=pi/(dw*stepskx*stepsky)*spectral11k(h,phi,w2,kx,ky);
						spectral33[j2]+=pi/(dw*stepskx*stepsky)*spectral33k(h,phi,w2,kx,ky);
						spectral13Re[j2]+=pi/(dw*stepskx*stepsky)*spectral13kRe(h,phi,w2,kx,ky);
						spectral13Im[j2]+=pi/(dw*stepskx*stepsky)*spectral13kIm(h,phi,w2,kx,ky);
						spectral31Re[j2]+=pi/(dw*stepskx*stepsky)*spectral31kRe(h,phi,w2,kx,ky);
						spectral31Im[j2]+=pi/(dw*stepskx*stepsky)*spectral31kIm(h,phi,w2,kx,ky);
					}
					if(j3<stepsw)
					{
						spectral11[j3]+=pi/(dw*stepskx*stepsky)*spectral11k(h,phi,w3,kx,ky);
						spectral33[j3]+=pi/(dw*stepskx*stepsky)*spectral33k(h,phi,w3,kx,ky);
						spectral13Re[j3]+=pi/(dw*stepskx*stepsky)*spectral13kRe(h,phi,w3,kx,ky);
						spectral13Im[j3]+=pi/(dw*stepskx*stepsky)*spectral13kIm(h,phi,w3,kx,ky);
						spectral31Re[j3]+=pi/(dw*stepskx*stepsky)*spectral31kRe(h,phi,w3,kx,ky);
						spectral31Im[j3]+=pi/(dw*stepskx*stepsky)*spectral31kIm(h,phi,w3,kx,ky);
					}
					if(j4<stepsw)
					{
						spectral11[j4]+=pi/(dw*stepskx*stepsky)*spectral11k(h,phi,w4,kx,ky);
						spectral33[j4]+=pi/(dw*stepskx*stepsky)*spectral33k(h,phi,w4,kx,ky);
						spectral13Re[j4]+=pi/(dw*stepskx*stepsky)*spectral13kRe(h,phi,w4,kx,ky);
						spectral13Im[j4]+=pi/(dw*stepskx*stepsky)*spectral13kIm(h,phi,w4,kx,ky);
						spectral31Re[j4]+=pi/(dw*stepskx*stepsky)*spectral31kRe(h,phi,w4,kx,ky);
						spectral31Im[j4]+=pi/(dw*stepskx*stepsky)*spectral31kIm(h,phi,w4,kx,ky);
					}
				}
			}
			
			spectral110=fopen("spectral110.txt","w");
			spectral330=fopen("spectral330.txt","w");
			spectral13Re0=fopen("spectral13Re0.txt","w");
			spectral13Im0=fopen("spectral13Im0.txt","w");
			spectral31Re0=fopen("spectral31Re0.txt","w");
			spectral31Im0=fopen("spectral31Im0.txt","w");
			for(i=0;i<stepsw;++i)
			{
				w=wmin+i*dw;
				fprintf(spectral110,"%lf \t %lf\n",w,spectral11[i]);
				fprintf(spectral330,"%lf \t %lf\n",w,spectral33[i]);
				fprintf(spectral13Re0,"%lf \t %lf\n",w,spectral13Re[i]);
				fprintf(spectral13Im0,"%lf \t %lf\n",w,spectral13Im[i]);
				fprintf(spectral31Re0,"%lf \t %lf\n",w,spectral13Re[i]);
				fprintf(spectral31Im0,"%lf \t %lf\n",w,spectral31Im[i]);
			}
			fclose(spectral110);
			fclose(spectral330);
			fclose(spectral13Re0);
			fclose(spectral13Im0);
			fclose(spectral31Re0);
			fclose(spectral31Im0);
			
			
			for(j=0;j<stepsw;++j)
			{
				for(i=0;i<j;++i)
				{
					G11[j]+=(dw/pi)*spectral11[i]/(dw*(j-i));
					G33[j]+=(dw/pi)*spectral33[i]/(dw*(j-i));
					G13Re[j]+=(dw/pi)*spectral13Re[i]/(dw*(j-i));
					G13Im[j]+=(dw/pi)*spectral13Im[i]/(dw*(j-i));
					G31Re[j]+=(dw/pi)*spectral31Re[i]/(dw*(j-i));
					G31Im[j]+=(dw/pi)*spectral31Im[i]/(dw*(j-i));
					//printf("%lf \t %lf \n",ReG11[j],ReG13[j]);
					//printf("%d,%d \t %lf \t %lf \n",j,i,(dw/pi)*ImG11[i]/(dw*(j-i)),(dw/pi)*ImG13[i]/(dw*(j-i)));
				}
				G11[j]+=spectral11[j]*(dw/pi)*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G33[j]+=spectral33[j]*(dw/pi)*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G13Re[j]+=spectral13Re[j]*(dw/pi)*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G13Im[j]+=spectral13Im[j]*(dw/pi)*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G31Re[j]+=spectral31Re[j]*(dw/pi)*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				G31Im[j]+=spectral31Im[j]*(dw/pi)*log((wmin+j*dw+dw*0.5)/(wmin+j*dw-dw*0.5));
				//printf("%lf \t %lf \n",ReG11[j],ReG13[j]);
				//printf("%d,%d \t %lf \t %lf \n",j,j,ImG11[j]*(dw/pi)*log((j*dw+dw*0.5)/(j*dw-dw*0.5)),ImG13[j]*(dw/pi)*log((j*dw+dw*0.5)/(j*dw-dw*0.5)));
				for(i=j+1;i<stepsw;i++)
				{
					G11[j]+=(dw/pi)*spectral11[i]/(dw*(j-i));
					G33[j]+=(dw/pi)*spectral33[i]/(dw*(j-i));
					G13Re[j]+=(dw/pi)*spectral13Re[i]/(dw*(j-i));
					G13Im[j]+=(dw/pi)*spectral13Im[i]/(dw*(j-i));
					G31Re[j]+=(dw/pi)*spectral31Re[i]/(dw*(j-i));
					G31Im[j]+=(dw/pi)*spectral31Im[i]/(dw*(j-i));
					//printf("%lf \t %lf \n",ReG11[j],ReG13[j]);
					//printf("%d,%d \t %lf \t %lf \n",j,i,(dw/pi)*ImG11[i]/(dw*(j-i)),(dw/pi)*ImG13[i]/(dw*(j-i)));
				}
				
				G13Re[j]+=spectral13Im[j];
				G13Im[j]-=spectral13Re[j];
				G31Re[j]+=spectral31Im[j];
				G31Im[j]-=spectral31Re[j];
			}
		
			G110=fopen("G110.txt","w");
			G330=fopen("G330.txt","w");
			G13Re0=fopen("G13Re0.txt","w");
			G13Im0=fopen("G13Im0.txt","w");
			G31Re0=fopen("G31Re0.txt","w");
			G31Im0=fopen("G31Im0.txt","w");
			for(i=0;i<stepsw;++i)
			{
				w=wmin+i*dw;
				fprintf(G110,"%lf \t %lf\n",w,G11[i]);
				fprintf(G330,"%lf \t %lf\n",w,G33[i]);
				fprintf(G13Re0,"%lf \t %lf\n",w,G13Re[i]);
				fprintf(G13Im0,"%lf \t %lf\n",w,G13Im[i]);
				fprintf(G31Re0,"%lf \t %lf\n",w,G31Re[i]);
				fprintf(G31Im0,"%lf \t %lf\n",w,G31Im[i]);
			}
			fclose(G110);				
			fclose(G330);
			fclose(G13Re0);		
			fclose(G13Im0);
			fclose(G31Re0);
			fclose(G31Im0);
			
		/*	for(i=0;i<stepsw;++i)
			{
				printf("%lf",i);
			}*/
			
			polesRe=fopen("polesRe.txt","w");
			for(i=0;i<stepsw;++i)
			{
				//printf("%d",i);
				PolesRe[i]=G11[i]*G33[i]-spectral11[i]*spectral33[i]-G13Re[i]*G31Re[i]+G13Im[i]*G31Im[i];
				w=wmin+i*dw;
				fprintf(polesRe,"%lf \t %lf\n",w,PolesRe[i]);
				/*printf("%lf \t %lf\n",w,PolesRe[i]);*/
			}
			fclose(polesRe);
			
			polesIm=fopen("polesIm.txt","w");
			for(i=0;i<stepsw;++i)
			{
				PolesIm[i]=G11[i]*spectral33[i]+spectral11[i]*G33[i]-G13Re[i]*G31Im[i]-G13Im[i]*G31Re[i];
				w=wmin+i*dw;
				fprintf(polesIm,"%lf \t %lf\n",w,PolesIm[i]);
			}
			fclose(polesIm);
			
			j1=int(-wmin/dw)+1;
			j2=int((Gap-wmin)/dw);
			printf("\n %d \t %d \n %lf \t %lf \n",j1,j2,PolesRe[j1],PolesRe[j2]);
			if(PolesRe[j1]*PolesRe[j2]<0)
			{
				while(abs(j1-j2)>1)
				{
					i=int((j1+j2)/2.0);
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
				i=int((j1+j2)/2.0);
				boundstate=wmin+i*dw;
				bstexists=1;
			}
			else
			{
				boundstate=0;
			}
			
			printf("boundstate=%lf \n",boundstate);
			//printf("existence=%d \n",bstexists);
			//fprintf(bst,"%lf \t %lf \t %lf \n",phi,h,boundstate);
			//fprintf(bsth,"%lf \t %lf \n",h,boundstate);
			//fprintf(bstphi,"%lf \t %lf \n",phi,boundstate);
			
			/*if(bstexists==1)
			{
				break;
			}*/
			
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
		}
	//	fprintf(bstexistence,"%lf \t %lf \t %d \n",K,G,bstexists);
//	}
	fclose(test);
	
	//fclose(gap);
	//fclose(gaph);
	//fclose(gapphi);
	
	//fclose(bst);
	//fclose(bsth);
	//fclose(bstphi);
	
/*	}
	}
	
	fclose(bstexistence);*/
	
	return 0;
}

double Ree1(double phi,double kx,double ky)
{
	double value;
	value=(S*(-2*(2*G + K)*sin((sqrt(3)*kx)/2.)*sin(ky/2.)*sin(2*phi)*(2*sqrt(3) + sqrt(3)*cos(2*phi) + 3*sin(2*phi)) + 24*J3*cos(sqrt(3)*kx)*cos(ky)*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) + 2*cos((sqrt(3)*kx)/2.)*cos(ky/2.)*(2*(G + 6*J + 2*K) + (2*G + K)*cos(2*phi))*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) - (-2*(G + 6*J + 2*K)*cos(ky) - 12*J3*cos(2*ky) + 2*(2*G + K)*cos(ky)*cos(2*phi))*(2 + cos(2*phi) + sqrt(3)*sin(2*phi))))/(8.*pow(3*cos(phi) + sqrt(3)*sin(phi),2));
	return value;
}

double Ime1(double phi,double kx,double ky)
{
	double value;
	value=(S*((2*G + K)*cos(ky/2.)*sin((sqrt(3)*kx)/2.)*sin(2*phi)*(2*sqrt(3) + sqrt(3)*cos(2*phi) + 3*sin(2*phi)) + cos((sqrt(3)*kx)/2.)*(2*(G + 6*J + 2*K) + (2*G + K)*cos(2*phi))*sin(ky/2.)*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) + (-G - 6*J - 2*K + 12*J3*(-cos(sqrt(3)*kx) + cos(ky)) + (2*G + K)*cos(2*phi))*sin(ky)*(2 + cos(2*phi) + sqrt(3)*sin(2*phi))))/(4.*pow(3*cos(phi) + sqrt(3)*sin(phi),2));
	return value;
}


double Ree2(double phi,double kx,double ky)
{
	double value;
	value=(S*(2*(4*G - K)*cos((sqrt(3)*kx + ky)/2.) + 4*(G - K)*cos((sqrt(3)*kx)/2.)*cos(ky/2.)*cos(2*phi) + (2*G + K)*cos((sqrt(3)*kx + ky)/2.)*cos(4*phi) + 2*sqrt(6)*(G - K)*cos(3*phi)*(-sin((sqrt(3)*kx - ky)/2.) + sin(ky)) - 2*sqrt(6)*(G - K)*cos(phi)*(2*sin((sqrt(3)*kx - ky)/2.) + sin(ky) + 3*sin((sqrt(3)*kx + ky)/2.)) - 6*sqrt(2)*(G - K)*(sin(ky) + sin((sqrt(3)*kx + ky)/2.))*sin(phi) + 2*sqrt(3)*(G - K)*cos((sqrt(3)*kx + ky)/2.)*sin(2*phi) + 2*cos(ky)*(3*G + (2*G + K)*cos(2*phi))*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) + cos((sqrt(3)*kx - ky)/2.)*(14*G + K - 2*(2*G + K)*cos(4*phi) + 2*sqrt(3)*(5*G + K)*sin(2*phi)) - 2*sqrt(2)*(G - K)*(sin((sqrt(3)*kx - ky)/2.) + sin(ky) + 2*sin((sqrt(3)*kx + ky)/2.))*sin(3*phi) - sqrt(3)*(2*G + K)*cos((sqrt(3)*kx + ky)/2.)*sin(4*phi)))/(8.*pow(3*cos(phi) + sqrt(3)*sin(phi),2));
	return value;
}

double Ime2(double phi,double kx,double ky)
{
	double value;
	value=(S*(-2*sqrt(2)*(G - K)*cos((sqrt(3)*kx - ky)/2.)*(1 + 2*cos(2*phi))*(sqrt(3)*cos(phi) + sin(phi)) - 2*(2*G + K)*cos(ky/2.)*sin((sqrt(3)*kx)/2.)*sin(2*phi)*(2*sqrt(3) + sqrt(3)*cos(2*phi) + 3*sin(2*phi)) - 2*cos((sqrt(3)*kx)/2.)*(-6*G + (2*G + K)*cos(2*phi))*sin(ky/2.)*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) - 2*((3*G + (2*G + K)*cos(2*phi))*sin(ky) + 2*sqrt(2)*(G - K)*cos(ky)*sin(phi))*(2 + cos(2*phi) + sqrt(3)*sin(2*phi)) + 2*sqrt(2)*(G - K)*cos((sqrt(3)*kx + ky)/2.)*(3*sqrt(3)*cos(phi) + 3*sin(phi) + 2*sin(3*phi))))/(8.*pow(3*cos(phi) + sqrt(3)*sin(phi),2));
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
	value=(e0(h)*e0(h)-w*w)*(e0(h)*e0(h)-w*w)-2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(e0(h)*e0(h)+w*w)
	-(Ree2(phi,kx,ky)*Ree2(phi,kx,ky)+Ime2(phi,kx,ky)*Ime2(phi,kx,ky)+Ree2(phi,-kx,-ky)*Ree2(phi,-kx,-ky)+Ime2(phi,-kx,-ky)*Ime2(phi,-kx,-ky))*(e0(h)*e0(h)-w*w)
	+pow(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)-Ime1(phi,kx,ky)*Ime1(phi,kx,ky)-Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky),2)+pow(2*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky)-Ree2(phi,-kx,-ky)*Ime2(phi,kx,ky),2);
	return value;
}

double B(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=4*w*(w*w-e0(h)*e0(h)-(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))
	+(Ree2(phi,kx,ky)*Ree2(phi,kx,ky)+Ime2(phi,kx,ky)*Ime2(phi,kx,ky)+Ree2(phi,-kx,-ky)*Ree2(phi,-kx,-ky)+Ime2(phi,-kx,-ky)*Ime2(phi,-kx,-ky))/2.0);
	return value;	
}

double C11(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-(e0(h)+w)*(e0(h)+w)*(e0(h)-w)+(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(e0(h)-w)
	+(Ree2(phi,-kx,-ky)*Ree2(phi,-kx,-ky)+Ime2(phi,-kx,-ky)*Ime2(phi,-kx,-ky))*(e0(h)+w);
	return value;
}

double D11(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-2*(e0(h)*e0(h)-w*w)+(e0(h)+w)*(e0(h)+w)+(Ree2(phi,-kx,-ky)*Ree2(phi,-kx,-ky)+Ime2(phi,-kx,-ky)*Ime2(phi,-kx,-ky))-(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky));
	return value;
}

double C33(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-(e0(h)+w)*(e0(h)-w)*(e0(h)-w)+(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(e0(h)+w)+(Ree2(phi,kx,ky)*Ree2(phi,kx,ky)+Ime2(phi,kx,ky)*Ime2(phi,kx,ky))*(e0(h)-w);
	return value;
}

double D33(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=2*(e0(h)*e0(h)-w*w)-(e0(h)-w)*(e0(h)-w)+(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))-(Ree2(phi,kx,ky)*Ree2(phi,kx,ky)+Ime2(phi,kx,ky)*Ime2(phi,kx,ky));
	return value;
}

double C13Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-(Ree1(phi,kx,ky)*Ree2(phi,kx,ky)+Ime1(phi,kx,ky)*Ime2(phi,kx,ky))*(e0(h)-w)-(Ree1(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ime1(phi,kx,ky)*Ime2(phi,-kx,-ky))*(e0(h)+w);
	return value;
}

double C13Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-(Ime1(phi,kx,ky)*Ree2(phi,kx,ky)-Ree1(phi,kx,ky)*Ime2(phi,kx,ky))*(e0(h)-w)-(-Ree1(phi,kx,ky)*Ime2(phi,-kx,-ky)-Ime1(phi,kx,ky)*Ree2(phi,-kx,-ky))*(e0(h)+w);
	return value;
}

double D13Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=Ree1(phi,kx,ky)*Ree2(phi,kx,ky)+Ime1(phi,kx,ky)*Ime2(phi,kx,ky)-(Ree1(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ime1(phi,kx,ky)*Ime2(phi,-kx,-ky));
	return value;
}

double D13Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=Ime1(phi,kx,ky)*Ree2(phi,kx,ky)-Ree1(phi,kx,ky)*Ime2(phi,kx,ky)-(-Ree1(phi,kx,ky)*Ime2(phi,-kx,-ky)-Ime1(phi,kx,ky)*Ree2(phi,-kx,-ky));
	return value;
}

double C31Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-(Ree1(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ime1(phi,kx,ky)*Ime2(phi,-kx,-ky))*(e0(h)+w)-(Ree1(phi,kx,ky)*Ree2(phi,kx,ky)+Ime1(phi,kx,ky)*Ime2(phi,kx,ky))*(e0(h)-w);
	return value;
}

double C31Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-(Ime1(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ree1(phi,kx,ky)*Ime2(phi,-kx,-ky))*(e0(h)+w)-(Ree1(phi,kx,ky)*Ime2(phi,kx,ky)-Ime1(phi,kx,ky)*Ree2(phi,kx,ky))*(e0(h)-w);
	return value;
}

double D31Re(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-(Ree1(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ime1(phi,kx,ky)*Ime2(phi,-kx,-ky)-(Ree1(phi,kx,ky)*Ree2(phi,kx,ky)+Ime1(phi,kx,ky)*Ime2(phi,kx,ky)));
	return value;
}

double D31Im(double h,double phi,double w, double kx, double ky)
{
	double value;
	value=-(Ime1(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ree1(phi,kx,ky)*Ime2(phi,-kx,-ky)-(Ree1(phi,kx,ky)*Ime2(phi,kx,ky)-Ime1(phi,kx,ky)*Ree2(phi,kx,ky)));
	return value;
}

double spectral11k(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C11(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D11(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}

double spectral33k(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C33(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D33(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
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

double spectral31kRe(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C31Re(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D31Re(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}

double spectral31kIm(double h,double phi,double w,double kx,double ky)
{
	double value;
	value=(B(h,phi,w,kx,ky)*C31Im(h,phi,w,kx,ky)-A(h,phi,w,kx,ky)*D31Im(h,phi,w,kx,ky))/(B(h,phi,w,kx,ky)*B(h,phi,w,kx,ky));
	return value;
}


double spectrum(double h,double phi,double kx,double ky)
{
	double value;
	value=sqrt(e0(h)*e0(h)+(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))/2.0
	-sqrt(4*e0(h)*e0(h)*(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))-(pow(Ree1(phi,kx,ky),2)+pow(Ime1(phi,kx,ky),2))*(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)+pow(Ree2(phi,-kx,-ky),2)+pow(Ime2(phi,-kx,-ky),2))
	+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
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
	+2*(Ree1(phi,kx,ky)*Ree1(phi,kx,ky)+Ime1(phi,kx,ky)*Ime1(phi,kx,ky))*(Ree2(phi,kx,ky)*Ree2(phi,-kx,-ky)+Ime2(phi,kx,ky)*Ime2(phi,-kx,-ky))+4*Ree1(phi,kx,ky)*Ime1(phi,kx,ky)*(Ime2(phi,kx,ky)*Ree2(phi,-kx,-ky)-Ree2(phi,kx,ky)*Ime2(phi,-kx,-ky))
	+pow(pow(Ree2(phi,kx,ky),2)+pow(Ime2(phi,kx,ky),2)-pow(Ree2(phi,-kx,-ky),2)-pow(Ime2(phi,-kx,-ky),2),2)/4.0);
	return value;
}
