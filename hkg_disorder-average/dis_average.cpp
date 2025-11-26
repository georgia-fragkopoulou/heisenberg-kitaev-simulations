#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <utility>      
#include <chrono>
#include <omp.h>

using namespace Eigen ;

const double pi=acos(-1);
const std::complex <double> I(0,1.0);

int L,Lx,Ly,N;
double J,K,G,J3;


struct neighbour{
	int jRed;
	int jGreen;
	int jBlue;
	int j3Red;
	int j3Green;
	int j3Blue;
};

struct impurity{
	int site1;
	int site2;
	Eigen::Matrix3d Himp;
};

double polar(double x,double y,double z);
double azimuthal(double x,double y,double z);
Eigen::Vector3d initializeSpins();
neighbour findNeighbours(int i,int Lx,int Ly,int N);
int findImpurity(int Lx,int Ly,int N);
Eigen::Vector2d recoverLattice(int i);
void structureFactor0(int L,Eigen::Vector3d spin[],Eigen::Vector2d position[]);
std::vector<std::vector<double>> structureFactorStatic(int L,Eigen::Vector3d spin[]);
Eigen::Matrix3d spinRotation(Eigen::Vector3d S);
void printMatrices(Eigen::MatrixXcd A,Eigen::MatrixXcd B,int N);
Eigen::MatrixXcd bogoliubovStep1(Eigen::MatrixXcd A,Eigen::MatrixXcd B,int N);
Eigen::MatrixXcd bogoliubovStep2(Eigen::MatrixXcd Z,std::vector<std::pair<double,int>> deg);
std::vector<std::pair<double,int>> findDegeneracies(Eigen::VectorXd Lambda);
std::pair<Eigen::VectorXd,Eigen::MatrixXcd> sortEigenvalues(Eigen::VectorXcd lambda, Eigen::MatrixXcd V, int N);
Eigen::MatrixXcd sortEigenvalues1(Eigen::MatrixXcd LV, int N);	//sort eigenvalues
void textureHeisenberg(int N,int imp1,Eigen::Vector2d position[],Eigen::Vector3d S[],double deviation[]);
void textureRuCl3(int N,int Lx,int imp1,Eigen::Vector2d position[],Eigen::Vector3d S[],double deviation[]);
void printEigendecomposition(Eigen::VectorXd eigenvalues,Eigen::MatrixXcd eigenvectorsZ,/*Eigen::MatrixXcd eigenvectorsT,*/int N,int n);
void densityOfStates0(Eigen::VectorXd Lambda,int n);
std::vector<double> densityOfStates(std::vector<double> spectrum);
bool sortComplex(std::complex<double> a, std::complex<double> b);
void removeColumn(Eigen::MatrixXcd& matrix, unsigned int colToRemove);
void removeElement(Eigen::VectorXd& vector, unsigned int elToRemove);
std::vector<std::vector<std::vector<double>>> structureFactorDynamical(Eigen::MatrixXcd U,Eigen::MatrixXcd V,Eigen::VectorXd Lambda, Eigen::Matrix3d R[],int L);
void averageDensityOfStates(std::vector<std::vector<double>> doss,int n);
void averageMagnetization(std::vector<Eigen::Vector3d> magnetizations,int numberOfImpurityConfigurations);
void averageStructureFactorStatic(std::vector<std::vector<std::vector<double>>> SSFs,int n,int L);
void averageStructureFactorDynamical(std::vector<std::vector<std::vector<std::vector<double>>>> DSFs,int n,int L);
std::vector<impurity> generateImpurities(int N,neighbour Snn[], double impConcentration, double dJ0, double dK0, double dG0, double width);

int main() 
{
	int i,j;
	srand(1.);

//	J=-0.5;
//	K=-5.0;
//	G=2.5;
//	J3=0.5;
	
	J=1;
	K=0;
	G=0;
	J3=0;

	
	std::cout << "J=" << J << std::endl << "K=" << K << std::endl << "G=" << G << std::endl << "J3=" << J3 << std::endl;
	
	//lattice
	L=24;
    Ly=L;
    Lx=2*L;
    N=Lx*Ly; //number of sites
    
    //impurity
	double dJ0=3;
	double dK0=0;
	double dG0=0;
	
	double impConcentration=0.1;
	double impDistributionWidth=0;
	int numberOfImpurityConfigurations=20;
	int impurityConfiguration;
	
	double convergence=pow(10,-12);
	int accuracy=-log10(convergence);
	std::cout<<"accuracy="<<accuracy<<std::endl;
//	double hcrit=2*J+K-0.5*G+6*J3+sqrt(K*K-G*K+9*G*G/4.0);
//	std::cout<<"Critical field: "<<hcrit<<std::endl;
   
    //magnetic field
    double h_magnitude=8;
    double theta_field, phi_field;
//    theta_field=atan(sqrt(2));	phi_field=pi;	//field along z-axis;
//    theta_field=pi/2.;	phi_field=pi/2.;	//in-plane field
    theta_field=0;		phi_field=0;		// field along [111] (perpendicular to plane)
    double hx=sin(theta_field)*cos(phi_field)/sqrt(6.)-sin(theta_field)*sin(phi_field)/sqrt(2.)+cos(theta_field)/sqrt(3.);
    double hy=sin(theta_field)*cos(phi_field)/sqrt(6.)+sin(theta_field)*sin(phi_field)/sqrt(2.)+cos(theta_field)/sqrt(3.);
    double hz=-2*sin(theta_field)*cos(phi_field)/sqrt(6.)+cos(theta_field)/sqrt(3.);
    Vector3d h_unit(hx,hy,hz);
    //h_unit << hx,hy,hz;
    Vector3d h;
    h=h_unit*h_magnitude;
	
	neighbour Snn[N];
	for(i=0;i<N;i++){
		Snn[i]=findNeighbours(i,Lx,Ly,N);
		//std::cout << "i=" << i << "\t jx=" << Snn[i].jx << "\t jy=" << Snn[i].jy << "\t jz=" << Snn[i].jz<<std::endl;
	}
	
	//Hamiltonian
	Matrix3d HRed;
	HRed << J+K,0,0,
   			0,J,G,
   			0,G,J;
   	Matrix3d HGreen;
    HGreen << J,0,G,
   			0,J+K,0,
   			G,0,J;	
	Matrix3d HBlue;
    HBlue << J,G,0,
   			G,J,0,
   			0,0,J+K;
	Matrix3d HJ3;
	HJ3 << J3,0,0,
			0,J3,0,
			0,0,J3;
	
	int row,column;
		
	std::vector<std::vector<double>> doss;
	std::vector<std::vector<std::vector<double>>> SSFs;
	std::vector<std::vector<std::vector<std::vector<double>>>> DSFs;
	std::vector<Eigen::Vector3d> magnetizations;
	std::vector<std::vector<impurity>> impuritiess;
	
  	auto beginTotal = std::chrono::high_resolution_clock::now();
	
	std::vector<impurity> imps;
	for(impurityConfiguration=0;impurityConfiguration<numberOfImpurityConfigurations;impurityConfiguration++){	
//		IMPURITY CONFIGURATION
		std::cout<< "Configuration "<< impurityConfiguration<< std::endl;	
		imps=generateImpurities(N,Snn,impConcentration,dJ0,dK0,dG0,impDistributionWidth);
		impuritiess.push_back(imps);
	//		test.open("impurities.dat");
		std::cout<< imps.size()<< " impurities for configuration "<< impurityConfiguration <<"at:"<< std::endl;
		for(i=0;i<imps.size();i++)
			std::cout<< imps[i].site1<< " - "<< imps[i].site2<< "\t dJ0="<< imps[i].Himp(0,0)<< std::endl;
	}
	
	#pragma omp parallel for num_threads(2) private(impurityConfiguration,i,j)
	for(impurityConfiguration=0;impurityConfiguration<numberOfImpurityConfigurations;impurityConfiguration++){	
//		IMPURITY CONFIGURATION
		std::cout<< "Configuration "<< impurityConfiguration<< std::endl;
		std::vector<impurity> impurities;
		#pragma omp critical
		{
			impurities=impuritiess[impurityConfiguration];
		}
		int numberOfImpurities=impurities.size();
	//		test.open("impurities.dat");
	   
		Eigen::Vector3d S[N];
		Eigen::Vector3d S_final[N];
		Vector3d H_MF[N];
		Vector3d S_new;
		Vector3d mean_field;
	   	double energy,energy_temp,energy_final;
	   	double w;
		double error[N];
		double error_total;
		int counter;
		double impurity_strength1,impurity_strength2;
	   	
	   	auto begin1 = std::chrono::high_resolution_clock::now();		//start measuring time
	   	
		for(i=0;i<N;i++){
    		S[i]=initializeSpins();
//    		std::cout<< "S("<<i<<")=("<< S[i](0)<< ","<< S[i](1)<< ","<< S[i](2)<< std::endl;
		}

		energy=0;
		for(i=0;i<N;i++){
				energy_temp=0.5*S[i].transpose()*(HRed*S[Snn[i].jRed]+HGreen*S[Snn[i].jGreen]+HBlue*S[Snn[i].jBlue]+J3*(S[Snn[i].j3Red]+S[Snn[i].j3Green]+S[Snn[i].j3Blue]))-S[i].dot(h);
				energy+=energy_temp;
			}
		for(i=0;i<numberOfImpurities;i++){
			energy+=S[impurities[i].site1].transpose()*impurities[i].Himp*S[impurities[i].site2];
		}
		energy_final=energy/N;
		
	   	for(j=0;j<20;j++){
	   		for(i=0;i<N;i++){
	    		S[i]=initializeSpins();
			}
			
			error_total=1.;	
			counter=0;	
			while(error_total>convergence&&counter<10000){
				for(i=0;i<N;i++){
					H_MF[i]=-(HRed*S[Snn[i].jRed]+HGreen*S[Snn[i].jGreen]+HBlue*S[Snn[i].jBlue]+J3*(S[Snn[i].j3Red]+S[Snn[i].j3Green]+S[Snn[i].j3Blue]))+h;
				}
				for(i=0;i<numberOfImpurities;i++){
					H_MF[impurities[i].site1]-=impurities[i].Himp*S[impurities[i].site2];
					H_MF[impurities[i].site2]-=impurities[i].Himp*S[impurities[i].site1];
				}
				for(i=0;i<N;i++){
					w=0.8;
//					w=0.2*rand()/RAND_MAX+0.6;
					mean_field=H_MF[i].normalized();
					S_new=(1-w)*S[i]+w*mean_field;
					S_new.normalize();
					
					error[i]=(S[i]-S_new).norm();					
					S[i]=S_new;
				}
				error_total=error[0];
				for(i=1;i<N;i++){
					if(error[i]>error_total){
						error_total=error[i];
					}
				}
				counter++;			
			}
			energy=0;
			for(i=0;i<N;i++){
				energy+=0.5*S[i].transpose()*(HRed*S[Snn[i].jRed]+HGreen*S[Snn[i].jGreen]+HBlue*S[Snn[i].jBlue]+J3*(S[Snn[i].j3Red]+S[Snn[i].j3Green]+S[Snn[i].j3Blue]))-S[i].dot(h);
			}
			for(i=0;i<numberOfImpurities;i++){
				energy+=S[impurities[i].site1].transpose()*impurities[i].Himp*S[impurities[i].site2];
			}
			
			energy_temp=energy/N;
			if(energy_temp<energy_final){
				for(i=0;i<N;i++){
					S_final[i]=S[i];
				}
				energy_final=energy_temp;
			}
		//	std::cout<<"Iteration: "<<j+1<<std::endl;
		//	std::cout<<"Number of steps: "<<counter-1<<std::endl;
//			std::cout<<"Energy: "<<energy_final<<std::endl;
//			test<< "---------------------------------------------------------------------------------------------------------------"<< std::endl;
//			test.close();
			
//			std::cout.precision(12);
//			testE<< "Initialization: "<< j<< std::endl;
//			for(i=0;i<N;i++)
//				testE<< "S("<< i<< ")="<< "("<< S[i](0)<<","<< S[i](1)<< ","<< S[i](2)<< ")"<< std::endl;
//			testE<< "Energy: "<<std::fixed<< energy_final<<"\n"<< std::endl;
		}
//		testE.close();
		
		std::cout<< "Classical state ok"<< std::endl;
		auto end1 = std::chrono::high_resolution_clock::now();
		auto elapsed1 = std::chrono::duration_cast<std::chrono::seconds>(end1 - begin1);
		std::cout<<"Time needed for classical minimization: "<< elapsed1.count()<< " seconds"<< std::endl;
	
		Eigen::Vector3d magnetization(0,0,0);                                                                                               
		for(i=0;i<N;i++){
			magnetization+=S_final[i];
		}
//		magnetization/=N;
		#pragma omp critical
		{
			magnetizations.push_back(magnetization);
		}
				
		std::vector<std::vector<double>> SFstatic;		
		SFstatic=structureFactorStatic(L,S_final);
		#pragma omp critical
		{
			SSFs.push_back(SFstatic);
		}

		auto begin2 = std::chrono::high_resolution_clock::now();

		Eigen::Matrix3d R0;	//switch to magnetic field reference frame
		Eigen::Vector3d x,y,z;
		Eigen::Vector3d z0(0,0,1);
		if(hz/h_magnitude<0.999){
			z=h_unit;
			x=z.cross(z0).normalized();
			y=z.cross(x);
			R0.col(0)=x;
			R0.col(1)=y;
			R0.col(2)=z;
		}
		else
			R0=Eigen::Matrix3d::Identity(3,3);
		
//		for(i=0;i<N;i++)
//			S_final[i]<< 1/sqrt(3),1/sqrt(3),1/sqrt(3);
		
//		std::cout<< "Field:\n"<< h<<std::endl;
		Eigen::Matrix3d R1[N], R[N];	//switch to local reference frame
		for(i=0;i<N;i++){
			R1[i]=spinRotation(R0.transpose()*S_final[i]);	
			R[i]=R0*R1[i];
//			R[i]=spinRotation(S_final[i]);
//			std::cout<<"i="<<i<<"\nSpin: \n"<< S_final[i]<<"\nRotation matrix: \n"<<R[i]<< "\nRotated spin: \n"<< R[i].transpose()*S_final[i]<<"\nfield:\n"<< R[i].transpose()*h<<std::endl;
		}
		
//		Eigen::Matrix3d R[N];
//		for(i=0;i<N;i++){
//			R[i]=spinRotation(S_final[i]);	
//			std::cout<<"i="<<i<<"\nSpin: \n"<< S_final[i]<<"\nRotation matrix: \n"<<R[i]<< "\nRotated spin: \n"<< R[i].transpose()*S_final[i]<<"\nfield:\n"<< R[i].transpose()*h<<std::endl;
//		}
		
		double Ssize=1;
		Eigen::Matrix3d mRed,mGreen,mBlue;
		Eigen::Matrix3d m3Red,m3Green,m3Blue;
		Eigen::Matrix3d mimp;
		Eigen::Vector3d hrot;
		
		//LINEAR HAMILTONIAN
		std::complex<double> linearCoeff[N];
		for(i=0;i<N;i++){
			linearCoeff[i]=0;
		}
		for(int n=0;n<N/2;n++){
			i=2*n;
			
			mRed=R[i].transpose()*HRed*R[Snn[i].jRed];
			mGreen=R[i].transpose()*HGreen*R[Snn[i].jGreen];
			mBlue=R[i].transpose()*HBlue*R[Snn[i].jBlue];
			m3Red=R[i].transpose()*HJ3*R[Snn[i].j3Red];
			m3Green=R[i].transpose()*HJ3*R[Snn[i].j3Green];
			m3Blue=R[i].transpose()*HJ3*R[Snn[i].j3Blue];
					
		//	std::cout<<i<<"\n"<<R[i]<<"\n"<<R[Snn[i].jRed]<<"\n"<<mRed<<std::endl;
			
			linearCoeff[i]+=Ssize*sqrt(Ssize/2.)*(mRed(0,2)-I*mRed(1,2)+mGreen(0,2)-I*mGreen(1,2)+mBlue(0,2)-I*mBlue(1,2)+m3Red(0,2)-I*m3Red(1,2)+m3Green(0,2)-I*m3Red(1,2)+m3Blue(0,2)-I*m3Blue(1,2));
			linearCoeff[Snn[i].jRed]+=Ssize*sqrt(Ssize/2.)*(mRed(2,0)-I*mRed(2,1));
			linearCoeff[Snn[i].jGreen]+=Ssize*sqrt(Ssize/2.)*(mGreen(2,0)-I*mGreen(2,1));
			linearCoeff[Snn[i].jBlue]+=Ssize*sqrt(Ssize/2.)*(mBlue(2,0)-I*mBlue(2,1));
			linearCoeff[Snn[i].j3Red]+=Ssize*sqrt(Ssize/2.)*(m3Red(2,0)-I*m3Red(2,1));
			linearCoeff[Snn[i].j3Green]+=Ssize*sqrt(Ssize/2.)*(m3Green(2,0)-I*m3Green(2,1));
			linearCoeff[Snn[i].j3Blue]+=Ssize*sqrt(Ssize/2.)*(m3Blue(2,0)-I*m3Blue(2,1));
		}
		
		for(i=0;i<numberOfImpurities;i++){		//impurities
			mimp=R[impurities[i].site1].transpose()*impurities[i].Himp*R[impurities[i].site2];
			linearCoeff[impurities[i].site1]+=(Ssize*sqrt(Ssize)/sqrt(2))*(mimp(0,2)-I*mimp(1,2));
			linearCoeff[impurities[i].site2]+=(Ssize*sqrt(Ssize)/sqrt(2))*(mimp(2,0)-I*mimp(2,1));
		}
		
		for(i=0;i<N;i++){
			hrot=R[i]*h;
			linearCoeff[i]+=sqrt(Ssize/2.)*(-hrot(0)+I*hrot(1));
		}
//		for(i=0;i<N;i++)
//			linearCoeff[i]+=sqrt(Ssize/2.)*(-h(0)+I*h(1));

//		for(i=0;i<N;i++)
//			linear<< dJ<< "\t"<<linearCoeff[i]<< std::endl;
		
//		BILINEAR HAMILTONIAN
		std::complex<double> test1;
		std::complex<double> test2;
		
		Eigen::MatrixXcd A(N,N);
		Eigen::MatrixXcd B(N,N);
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				A(i,j)=0;
				B(i,j)=0;				
			}
		}
		test1=A(0,1);
		test2=A(1,0);
		
		int n,l;
		double temp;
		for(n=0;n<N/2;n++){			//clean part
			i=2*n;	
			
			mRed=R[i].transpose()*HRed*R[Snn[i].jRed];
			mGreen=R[i].transpose()*HGreen*R[Snn[i].jGreen];
			mBlue=R[i].transpose()*HBlue*R[Snn[i].jBlue];
			m3Red=R[i].transpose()*HJ3*R[Snn[i].j3Red];
			m3Green=R[i].transpose()*HJ3*R[Snn[i].j3Green];
			m3Blue=R[i].transpose()*HJ3*R[Snn[i].j3Blue];
			
//			std::cout<< "i="<<i<<"mRed= \n"<< mRed<< "\nmGreen=\n"<< mGreen<< "\nmBlue=\n"<<mBlue<< std::endl;
			
			j=Snn[i].jRed;
			A(i,i)-=Ssize*mRed(2,2);
			A(j,j)-=Ssize*mRed(2,2);
			A(i,j)+=0.5*Ssize*(mRed(0,0)+mRed(1,1)-I*mRed(0,1)+I*mRed(1,0));
			A(j,i)+=0.5*Ssize*(mRed(0,0)+mRed(1,1)+I*mRed(0,1)-I*mRed(1,0));
			B(i,j)+=0.5*Ssize*(mRed(0,0)-mRed(1,1)-I*mRed(0,1)-I*mRed(1,0));
			B(j,i)+=0.5*Ssize*(mRed(0,0)-mRed(1,1)-I*mRed(0,1)-I*mRed(1,0));
			test1=A(0,1);
			test2=A(1,0);
			
			j=Snn[i].jGreen;
			A(i,i)-=Ssize*mGreen(2,2);
			A(j,j)-=Ssize*mGreen(2,2);
			A(i,j)+=0.5*Ssize*(mGreen(0,0)+mGreen(1,1)-I*mGreen(0,1)+I*mGreen(1,0));
			A(j,i)+=0.5*Ssize*(mGreen(0,0)+mGreen(1,1)+I*mGreen(0,1)-I*mGreen(1,0));
			B(i,j)+=0.5*Ssize*(mGreen(0,0)-mGreen(1,1)-I*mGreen(0,1)-I*mGreen(1,0));
			B(j,i)+=0.5*Ssize*(mGreen(0,0)-mGreen(1,1)-I*mGreen(0,1)-I*mGreen(1,0));
			test1=A(0,1);
			test2=A(1,0);
			
			j=Snn[i].jBlue;
			A(i,i)-=Ssize*mBlue(2,2);
			A(j,j)-=Ssize*mBlue(2,2);
			A(i,j)+=0.5*Ssize*(mBlue(0,0)+mBlue(1,1)-I*mBlue(0,1)+I*mBlue(1,0));
			A(j,i)+=0.5*Ssize*(mBlue(0,0)+mBlue(1,1)+I*mBlue(0,1)-I*mBlue(1,0));
			B(i,j)+=0.5*Ssize*(mBlue(0,0)-mBlue(1,1)-I*mBlue(0,1)-I*mBlue(1,0));
			B(j,i)+=0.5*Ssize*(mBlue(0,0)-mBlue(1,1)-I*mBlue(0,1)-I*mBlue(1,0));
			test1=A(0,1);
			test2=A(1,0);
			
			j=Snn[i].j3Red;
			A(i,i)-=Ssize*m3Red(2,2);
			A(j,j)-=Ssize*m3Red(2,2);
			A(i,j)+=0.5*Ssize*(m3Red(0,0)+m3Red(1,1)-I*m3Red(0,1)+I*m3Red(1,0));
			A(j,i)+=0.5*Ssize*(m3Red(0,0)+m3Red(1,1)+I*m3Red(0,1)-I*m3Red(1,0));
			B(i,j)+=0.5*Ssize*(m3Red(0,0)-m3Red(1,1)-I*m3Red(0,1)-I*m3Red(1,0));
			B(j,i)+=0.5*Ssize*(m3Red(0,0)-m3Red(1,1)-I*m3Red(0,1)-I*m3Red(1,0));
			test1=A(0,1);
			test2=A(1,0);
			
			j=Snn[i].j3Green;
			A(i,i)-=Ssize*m3Green(2,2);
			A(j,j)-=Ssize*m3Green(2,2);
			A(i,j)+=0.5*Ssize*(m3Green(0,0)+m3Green(1,1)-I*m3Green(0,1)+I*m3Green(1,0));
			A(j,i)+=0.5*Ssize*(m3Green(0,0)+m3Green(1,1)+I*m3Green(0,1)-I*m3Green(1,0));
			B(i,j)+=0.5*Ssize*(m3Green(0,0)-m3Green(1,1)-I*m3Green(0,1)-I*m3Green(1,0));
			B(j,i)+=0.5*Ssize*(m3Green(0,0)-m3Green(1,1)-I*m3Green(0,1)-I*m3Green(1,0));
			test1=A(0,1);
			test2=A(1,0);
			
			j=Snn[i].j3Blue;
			A(i,i)-=Ssize*m3Blue(2,2);
			A(j,j)-=Ssize*m3Blue(2,2);
			A(i,j)+=0.5*Ssize*(m3Blue(0,0)+m3Blue(1,1)-I*m3Blue(0,1)+I*m3Blue(1,0));
			A(j,i)+=0.5*Ssize*(m3Blue(0,0)+m3Blue(1,1)+I*m3Blue(0,1)-I*m3Blue(1,0));
			B(i,j)+=0.5*Ssize*(m3Blue(0,0)-m3Blue(1,1)-I*m3Blue(0,1)-I*m3Blue(1,0));
			B(j,i)+=0.5*Ssize*(m3Blue(0,0)-m3Blue(1,1)-I*m3Blue(0,1)-I*m3Blue(1,0));
			test1=A(0,1);
			test2=A(1,0);
			
			test1=A(0,1);
			test2=A(1,0);
		}
		
		for(i=0;i<numberOfImpurities;i++){		//impurities
			mimp=R[impurities[i].site1].transpose()*impurities[i].Himp*R[impurities[i].site2];
			A(impurities[i].site1,impurities[i].site1)-=Ssize*mimp(2,2);
			A(impurities[i].site2,impurities[i].site2)-=Ssize*mimp(2,2);
			A(impurities[i].site1,impurities[i].site2)+=0.5*Ssize*(mimp(0,0)+mimp(1,1)-I*mimp(0,1)+I*mimp(1,0));
			A(impurities[i].site2,impurities[i].site1)+=0.5*Ssize*(mimp(0,0)+mimp(1,1)+I*mimp(0,1)-I*mimp(1,0));
			B(impurities[i].site1,impurities[i].site2)+=0.5*Ssize*(mimp(0,0)-mimp(1,1)-I*mimp(0,1)-I*mimp(1,0));
			B(impurities[i].site2,impurities[i].site1)+=0.5*Ssize*(mimp(0,0)-mimp(1,1)-I*mimp(0,1)-I*mimp(1,0));
		}
		
		for(i=0;i<N;i++){			//magnetic field
			Eigen::Vector3d hrot;
			hrot=R[i].transpose()*h;	
			A(i,i)+=Ssize*hrot(2);
//			A(i,i)+=h.dot(S_final[i]);
		}
		
//		for(i=0;i<N;i++)
//			A(i,i)+=h(2);
		
		printMatrices(A,B,N);
		
		Eigen::MatrixXcd eigendecomposition;
		eigendecomposition=bogoliubovStep1(A,B,N);

		Eigen::VectorXd Lambda;
		Eigen::MatrixXcd Z;
		Lambda=eigendecomposition.topRows(1).real().transpose();
		Z=eigendecomposition.bottomRows(2*N);
		std::vector<std::pair<double,int>> deg=findDegeneracies(Lambda);
		Eigen::MatrixXcd T=bogoliubovStep2(Z,deg);
		std::cout<< "Bogoliubov done!"<< std:: endl;
		
		auto end2 = std::chrono::high_resolution_clock::now();
    	auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2);
		std::cout<<"Time needed for set-up and diagonalization of LSW Hamilotonian: "<< elapsed2.count()*1e-9<< " seconds"<< std::endl;
		
		n=Lambda.size()/2;
		std::vector<double> spectrum;
		for(int i=0;i<n;i++)
			spectrum.push_back(Lambda(i)-Lambda(i+n));
		if(n<N){
			for(int i=n;i<N;i++)
				spectrum.push_back(0);
		}
			
		
		std::vector<double> dos;
		dos=densityOfStates(spectrum);
		#pragma omp critical 
		{
			doss.push_back(dos);
		}

		Eigen::MatrixXcd U=T.topLeftCorner(N,n);
		Eigen::MatrixXcd V=T.topRightCorner(N,n);
		std::vector<std::vector<std::vector<double>>> SFdynamical;
		SFdynamical=structureFactorDynamical(U,V,Lambda,R,L);
		#pragma omp critical
		{
			DSFs.push_back(SFdynamical);
		}
	}
	
	
	averageMagnetization(magnetizations,numberOfImpurityConfigurations);
	averageDensityOfStates(doss,numberOfImpurityConfigurations);
	averageStructureFactorStatic(SSFs,numberOfImpurityConfigurations,L);
	averageStructureFactorDynamical(DSFs,numberOfImpurityConfigurations,L);
	
	auto endTotal = std::chrono::high_resolution_clock::now();
   	auto elapsedTotal = std::chrono::duration_cast<std::chrono::nanoseconds>(endTotal - beginTotal);
	std::cout<<"Total time: "<< elapsedTotal.count()*1e-9*60<< " minutes"<< std::endl;

    return 0;
}

double polar(double x,double y,double z){
	double theta;
	theta=acos(z/sqrt(x*x+y*y+z*z));
	return theta;
}

double azimuthal(double x,double y,double z){
	double phi;
	if(x>0){
		phi=atan(y/x);
	}
	else if(x<0){
		if(y>=0)
			phi=atan(y/x)+pi;
		else
			phi=atan(y/x)-pi;
	}
	else{
		if(y>0)
			phi=pi/2.;
		else
			phi=-pi/2.;
	}
	return phi;
}

double maximum(double a,double b,double c){
	double maximum;
	if(a>b){
		if(a>c){
			maximum=a;
		}
		else maximum=c;
	}
	else{
		if(b>c){
			maximum=b;
		}
		else maximum=c;
	}
	return maximum;
}

Eigen::Vector3d initializeSpins(){
	//srand((unsigned)time(0));
	Eigen::Vector3d spin;

	spin(0)=2.*rand()/RAND_MAX-1.;
	spin(1)=2.*rand()/RAND_MAX-1.;
	spin(2)=2.*rand()/RAND_MAX-1.;
	spin.normalize();

	return spin;
}

neighbour findNeighbours(int i,int Lx,int Ly,int N){
	neighbour nn;
	int row=i/Lx;
	int column=i%Lx;
	
	if(i%2==0){ //even
	nn.jRed=i+1;
	nn.jGreen=i-1;
	nn.jBlue=i-Lx+1;
	nn.j3Red=i-Lx-1;
	nn.j3Green=i-Lx+3;
	nn.j3Blue=i+Lx-1;
		if(column==0){ //Left boundary
			nn.jGreen+=Lx;
			nn.j3Red+=Lx;
			nn.j3Blue+=Lx;
		}
		if(row==0){ //lower boundary
			nn.jBlue+=N;
			nn.j3Red+=N;
			nn.j3Green+=N;
		}
		if(column==Lx-2){ //right boundary
			nn.j3Green-=Lx;
		}
		if(row==Ly-1){ //upper boundary
			nn.j3Blue-=N;
		}
	}
	else{ //odd
		nn.jRed=i-1;
		nn.jGreen=i+1;
		nn.jBlue=i+Lx-1;
		nn.j3Red=i+Lx+1;
		nn.j3Green=i+Lx-3;
		nn.j3Blue=i-Lx+1;
		if(column==1){ //left boundary
			nn.j3Green+=Lx;
		}
		if(column==Lx-1){ //right boundary
			nn.jGreen-=Lx;
			nn.j3Red-=Lx;
			nn.j3Blue-=Lx;
		}
		if(row==0){ //lower boundary
			nn.j3Blue+=N;
		}
		if(row==Ly-1){ //upper boundary
			nn.jBlue-=N;
			nn.j3Red-=N;
			nn.j3Green-=N;
		}
	}
	return nn;
}

std::vector<impurity> generateImpurities(int N,neighbour Snn[], double impConcentration, double dJ0, double dK0, double dG0, double width){
	std::vector<impurity> impurities;
	impurity imp;
	double disorderProbX,disorderProbY,disorderProbZ;
	double dJ,dK,dG;
	Eigen::Matrix3d Himp;
	
	for(int n=0;n<N/2;n++){
		int i=2*n;
		disorderProbX=1.0*rand()/RAND_MAX;
		disorderProbY=1.0*rand()/RAND_MAX;
		disorderProbZ=1.0*rand()/RAND_MAX;
		if(disorderProbX<=impConcentration){
			imp.site1=i;
			imp.site2=Snn[i].jRed;
			dJ=dJ0-width/2.+width*rand()/RAND_MAX;
			dK=dK0-width/2.+width*rand()/RAND_MAX;
			dG=dG0-width/2.+width*rand()/RAND_MAX;
			imp.Himp<< dJ+dK,0,0,
   				0,dJ,dG,
   				0,dG,dJ;
			impurities.push_back(imp);
		}
		if(disorderProbY<=impConcentration){
			imp.site1=i;
			imp.site2=Snn[i].jGreen;
			dJ=dJ0-width/2.+width*rand()/RAND_MAX;
			dK=dK0-width/2.+width*rand()/RAND_MAX;
			dG=dG0-width/2.+width*rand()/RAND_MAX;
			imp.Himp<< dJ,0,dG,
   				0,dJ+dK,0,
   				dG,0,dJ;
			impurities.push_back(imp);
		}
		if(disorderProbZ<=impConcentration){
			imp.site1=i;
			imp.site2=Snn[i].jBlue;
			dJ=dJ0-width/2.+width*rand()/RAND_MAX;
			dK=dK0-width/2.+width*rand()/RAND_MAX;
			dG=dG0-width/2.+width*rand()/RAND_MAX;
			imp.Himp<< dJ,dG,0,
   				dG,dJ,0,
   				0,0,dJ+dK;
			impurities.push_back(imp);
		}
	}
	return impurities;
}

int findImpurity(int Lx,int Ly,int N){
	int row,column;
	int imp,i;
	for(i=0;i<N;i++){
		row=i/Lx;
		column=i%Lx;
		if(row==Ly/2 && column==Lx/2){
			imp=i;
			break;
		}
	}
	return imp;
}

Eigen::Vector2d recoverLattice(int i){
	Eigen::Vector2d position;
	Eigen::Vector2d a1(sqrt(3.),0); 
	Eigen::Vector2d a2(sqrt(3.)/2.,3./2.);
	Eigen::Vector2d b(sqrt(3.)/2.,1./2.);
	int row=i/Lx;
	int column=i%(Lx/2);
	if(i%2==0){//even
		int row=i/Lx;
		int column=(i%Lx)/2;
		position=column*a1+row*a2;
	}
	else{//odd
		int row=i/Lx;
		int column=(i%Lx-1)/2;
		position=column*a1+row*a2+b;
	}
	return position;	
}

void textureHeisenberg(int N,int imp1,Eigen::Vector2d position[],Eigen::Vector3d S[],double deviation[]){
	std::ofstream resultsA;
	std::ofstream resultsB;
	double a,b,c;
	int i;
	
	resultsA.open("textureSubA.dat");
	resultsB.open("textureSubB.dat");
		
	for(i=0;i<N;i++){
		if(deviation[i]>0.05){	
			a=(S[i](0)+S[i](1)-2*S[i](2))/sqrt(6.);
			b=(-S[i](0)+S[i](1))/sqrt(2.);
			c=(S[i](0)+S[i](1)+S[i](2))/sqrt(3.);
			if(i%2==0){	
				std::cout<<"index: "<< i<< "position: "<< position[i]<< std::endl;
				resultsA<<position[i](0)<<"\t"<< position[i](1)<< "\t 0"<<std::endl;
				resultsA<< a << "\t" << b << "\t" << c << std::endl;
				resultsA<< deviation[i]<< std::endl;
			}
			else{
				std::cout<<"index: "<< i<< "position: "<< position[i]<< std::endl;
				resultsB<<position[i](0)<<"\t"<< position[i](1)<< "\t 0"<<std::endl;
				resultsB<< a << "\t" << b << "\t" << c << std::endl;
				resultsB<< deviation[i]<< std::endl;
			}
		}
	}
	if(imp1%2==0){
		resultsA<<position[imp1](0)<<"\t"<< position[imp1](1)<< "\t 0"<<std::endl;
		resultsA<< 0 << "\t" << 0 << "\t" << 0 << std::endl;
		resultsA<< 0<< std::endl;
	}
	else{
		resultsB<<position[imp1](0)<<"\t"<< position[imp1](1)<< "\t 0"<<std::endl;
		resultsB<< 0 << "\t" << 0 << "\t" << 0 << std::endl;
		resultsB<< 0<< std::endl;
	}
			
	resultsA.close();
	resultsB.close();
}

void textureRuCl3(int N,int Lx,int imp1,Eigen::Vector2d position[],Eigen::Vector3d S[],double deviation[]){
	std::ofstream results1A;
	std::ofstream results1B;
	std::ofstream results2A;
	std::ofstream results2B;
	double a,b,c;
	int i,chain;
	
	results1A.open("textureChain1A.dat");
	results1B.open("textureChain1B.dat");
	results2A.open("textureChain2A.dat");
	results2B.open("textureChain2B.dat");
			
	for(i=0;i<N;i++){
		if(deviation[i]>0.05){	
			a=(S[i](0)+S[i](1)-2*S[i](2))/sqrt(6.);
			b=(-S[i](0)+S[i](1))/sqrt(2.);
			c=(S[i](0)+S[i](1)+S[i](2))/sqrt(3.);
			chain=i/Lx;
			if(chain%2==0){
				if(i%2==0){
					results1A<<position[i](0)<<"\t"<< position[i](1)<< "\t 0"<<std::endl;
					results1A<< a << "\t" << b << "\t" << c << std::endl;
					results1A<< deviation[i]<< std::endl;
				}
				else{
					results1B<<position[i](0)<<"\t"<< position[i](1)<< "\t 0"<<std::endl;
					results1B<< a << "\t" << b << "\t" << c << std::endl;
					results1B<< deviation[i]<< std::endl;
				}
			}
			else{
				if(i%2==0){
					results2A<<position[i](0)<<"\t"<< position[i](1)<< "\t 0"<<std::endl;
					results2A<< a << "\t" << b << "\t" << c << std::endl;
					results2A<< deviation[i]<< std::endl;
				}
				else{
					results2B<<position[i](0)<<"\t"<< position[i](1)<< "\t 0"<<std::endl;
					results2B<< a << "\t" << b << "\t" << c << std::endl;
					results2B<< deviation[i]<< std::endl;
				}
			}
		}
	}
	if(imp1%2==0){
		results1A<<position[imp1](0)<<"\t"<< position[imp1](1)<< "\t 0"<<std::endl;
		results1A<< 0 << "\t" << 0 << "\t" << 0 << std::endl;
		results1A<< 0<< std::endl;
	}
	else{
		results1B<<position[imp1](0)<<"\t"<< position[imp1](1)<< "\t 0"<<std::endl;
		results1B<< 0 << "\t" << 0 << "\t" << 0 << std::endl;
		results1B<< 0<< std::endl;
	}

	results1A.close();
	results1B.close();
	results2A.close();
	results2B.close();
}

Eigen::Vector2d distance(int i,int j,int Lx,int Ly){	// distance between two sites i and j on a lattice with size Lx x Ly
	Eigen::Vector2d a1(sqrt(3.),0);
	Eigen::Vector2d a2(sqrt(3.)/2.,3./2.);	
	
	Eigen::Vector2d xi=recoverLattice(i);
	Eigen::Vector2d xj=recoverLattice(j);
	Eigen::Vector2d ri[4],rj[4];

	ri[0]=xi;
	ri[1]=xi+Lx*a1;
	ri[2]=xi+Ly*a2;
	ri[3]=xi+Lx*a1+Ly*a2;
	rj[0]=xj;
	rj[1]=xj+Lx*a1;
	rj[2]=xj+Ly*a2;
	rj[3]=xj+Lx*a1+Ly*a2;
	
	Eigen::Vector2d distance;
	Eigen::Vector2d d[4][4];
	double dnorm;
	dnorm=(Lx*a1+Ly*a2).norm();
	for(int n=0;n<4;n++){
		for(int m=0;m<4;m++){
			d[n][m]=ri[n]-rj[m];
			if(d[n][m].norm()<dnorm){
				distance=d[n][m];
				dnorm=d[n][m].norm();
			}
		}
	}
	
	return distance;
}

void structureFactor0(int L,Eigen::Vector3d spin[],Eigen::Vector2d position[]){
	int i,j,n,n1,n2,m;
	int N=2*L*L;
	std::complex <double> I(0,1.0);
	Eigen::Vector2d a1(sqrt(3.),0);
	Eigen::Vector2d a2(sqrt(3.)/2.,3./2.);
	Eigen::Vector2d b1(2*pi/sqrt(3.),-2*pi/3.);
	Eigen::Vector2d b2(0,4*pi/3.);
	
	Eigen::Vector2d r[N][4];
	for(i=0;i<N;i++){
		r[i][0]=position[i];
		r[i][1]=position[i]+L*a1;
		r[i][2]=position[i]+L*a2;
		r[i][3]=position[i]+L*(a1+a2);
	}
	
	//calculate distance with PBC
	Eigen::Vector2d distance[N][N];
	Eigen::Vector2d d[4][4];
	double dnorm;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			dnorm=(L*a1+L*a2).norm();
			for(n=0;n<4;n++){
				for(m=0;m<4;m++){
					d[n][m]=r[i][n]-r[j][m];
					if(d[n][m].norm()<dnorm){
						distance[i][j]=d[n][m];
						dnorm=d[n][m].norm();
					}
				}
			}
		}
	}
	
	double temp;
	Eigen::Vector2d q[(2*L+1)*(2*L+1)];
	double S[(2*L+1)*(2*L+1)];
	for(m=0;m<(2*L+1)*(2*L+1);m++){
		S[m]=0;
	}
	
	m=0;
	for(n1=-L;n1<=L;n1++){
		for(n2=-L;n2<=L;n2++){
			if(n1-n2<=L&&n1-n2>=-L){
				q[m]=(n1*b1+n2*b2)/L;
				for(i=0;i<N;i++){
					for(j=0;j<N;j++){
						temp=(spin[i].dot(spin[j])*exp(I*q[m].dot(distance[i][j]))).real();
						S[m]+=temp/N;  	
					}
				}
				m++;	
			}
		}
	}
			
	std::fstream resultss;
	resultss.open("structure_factor.dat",std::ios::out);
	for(i=0;i<m;i++){
		resultss<< q[i](0)<< "\t"<< q[i](1)<< "\t"<< S[i]<< std::endl;
	}
	resultss.close();
}

std::vector<std::vector<double>> structureFactorStatic(int L,Eigen::Vector3d spin[]){
	int N=2*L*L;
	std::complex <double> I(0,1.0);
	Eigen::Vector2d b1(2*pi/sqrt(3.),-2*pi/3.);
	Eigen::Vector2d b2(0,4*pi/3.);
	
	
	Eigen::Vector2d r;
	std::complex<double> fourrier0,fourrier1,fourrierCoeff,fourrierdq1,fourrierdq2;
	std::vector<std::vector<Eigen::Vector3cd>> Sq;
	std::vector<Eigen::Vector3cd> Sq2;
	Sq2.assign(2*L+1,Eigen::Vector3cd::Zero());
	Sq.assign(2*L+1,Sq2);
	std::vector<std::vector<double>>S;
	std::vector<double> S2;
	S2.assign(2*L+1,0);
	S.assign(2*L+1,S2);
	
	auto begin = std::chrono::high_resolution_clock::now();
	for(int i=0;i<N;i++){
		r=recoverLattice(i);
		fourrier0=exp(I*r.dot(-b1-b2));
		fourrierdq1=exp(I*r.dot(b1)/(1.0*L));
		fourrierdq2=exp(I*r.dot(b2)/(1.0*L));
		fourrier1=fourrier0;
		for(int n1=-L;n1<L+1;n1++){
			fourrierCoeff=fourrier1;
			for(int n2=-L;n2<L+1;n2++){
				Sq[n1+L][n2+L]+=fourrierCoeff*spin[i];
//				S[n1+L][n2+L]+=(std::norm(Sq[n1+L][n2+L](0))+std::norm(Sq[n1+L][n2+L](1))+std::norm(Sq[n1+L][n2+L](2)))/(1.0*N);
				fourrierCoeff*=fourrierdq2;	
			}
			fourrier1*=fourrierdq1;
		}
	}
	for(int n1=-L;n1<L+1;n1++){
		for(int n2=-L;n2<L+1;n2++){
			S[n1+L][n2+L]=(std::norm(Sq[n1+L][n2+L](0))+std::norm(Sq[n1+L][n2+L](1))+std::norm(Sq[n1+L][n2+L](2)))/(1.0*N);
		}
		fourrier1*=fourrierdq1;
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
	std::cout<<"Time needed for the static structure factor: "<< elapsed.count()*1e-9<< " seconds"<< std::endl;
		
	std::ofstream results;
	results.open("SFstatic.dat");
	Eigen::Vector2d q;
	for(int n1=-L;n1<=L;n1++){
		for(int n2=-L;n2<=L;n2++){
			if(n1-n2<=L&&n1-n2>=-L){
				q=(n1*b1+n2*b2)/L;
				results<< q(0)<< "\t"<< q(1)<< "\t"<< S[n1+L][n2+L]<< std::endl;
			}
		}
	}
	results.close();
	
	return S;
}


std::vector<std::vector<std::vector<double>>> structureFactorDynamical(Eigen::MatrixXcd U,Eigen::MatrixXcd V,Eigen::VectorXd Lambda, Eigen::Matrix3d R[],int L){
    std::cout<< "SF calculation started" << std::endl;
    Eigen::Vector2d b1(2*pi/sqrt(3.),-2*pi/3.);
	Eigen::Vector2d b2(0,4*pi/3.);
	
	std::vector<double> omegak;
	int n=Lambda.size()/2;
	for(int i=0;i<n;i++)
		omegak.push_back(Lambda(i)-Lambda(i+n));
	double omegakMax=*max_element(omegak.begin(),omegak.end());
	double epsilon=0.1;
	int N=2*L*L;
	
	Eigen::Vector2d r;
	std::complex<double> fourrier0,fourrier1,fourrierCoeff,fourrierdq1,fourrierdq2;
	Eigen::Matrix3cd ui,vi;
	std::vector<std::vector<Eigen::Matrix3cd>> uq,vq;
	std::vector<Eigen::Matrix3cd> uq2,vq2;
	uq2.assign(2*L+1,Eigen::Matrix3cd::Zero());
	vq2.assign(2*L+1,Eigen::Matrix3cd::Zero());
	uq.assign(2*L+1,uq2);
	vq.assign(2*L+1,vq2);
	std::vector<std::vector<double>>Sk;
	std::vector<double> Sk2;
	Sk2.assign(2*L+1,0);
	Sk.assign(2*L+1,Sk2);
	std::vector<std::vector<std::vector<double>>>S;
	S.assign(110*(int)omegakMax,Sk);
	Eigen::Vector3cd T;
	double Lorentzian;
	
	auto begin3 = std::chrono::high_resolution_clock::now();
	
	for(int k=0;k<n;k++){
  //      std::cout<<"k="<<k<<std::endl;
        uq2.assign(2*L+1,Eigen::Matrix3cd::Zero());
        vq2.assign(2*L+1,Eigen::Matrix3cd::Zero());
        uq.assign(2*L+1,uq2);
        vq.assign(2*L+1,vq2);
        Sk2.assign(2*L+1,0);
        Sk.assign(2*L+1,Sk2);
		for(int i=0;i<N;i++){
        //    std::cout<< "i="<< i<< std::endl;
			r=recoverLattice(i);
			ui=R[i]*(U(i,k)+conj(V(i,k)));
			vi=R[i]*(U(i,k)-conj(V(i,k)));
			fourrier0=exp(I*r.dot(-b1-b2));
			fourrierdq1=exp(I*r.dot(b1)/(1.0*L));
			fourrierdq2=exp(I*r.dot(b2)/(1.0*L));
			fourrier1=fourrier0;
			for(int n1=-L;n1<L+1;n1++){
    //            std::cout<<"n1="<<n1<<std::endl;
				fourrierCoeff=fourrier1;
				for(int n2=-L;n2<L+1;n2++){
					uq[n1+L][n2+L]+=fourrierCoeff*ui;
					vq[n1+L][n2+L]+=fourrierCoeff*vi;
//					T=uq[n1+L][n2+L].col(0)-I*vq[n1+L][n2+L].col(1);
//					Sk[n1+L][n2+L]+=(std::norm(T(0))+std::norm(T(1))+std::norm(T(2)))/(1.0*N);
					fourrierCoeff*=fourrierdq2;	
				}
				fourrier1*=fourrierdq1;
			}
		}
		for(int n1=-L;n1<L+1;n1++){
			for(int n2=-L;n2<L+1;n2++){
				T=uq[n1+L][n2+L].col(0)-I*vq[n1+L][n2+L].col(1);
				Sk[n1+L][n2+L]=(std::norm(T(0))+std::norm(T(1))+std::norm(T(2)))/(1.0*N);
			}
		}
		for(int omega=0;omega<110*(int)omegakMax;omega++){
    //        std::cout<< "omega="<< omega<< std::endl;
			Lorentzian=epsilon/(pi*(pow(0.01*omega-omegak[k],2)+pow(epsilon,2)));
			for(int n1=-L;n1<L+1;n1++){
				for(int n2=-L;n2<L+1;n2++){
					S[omega][n1+L][n2+L]+=Sk[n1+L][n2+L]*Lorentzian;
				}
			}
		}
	}
	std::cout<<"Loop done"<<std::endl;
	
	auto end3 = std::chrono::high_resolution_clock::now();
	auto elapsed3 = std::chrono::duration_cast<std::chrono::nanoseconds>(end3 - begin3);
	std::cout<<"Time needed for the dynamical structure factor: "<< elapsed3.count()*1e-9<< " seconds"<< std::endl;
	
	std::ofstream result;
	int nn,n1,n2;
	Eigen::Vector2d q;
	result.open("SFdynamical.dat");
	for(int omega=0;omega<110*(int)omegakMax;omega++){
		nn=0;
		n1=-L;												
		n2=-L/2;						
		q=(n1*b1+n2*b2)/L;	
//        std::cout<<"momentum:"<< nn<< std::endl;
		result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< S[omega][n1+L][n2+L]<< std::endl;	
		while(n1<0){					
			n1+=2;						
			n2++;						
			q=(n1*b1+n2*b2)/L;
			nn++;
//            std::cout<<"momentum:"<< nn<< std::endl;
			result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< S[omega][n1+L][n2+L]<< std::endl;	
		}
		while(n2<L/2){
			n2++;
			q=(n1*b1+n2*b2)/L;
			nn++;
//            std::cout<<"momentum:"<< nn<< std::endl;
			result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< S[omega][n1+L][n2+L]<< std::endl;
		}
		while(n1<L){
			n1+=2;
			n2++;
			q=(n1*b1+n2*b2)/L;
			nn++;
//            std::cout<<"momentum:"<< nn<< std::endl;
			result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< S[omega][n1+L][n2+L]<< std::endl;
		}
		while(n1>0){
			n1--;
			n2--;
			q=(n1*b1+n2*b2)/L;
			nn++;
//            std::cout<<"momentum:"<< nn<< std::endl;
			result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< S[omega][n1+L][n2+L]<< std::endl;
		}
	}
	result.close();
	
	return S;
}

//Eigen::Matrix3d spinRotation(Eigen::Vector3d S){
//	double phi=azimuthal(S(0),S(1),S(2));
//	double theta=polar(S(0),S(1),S(2));
//	
//	Eigen::Matrix3d R;
//	R(0,0)=cos(phi)*cos(theta);
//	R(0,1)=-sin(phi)*cos(theta);
//	R(0,2)=sin(theta);
//	R(1,0)=sin(phi);
//	R(1,1)=cos(phi);
//	R(1,2)=0;
//	R(2,0)=-cos(phi)*sin(theta);
//	R(2,1)=sin(phi)*sin(theta);
//	R(2,2)=cos(theta);
//	
//	return R;
//}

//Eigen::Matrix3d spinRotation(Eigen::Vector3d S){
//	double phi=azimuthal(S(0),S(1),S(2));
//	double theta=polar(S(0),S(1),S(2));
//	
//	Eigen::Matrix3d R;
//	R(0,0)=cos(phi)*cos(phi)*cos(theta)+sin(phi)*sin(phi);
//	R(0,1)=0.5*sin(2*phi)*(cos(theta)-1);
//	R(0,2)=cos(phi)*sin(theta);
//	R(1,0)=0.5*sin(2*phi)*(cos(theta)-1);
//	R(1,1)=sin(phi)*sin(phi)*cos(theta)+cos(phi)*cos(phi);
//	R(1,2)=sin(phi)*sin(theta);
//	R(2,0)=-cos(phi)*sin(theta);
//	R(2,1)=-sin(phi)*sin(theta);
//	R(2,2)=cos(theta);
//	
//	return R;
//}

Eigen::Matrix3d spinRotation(Eigen::Vector3d S){
	double phi=azimuthal(S(0),S(1),S(2));
	double theta=polar(S(0),S(1),S(2));
	
	Eigen::Matrix3d R;
	R(0,0)=-sin(phi);
	R(0,1)=-cos(phi)*cos(theta);
	R(0,2)=cos(phi)*sin(theta);
	R(1,0)=cos(phi);
	R(1,1)=-sin(phi)*cos(theta);
	R(1,2)=sin(phi)*sin(theta);
	R(2,0)=0;
	R(2,1)=sin(theta);
	R(2,2)=cos(theta);
	
	return R;
}

//Eigen::Matrix3cd spinRotation(Eigen::Vector3d S){
//	double phi=azimuthal(S(0),S(1),S(2));
//	double theta=polar(S(0),S(1),S(2));
//	
//	Eigen::Matrix3cd R;
//	R(0,0)=0.5*(cos(phi)*cos(theta)+I*sin(phi));
//	R(0,1)=0.5*(cos(phi)*cos(theta)-I*sin(phi));
//	R(0,2)=cos(phi)*sin(theta);
//	R(1,0)=0.5*(sin(phi)*cos(theta)+I*cos(phi));
//	R(1,1)=0.5*(sin(phi)*cos(theta)-I*cos(phi));
//	R(1,2)=sin(phi)*sin(theta);
//	R(2,0)=-0.5*sin(theta);
//	R(2,1)=-0.5*sin(theta);
//	R(2,2)=cos(theta);
//	
//	return R;
//}

void printMatrices(Eigen::MatrixXcd A,Eigen::MatrixXcd B,int N){
	std::ofstream fileA;
	std::ofstream fileB;
	int i,j;
	
	fileA.open("matrixA.dat");
	fileB.open("matrixB.dat");
	
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			fileA<< A(i,j)<< "\t";
			fileB<< B(i,j)<< "\t";
		}
		fileA<< std::endl;
		fileB<< std::endl;
	}
	fileA.close();
	fileB.close();
}

Eigen::MatrixXcd bogoliubovStep1(Eigen::MatrixXcd A,Eigen::MatrixXcd B,int N){
	Eigen::MatrixXd Sigma(2*N,2*N);
	Sigma<< MatrixXd::Identity(N,N),
		MatrixXd::Zero(N,N),
		MatrixXd::Zero(N,N),
		-MatrixXd::Identity(N,N);
	
	Eigen::MatrixXcd M(2*N,2*N);
	M << 0.5*A, 0.5*B.conjugate().transpose(),
		0.5*B, 0.5*A.transpose();
//	std::cout<< M<< std::endl;	
	
//	STEP ONE
	Eigen::MatrixXcd eigendecomposition0(2*N+1,2*N);
	Eigen::MatrixXcd eigendecomposition;
	Eigen::VectorXd Lambda;
	Eigen::MatrixXcd Z;
	
	ComplexEigenSolver<MatrixXcd> diagonalization1;
	diagonalization1.compute(Sigma*M);
//	std::cout<< "Eigenvalues: \n"<< diagonalization1.eigenvalues()<<std::endl; //"\n Eigenvectors: \n"<< diagonalization1.eigenvectors()<< std::endl;
	eigendecomposition0.topRows(1)=diagonalization1.eigenvalues().transpose();		//eigenvalues on the first row
	eigendecomposition0.bottomRows(2*N)=diagonalization1.eigenvectors();			//eigenvectors underneath
//	Lambda=diagonalization1.eigenvalues().real();
	eigendecomposition=sortEigenvalues1(eigendecomposition0,N);
	
	return eigendecomposition;
}

std::vector<std::pair<double,int>> findDegeneracies(Eigen::VectorXd Lambda){		
	std::vector <double>  uniqueEigenvalues;
	std::vector <int> degeneracy;
	
	double uniqueValue=Lambda(0);
	uniqueEigenvalues.assign(1,uniqueValue);
	int j=0;
	int dj=1;
	for(int i=1;i<Lambda.size();i++){
		if(abs(Lambda(i)-uniqueValue)<0.001){
			dj++;
		}
		else{
			degeneracy.push_back(dj);
			uniqueEigenvalues.push_back(Lambda(i));
			uniqueValue=Lambda(i);
			j++;
			dj=1;
		}	
	}
	degeneracy.push_back(dj);
	int m=j+1;
	
//	std:: cout<< "\n";
//	for(int i=0;i<uniqueEigenvalues.size();i++)
//		std::cout<< "Eigenvalue "<< uniqueEigenvalues[i]<< " with multiplicity "<< degeneracy[i]<< std:: endl; 
//	std:: cout<< "Number of unique eigenvalues "<< m<< "\n" << std::endl;
	
	std::vector<std::pair<double,int>> deg;
	for(int i=0;i<m;i++)
		deg.push_back(std::make_pair(uniqueEigenvalues[i],degeneracy[i]));
		
	return deg;
	
//	std::ofstream DoS;
//	DoS.open("densityOfStates.dat");
//	for(int i=0;i<m/2;i++)
//		DoS<<std::setprecision(16)<< uniqueEigenvalues[i]<< "\t"<< degeneracies[i]<< std::endl;
//	DoS.close();
}
	
Eigen::MatrixXcd bogoliubovStep2(Eigen::MatrixXcd Z,std::vector<std::pair<double,int>> deg){
	Eigen::ComplexEigenSolver<MatrixXcd> diagonalization2;
	int m=deg.size();
//	std::cout<< "deg="<< std::endl;
//	for(int i=0;i<m;i++)
//		std::cout<< deg[i].first<< "\t"<< deg[i].second<< std::endl; 
	int n=Z.cols()/2;
	int N=Z.rows()/2;
	Eigen::MatrixXcd ZBlock[m];
	Eigen::MatrixXcd LBlock[m];
	Eigen::MatrixXcd UBlock[m];
	Eigen::ArrayXcd lBlock[m];
	Eigen::MatrixXcd lMinusHalfBlock[m];
	Eigen::MatrixXcd lMinusHalf;
	Eigen::MatrixXcd U;
	Eigen::MatrixXcd T;
	Eigen::MatrixXd Sigma(2*N,2*N);
	Sigma<< MatrixXd::Identity(N,N),
		MatrixXd::Zero(N,N),
		MatrixXd::Zero(N,N),
		-MatrixXd::Identity(N,N);
	std::vector<int> degeneracies(m);
	for(int i=0;i<m;i++)
		degeneracies[i]=deg[i].second;	//first element of deg is the eigenvalue and second the degree of degeneracy

	int d=0;
	for(int j=0;j<m;j++){	
		ZBlock[j]=Z.block(0,d,2*N,degeneracies[j]);
		LBlock[j]=ZBlock[j].adjoint()*Sigma*ZBlock[j];	//dxd
		diagonalization2.compute(LBlock[j]);
		UBlock[j]=diagonalization2.eigenvectors();	//dxd
		lBlock[j]=diagonalization2.eigenvalues().array();
		lMinusHalfBlock[j]=(lBlock[j].abs().rsqrt()).matrix().asDiagonal();	//rsqrt = reciprocal sqrt
		d+=degeneracies[j];
		
//		std::cout<< "j="<<j<<std::endl;
//		std::cout<< "ZBlock["<<j<< "]=\n"<< ZBlock[j]<< std::endl;
//		std::cout<<"LBlock["<< j<< "]="<< LBlock[j]<< std::endl;
//		std::cout<<"lBlock["<< j<< "]="<< lBlock[j]<< std::endl;
//		std::cout<< "lMinusHalf["<< j<< "]=\n"<< lMinusHalfBlock[j]<< std::endl;
	}	

	lMinusHalf=MatrixXd::Zero(2*n,2*n);
	U=MatrixXd::Zero(2*n,2*n);
	d=0;
	for(int j=0;j<m;j++){
		lMinusHalf.block(d,d,degeneracies[j],degeneracies[j])=lMinusHalfBlock[j];
		U.block(d,d,degeneracies[j],degeneracies[j])=UBlock[j];
		d+=degeneracies[j];
	}
	
	
	T=Z*U*lMinusHalf;
	
//	std::cout<< "U=\n"<< U<< "\n"<< std::endl;
//	std::cout<< "l=\n"<< lMinusHalf<< "\n"<< std::endl;
//	std::cout<<"T=\n"<< T<< "\n"<< std::endl;
//	std::cout<< "T_dagger*Sigma*T=\n"<< T.adjoint()*Sigma*T<< "\n"<< std::endl;
//	std::cout<< "T_dagger*M*T="<< T.adjoint()*M*T<< "\n"<<std::endl;
	
	
//	printEigendecomposition(Lambda,Z,/*T,*/N,n);
		
	//std::complex <double> lamda[2*N]=eigenDecomp.eigenvalues();
	
	return T;
}

bool compareAscending (std::pair<double, Eigen::VectorXcd> x,std::pair<double, Eigen::VectorXcd> y) 	//compare pairs wrt first element in ascending order
    { return x.first < y.first;}
    
bool compareDescending (std::pair<double, Eigen::VectorXcd> x,std::pair<double, Eigen::VectorXcd> y)	//compare pairs wrt first element in descending order
    { return x.first > y.first;}

Eigen::MatrixXcd sortEigenvalues1(Eigen::MatrixXcd LV, int N)	//sort eigenvalues
{
	Eigen::VectorXd lambda=LV.topRows(1).transpose().real();
//	std::cout<< lambda<< std::endl;
	Eigen::MatrixXcd V=LV.bottomRows(2*N);
	int numZeroModes=0;
	int nn=2*N-numZeroModes;
	double test;
	for(int i=0;i<nn;i++){	//After removing elemnert i, i+1 now becomes i, therefore I am skipping it. Fix it!
//		std::cout<< "i="<< i<< "\t"<< llambda(i)<< std::endl;
		test=abs(lambda(i));
		if(test<0.001){
			removeElement(lambda,i);
			removeColumn(V,i);
			numZeroModes++;
			nn=2*N-numZeroModes;
			i--;
			std::cout<< lambda<< std::endl;
		}
	}
	int n=nn/2;
	
//	std::cout<< "System size "<< N<< ", non-zero eigenvalues "<<n << std::endl;
    std::pair<double, Eigen::VectorXcd> positives[n];
    std::pair<double, Eigen::VectorXcd> negatives[n];
  
  	int j1=0;
  	int j2=0;
    for (int i = 0; i < 2*n; i++) 
    {
    	if(lambda(i)>0){
    		positives[j1].first = lambda(i);
    		positives[j1].second = V.col(i);
    		j1++;
		}
		else if(lambda(i)<0)
		{
			negatives[j2].first = lambda(i);
        	negatives[j2].second = V.col(i);			
        	j2++;
		}
    }
        
	// Sorting the pair array.
    std::sort(positives,positives+n,compareDescending);
    std::sort(negatives,negatives+n,compareAscending);
        
    Eigen::VectorXd eigenvalues(2*n);
    Eigen::MatrixXcd eigenvectors(2*N,2*n);
    
    for(int i=0;i<n;i++){
    	eigenvalues(i)=positives[i].first;
    	eigenvectors.col(i)=positives[i].second;
    	eigenvalues(i+n)=negatives[i].first;
    	eigenvectors.col(i+n)=negatives[i].second;
 //   	std::cout<<eigenvalues(i)<< "\t"<< eigenvalues(i+n)<< std::endl;
	}
	
	Eigen::MatrixXcd eigenset(2*N+1,2*n);
	eigenset.topRows(1)=eigenvalues.transpose();
	eigenset.bottomRows(2*N)=eigenvectors;
      
  	return eigenset;
}
/*
std::pair<Eigen::VectorXd,Eigen::MatrixXcd> sortEigenvalues(Eigen::VectorXcd lambda, Eigen::MatrixXcd V, int N)	//sort eigenvalues
{
//	int nn=0; 
//	Eigen::VectorXd lambda;
//	Eigen::MatrixXcd V;
	int numZeroModes=0;
	double test;
	for(int i=0;i<2*N;i++){
//	std::cout<< "i="<< i<< "\t"<< llambda(i)<< std::endl;
	test=abs(real(lambda(i)));
//		if(test>0.0001){
//			lambda(nn)=real(llambda(i));
//			V.col(nn)=VV.col(i);
//			nn++;
//		}
		if(test<0.0001){
			removeElement(lambda,i);
			removeColumn(V,i);
			numZeroModes++;
		}
	}
	
	int n=N-numZeroModes/2;
//	std::cout<< "System size "<< N<< ", non-zero eigenvalues "<<n << std::endl;
    std::pair<double, Eigen::VectorXcd> positives[n];
    std::pair<double, Eigen::VectorXcd> negatives[n];
  
  	int j1=0;
  	int j2=0;
    for (int i = 0; i < 2*n; i++) 
    {
    	if(real(lambda(i))>0){
    		positives[j1].first = real(lambda(i));
    		positives[j1].second = V.col(i);
    		j1++;
		}
		else if(real(lambda(i))<0)
		{
			negatives[j2].first = real(lambda(i));
        	negatives[j2].second = V.col(i);			
        	j2++;
		}
    }
        
	// Sorting the pair array.
    std::sort(positives,positives+n,compareDescending);
    std::sort(negatives,negatives+n,compareAscending);
        
    Eigen::VectorXd eigenvalues(2*n);
    Eigen::MatrixXcd eigenvectors(2*N,2*n);
    
    for(int i=0;i<n;i++){
    	eigenvalues(i)=positives[i].first;
    	eigenvectors.col(i)=positives[i].second;
    	eigenvalues(i+n)=negatives[i].first;
    	eigenvectors.col(i+n)=negatives[i].second;
 //   	std::cout<<eigenvalues(i)<< "\t"<< eigenvalues(i+n)<< std::endl;
	}
	
	std::pair<Eigen::VectorXd,Eigen::MatrixXcd> eigenset;
	eigenset.first=eigenvalues;
	eigenset.second=eigenvectors;
      
  	return eigenset;
}*/

void removeColumn(Eigen::MatrixXcd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

void removeElement(Eigen::VectorXd& vector, unsigned int elToRemove)
{
    unsigned int numElements = vector.size()-1;

    if( elToRemove < numElements )
        vector.segment(elToRemove,numElements-elToRemove) = vector.tail(numElements-elToRemove);

    vector.conservativeResize(numElements);
}

void printEigendecomposition(Eigen::VectorXd eigenvalues,Eigen::MatrixXcd eigenvectorsZ,/*Eigen::MatrixXcd eigenvectorsT,*/int N,int n){
	std::ofstream file1;
	std::ofstream file2;
	std::ofstream file3;
//	std::ofstream file4;
	int i,j;
	
	file1.open("eigenvalues.dat");
	file2.open("spectrum.dat");
	file3.open("eigenvectorsZ.dat");
//	file4.open("eigenvectorsT.dat");
	
	for(i=0;i<2*n;i++){
		file1<< std::setprecision(16)<< eigenvalues(i)<< std::endl;
		if(i<n)
			file2<< std::setprecision(16)<< eigenvalues(i)-eigenvalues(i+n)<< std::endl;
	}
	
	for(i=0;i<2*N;i++){
		for(j=0;j<2*n;j++){
			file3<< eigenvectorsZ(i,j)<< "\t";
//			file4<< eigenvectorsT(i,j)<< "\t";
		}
		file3<< std::endl;
//		file4<< std::endl;
	}
	
	file1.close();
	file2.close();
	file3.close();
//	file4.close();
}

std::vector<double> densityOfStates(std::vector<double> spectrum){
	double omegakMax=*max_element(spectrum.begin(),spectrum.end());
	int N=spectrum.size();
	double epsilon=0.1;
	std::vector<double> dos;
	dos.assign(110*(int)omegakMax,0);
	
	for(int k=0;k<N;k++){
		for(int omega=0;omega<110*(int)omegakMax;omega++){
			dos[omega]+=epsilon/(pi*(pow(0.01*omega-spectrum[k],2)+pow(epsilon,2)));
		}
	}
	std::ofstream result;
	result.open("DOS.dat");
	for(int omega=0;omega<110*(int)omegakMax;omega++)
		result<< 0.01*omega<< "\t"<< dos[omega]<< std::endl;
	result.close();
	
	return dos;
}

void densityOfStates0(Eigen::VectorXd Lambda,int n){
	int i,j;	
	std::vector<double> energies(n);

	for(i=0;i<n;i++)
		energies[i]=Lambda(i);
		
	double max=*max_element(energies.begin(),energies.end());
	double min=*min_element(energies.begin(),energies.end());
	int intervals=round(sqrt(n));
	if(intervals%2==0)
		intervals+=1;
//.	int intervals=n;
	double interval=(max-min)/intervals;
	
	std::vector<double> densityOfStates(intervals);
	densityOfStates.assign(intervals,0);
	
	double omega;
	for(i=0;i<n;i++){	
		for(j=0;j<intervals;j++){
			omega=min+(j+0.5)*interval;
			if(abs(energies[i]-omega)<=0.5*interval){
				densityOfStates[j]++;
				break;
			}
		}
	}
	
	std::ofstream file;
	file.open("densityOfStates.dat");
	
	for(j=0;j<intervals;j++){
		omega=min+(j+0.5)*interval;
		file<< omega<< "\t"<< densityOfStates[j]/n<< std::endl;
	}
	file.close();
}

bool sortComplex(std::complex<double> a, std::complex<double> b)
{
    if (real(a) == real(b))
        return imag(a) < imag(b);
    return real(a) < real(b);
}

void averageDensityOfStates(std::vector<std::vector<double>> doss,int n){
	int omegaMax;
	std::vector<int> sizes;
	for(int i=0;i<n;i++)
		sizes.push_back(doss[i].size());	
	omegaMax=*min_element(sizes.begin(),sizes.end());
	
	double sum;
	std::vector<double> averageDOS;
	for(int omega=0;omega<omegaMax;omega++){
		sum=0;
		for(int i=0;i<n;i++){
			sum+=doss[i][omega];
		}
		averageDOS.push_back(sum/(1.0*n));
	}
	
	std::ofstream file;
	file.open("averageDOS_h=2,22.dat");
	for(int omega=0;omega<omegaMax;omega++)
		file<< 0.01*omega<< "\t"<< averageDOS[omega]<< std::endl;
	file.close();
}

void averageMagnetization(std::vector<Eigen::Vector3d> magnetizations,int numberOfImpurityConfigurations){
	Eigen::Vector3d averageM(0,0,0);
	for(int i=0;i<numberOfImpurityConfigurations;i++)
		averageM+=magnetizations[i];
	averageM/=(1.0*numberOfImpurityConfigurations);
	
	Eigen::Vector3d a,b,c;
    a<< 1/sqrt(6.),1/sqrt(6.),-2/sqrt(6.);
    b<< -1/sqrt(2.),1/sqrt(2.),0;
    c<< 1/sqrt(3.),1/sqrt(3.),1/sqrt(3.);
    
    std::ofstream results;
    results.open("averageMagnetization.dat");
    results<< std::setprecision(12) << averageM.dot(a)<< "\t"<< averageM.dot(b)<< "\t" << averageM.dot(c)<<std::endl;
    results.close();
}

void averageStructureFactorStatic(std::vector<std::vector<std::vector<double>>> SSFs,int n,int L){
	std::vector<std::vector<double>> averageSSF;
	std::vector<double> averageSSF2;
	averageSSF2.assign(2*L+1,0);
	averageSSF.assign(2*L+1,averageSSF2);
	double sum;
	
	for(int n1=-L;n1<L+1;n1++){
		for(int n2=-L;n2<L+1;n2++){
			sum=0;
			for(int i=0;i<n;i++){
				sum+=SSFs[i][n1+L][n2+L];
			}
			averageSSF[n1+L][n2+L]+=sum/(1.0*n);
		}
	}
	
	Eigen::Vector2d b1(2*pi/sqrt(3.),-2*pi/3.);
	Eigen::Vector2d b2(0,4*pi/3.);
	std::ofstream results;
	results.open("averageSFstatic_h=2,22.dat");
	Eigen::Vector2d q;
	for(int n1=-L;n1<=L;n1++){
		for(int n2=-L;n2<=L;n2++){
			if(n1-n2<=L&&n1-n2>=-L){
				q=(n1*b1+n2*b2)/L;
				results<< q(0)<< "\t"<< q(1)<< "\t"<< averageSSF[n1+L][n2+L]<< std::endl;
			}
		}
	}
	results.close();
}

void averageStructureFactorDynamical(std::vector<std::vector<std::vector<std::vector<double>>>> DSFs,int n,int L){
	
	int omegaMax;
	std::vector<int> sizes;
	for(int i=0;i<n;i++)
		sizes.push_back(DSFs[i].size());	
	omegaMax=*min_element(sizes.begin(),sizes.end());
	
	std::vector<std::vector<std::vector<double>>> averageDSF;
	std::vector<std::vector<double>> averageDSF1;
	std::vector<double> averageDSF2;
	averageDSF2.assign(2*L+1,0);
	averageDSF1.assign(2*L+1,averageDSF2);
	averageDSF.assign(omegaMax,averageDSF1);
	
	double sum;
	for(int omega=0;omega<omegaMax;omega++){
		for(int n1=-L;n1<L+1;n1++){
			for(int n2=-L;n2<L+1;n2++){
				sum=0;
				for(int i=0;i<n;i++){
					sum+=DSFs[i][omega][n1+L][n2+L];
				}
				averageDSF[omega][n1+L][n2+L]+=sum/(1.0*n);
			}
		}
	}
	
	Eigen::Vector2d b1(2*pi/sqrt(3.),-2*pi/3.);
	Eigen::Vector2d b2(0,4*pi/3.);
	std::ofstream result;
	int nn,n1,n2;
	Eigen::Vector2d q;
	result.open("averageSFdynamical_h=2,22.dat");
	for(int omega=0;omega<omegaMax;omega++){
		nn=0;
		n1=-L;												
		n2=-L/2;						
		q=(n1*b1+n2*b2)/L;	
//        std::cout<<"momentum:"<< nn<< std::endl;
		result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< averageDSF[omega][n1+L][n2+L]<< std::endl;	
		while(n1<0){					
			n1+=2;						
			n2++;						
			q=(n1*b1+n2*b2)/L;
			nn++;
//            std::cout<<"momentum:"<< nn<< std::endl;
			result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< averageDSF[omega][n1+L][n2+L]<< std::endl;	
		}
		while(n2<L/2){
			n2++;
			q=(n1*b1+n2*b2)/L;
			nn++;
//            std::cout<<"momentum:"<< nn<< std::endl;
			result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< averageDSF[omega][n1+L][n2+L]<< std::endl;
		}
		while(n1<L){
			n1+=2;
			n2++;
			q=(n1*b1+n2*b2)/L;
			nn++;
//            std::cout<<"momentum:"<< nn<< std::endl;
			result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< averageDSF[omega][n1+L][n2+L]<< std::endl;
		}
		while(n1>0){
			n1--;
			n2--;
			q=(n1*b1+n2*b2)/L;
			nn++;
//            std::cout<<"momentum:"<< nn<< std::endl;
			result<< nn<< "\t"<< q(0)<< "\t"<< q(1)<< "\t"<< 0.01*omega<< "\t"<< averageDSF[omega][n1+L][n2+L]<< std::endl;
		}
	}
	result.close();
}
