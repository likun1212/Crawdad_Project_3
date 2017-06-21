//The mass header file

#include "mass.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <armadillo>
#include <complex>

#include <array>   // to be able to use this library compile using this g++ -std=c++0x FileNames


int index1(int i, int j);
int index_call (int i, int j, int k, int l);
using namespace std;

int main()
{


/**************************************************************/	
//Reading nuclear repulsion energy from the file 
/**************************************************************/	


	double enuc=0.0;
	ifstream enuc1 ("enuc.dat", ios::in);
	if (enuc1.is_open())
	{
		 cout<<"\nReading from the nuclear repulsion energy file\n";
		 enuc1>>enuc;
		 cout<<"\nThe nuclear repulsion energy is:\n"<<enuc<<endl;
	}

/**************************************************************/	
//Reading the one electron integrals 
/**************************************************************/	

//Reading the overlap matrix elements
	arma::mat S_file;
	S_file.load("s.dat");
	int last_s=S_file(S_file.n_rows-1,0);
	arma::mat S=arma::zeros(last_s,last_s);
	for(int i=0; i<S_file.n_rows; i++)
	{
		int i1=S_file(i,0);
		int i2=S_file(i,1);
		double val=S_file(i,2);
		S(i1-1,i2-1)=val;
		S(i2-1,i1-1)=val;
	}	
	S.print("\nThe overlap intergral elements are\n ");

	arma::mat T_file;
	T_file.load("t.dat");
	int last_t=T_file(T_file.n_rows-1,0);
	arma::mat T=arma::zeros(last_t,last_t);
	for(int i=0; i<T_file.n_rows; i++)
	{
		int i1=T_file(i,0);
		int i2=T_file(i,1);
		double val=T_file(i,2);
		T(i1-1,i2-1)=val;
		T(i2-1,i1-1)=val;
	}	
	T.print("\nThe kinetic energy intergals are\n");

	arma::mat V_file;
	V_file.load("v.dat");
	int last_v=V_file(V_file.n_rows-1,0);
	arma::mat V=arma::zeros(last_v,last_v);
	for(int i=0; i<V_file.n_rows; i++)
	{
		int i1=V_file(i,0);
		int i2=V_file(i,1);
		double val=V_file(i,2);
		V(i1-1,i2-1)=val;
		V(i2-1,i1-1)=val;
	}	
	V.print("\nThe nuclear attaraction intergral elements are\n ");

	arma::mat H_c=T+V;
	real(H_c).print("\nThis is the core Hamiltonian\n");

/**************************************************************/	
//Reading the two electron integrals 
/**************************************************************/	

//here we have to store the 4D matrix involving the 2 electron integrals
//into a one dimensional matrix
	arma::mat d_file;
	d_file.load("eri.dat");
	arma::vec d_I =arma::zeros(1000);	
	for(int i1=0; i1<d_file.n_rows; i1++)
	{
		int i=d_file(i1,0);
		int j=d_file(i1,1);
		int k=d_file(i1,2); 
		int l=d_file(i1,3);
		double val=d_file(i1,4);
		int ij, kl,ijkl;
		if(i>j)
		{
			 ij=index1(i,j);
		}
		else    
		{
			 ij=index1(j,i);
		}
		if (k>l)
		{
			 kl=index1(k,l);
		}
		else
		{
			kl=index1(l,k);
		}
		if(ij > kl)
		{
			 ijkl=ij*(ij+1)/2+(kl);
		}
		else //if (ij < kl)

		{
			 ijkl=kl*(kl+1)/2+(ij);
		}
		d_I(ijkl)=val;
	}	


/**************************************************************/
//Building the orthogonalization matrix
/**************************************************************/	

	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym (eigval,eigvec,S);
	real(eigvec).print("\nPrinting the eigenvectors of the overlap matrix\n");
	real(eigval).print("\nPrinting the eigenvalues of the overlap matrix\n");
	arma::mat tran_eigvec=eigvec.t();
	arma::mat prod=tran_eigvec*eigvec;
	real(prod).print("\nThe product of the eigenvec and its transpose\n");
	arma::mat sqrt_eigen=arma::pow(eigval,-0.5);
	real(sqrt_eigen).print("\nThe square root of the eigenvalue matrix\n");
	arma::mat diag_ee=arma::diagmat(sqrt_eigen);
	real(diag_ee).print("\nThe diagonal matrix of the eigenvalues\n");
	arma::mat S_sqrt=eigvec*diag_ee*tran_eigvec;
	arma::mat r_S_sqrt=arma::real(S_sqrt);
	real(S_sqrt).print("\nThe orthogonalized symmetric overlap matrix\n");


/**************************************************************/
//Building the initial guess density matrix 
/**************************************************************/	

	arma::mat I_Gs=S_sqrt.t()*H_c*S_sqrt;
	I_Gs.print("\nThis is the initial Fock matrix(orthogonal basis)\n");
//Note that this contains the inital orbital range
	arma::vec eig11; 
	arma::mat eigenvec;
	arma::eig_sym(eig11,eigenvec,I_Gs);
	real(eig11).print("\nThe eigenvalues (initially)\n");
	real(eigenvec).print("\nThe eigenvectors (initially)\n");
	//eigvec1.print("\nThe initial MO coefficient Matrix is \n");
//Transforming the eigenvectors into the original basis
	arma::mat eig_AO = S_sqrt*eigenvec;
	real(eig_AO).print("\nThe initial MO cofficients\n");
//The density matrix
//Since only the first MO of water are occupied so the lowest five 
//MO are taken into consideration and they are effectively squared.	
	eig_AO.shed_cols(5,6);
	arma::mat den_I = eig_AO*eig_AO.t();	
	real(den_I).print("\nThe initial guess density matrix\n");

/**************************************************************/
//Computing the initial SCF energy matrix and the Fock Matrix
/**************************************************************/	

//Computing the SCF electronic energy using the density matrix

	double energy=0.0;
	for(int i=0; i< H_c.n_rows; i++)
	{
		for(int j=0; j<H_c.n_cols; j++)
		{
			energy+=den_I(i,j)*(H_c(i,j)+I_Gs(i,j));
		}
	}
	
	printf("\n %s\t %s\t %s\t %s\t\t %s\n %d %19.5f \n","Iteration", "Electronic energy","E_total","Delta(E)","RMS(D)",0,energy);

/**************************************************************/
//Computing the new Fock Matrix and starting the SCF loop
/**************************************************************/	
	arma::mat density_initial=den_I;
	double energy_last=energy;
//setting the precision here
	cout.precision(15);
	arma::mat F_mv=H_c;
        arma::mat coeff_car=den_I;	
//the SCF loop starts here
	for(int ii=1; ii<100000; ii++)
	{	
		F_mv=H_c;
		for(int i=0;i<H_c.n_rows;i++)
		{
			for(int j=0; j<H_c.n_cols; j++)
			{
				for(int k=0;k<H_c.n_cols;k++)
				{
					for(int l=0;l<H_c.n_cols;l++)
					{
						F_mv(i,j)+=den_I(k,l)*(2*d_I(index_call(i+1,j+1,k+1,l+1))-d_I(index_call(i+1,k+1,j+1,l+1)));
					}
				}
			}
		}
		//F_mv.print("\nThe Fock matrix with 2 electron integrals \n");
		arma::mat F_store=S_sqrt.t()*F_mv*S_sqrt;
		arma::vec eig_store;
		arma::mat eigenvec_store;
		arma::eig_sym(eig_store,eigenvec_store,F_store);
		arma:: mat eig_AO_store=S_sqrt*eigenvec_store;
		coeff_car=eig_AO_store;
		//eig_AO_store.print("\n");
		eig_AO_store.shed_cols(5,6);
		arma:: mat density_store=eig_AO_store*eig_AO_store.t();
		double energy1=0.0;
		for(int i=0; i< H_c.n_rows; i++)
		{
			for(int j=0; j<H_c.n_cols; j++)
			{
				energy1+=density_store(i,j)*(H_c(i,j)+F_mv(i,j));
			}
		}
		double enery_last=energy1;
		double d_ene=energy_last-energy1;
		double rms=0.0;
		double rms_last=0.0;
		for(int i=0; i<density_store.n_rows;i++)
		{
			for (int j=0;j<density_store.n_cols;j++)
			{
				rms+=pow(density_store(i,j)-density_initial(i,j),2);
			}
		}
		rms=pow(rms,0.5);
		printf("\n %d\t  %10.12f \t %10.12f \t %10.12f \t %10.12f \n ",ii,energy1,energy-enuc,d_ene,rms);	
		//F_mv.print("\n");
		if (abs(d_ene)<1e-12 && abs(rms_last-rms)< 1e-12 )
			break;
//Changing the stuff here
		density_initial=density_store;
		energy_last=energy1;
		rms_last=rms;
		den_I=density_store;
	}		
/**************************************************************/
//The MO-basis Fock-Fatrix
/**************************************************************/
	//F_mv.print("\n");
	//coeff_car.print("\n");	
	arma::mat F_MO2=coeff_car.t()*F_mv*coeff_car;
	cout<<"\n\nChecking for the Diagonal form of Fork Matrix in the basis of the Molecular Orbitals\n\n";

	for (int i=0; i<F_MO2.n_rows;i++)
	{
		for(int j=0; j<F_MO2.n_cols;j++)
		{
			printf("%15.7f",F_MO2(i,j));
		}
		cout<<endl;
	}
			
/********************************************************************/
//One Electron Properties- electronic contribution to dipole operator
/*********************************************************************/
	
}



int index1(int i, int j)
{
	int ij=i*(i+1)/2+j;
	return ij;	
}


int index_call(int i, int j , int k, int l)
{
	int ij,kl,ijkl;
	if(i>j)
	{
		 ij=index1(i,j);
	}
	else    
	{
		 ij=index1(j,i);
	}
	if (k>l)
	{
		 kl=index1(k,l);
	}
	else
	{
		kl=index1(l,k);
	}
	if(ij > kl)
	{
		 ijkl=ij*(ij+1)/2+(kl);
	}
	else //if (ij < kl)

	{
		 ijkl=kl*(kl+1)/2+(ij);
	}

	return ijkl;
}








