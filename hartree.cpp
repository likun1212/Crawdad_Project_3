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
	//cout<<S_file(27,0);
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
	//cout<<S_file(27,0);
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
	H_c.print("\nThis is the core Hamiltonian\n");

/**************************************************************/	
//Reading the two electron integrals 
/**************************************************************/	

//here we have to store the 4D matrix involving the 2 electron integrals
//into a one dimensional matrix
	arma::mat d_file;
	d_file.load("eri.dat");
	arma::vec d_I =arma::zeros();	
	for(int i1=0; i1<d_file.n_rows; i1++)
	{
		int i=d_file(i1,0);
		int j=d_file(i1,1);
		int k=d_file(i1,2); 
		int l=d_file(i1,3);
		int ij=i*10+j;
		int kl=k*10=l;
		double val=d_file(i1,4);
		int ijkl=0;
		if(ij > kl)
		{
			 ijkl=ij*(ij+1)/2+(kl);
		}
		else 
		{
			 ijkl=kl*(kl+1)/2+(ij);
		}
		d_I(ijkl)=val;
	}	
	//d_I.print("\nThe 2 electron  intergral elements are\n ");
	for(int i=0; i<d_I.n_rows; i++)
	{
		cout<<i<<" "<<d_I(i)<<endl;
	}






/**************************************************************/
//Building the orthogonalization matrix
/**************************************************************/	






/**************************************************************/
//Building the initial guess density matrix 
/**************************************************************/	





}
