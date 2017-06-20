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

int index1(int i,int j);
int index2(int i,int j);
int index_call(int i, int j , int k, int l);

using namespace std;

int main()
{
	arma::mat d_file;
	d_file.load("eri.dat");
	arma::vec index(d_file.n_rows);
	arma::vec values(1000);
	for(int i1=0; i1<d_file.n_rows; i1++)
	{
		int i=d_file(i1,0);
		int j=d_file(i1,1);
		int k=d_file(i1,2); 
		int l=d_file(i1,3);
		double val=d_file(i1,4);
		int ij, kl;
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
		int ijkl;
		cout<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<ij<<" "<<kl<<" "<<endl;
		if(ij > kl)
		{
			ijkl=ij*(ij+1)/2+(kl);
		}
		else 
		{
			ijkl=kl*(kl+1)/2+(ij);
		}
		index(i1)=ijkl;
		values(ijkl)=val;
	}
	arma::vec sindex=arma::sort(index);
	for(int i=0; i<sindex.n_rows; i++)
	{
		cout<<int(sindex(i))<<endl;
	}
	for(int i=0; i<values.n_rows; i++)
	{
		printf("\n %d %8.5f\n ",i,values(i));
	}

	for(int i=0;i<7;i++)
	{
		for(int j=0; j<7; j++)
		{
			for(int k=0;k<7;k++)
			{
				for(int l=0;l<7;l++)
				{
					printf("\n %d %d %d %d %s %d %8.5f %d %8.5f \n",i+1,j+1,k+1,l+1,"--",index_call(i+1,j+1,k+1,l+1),values(index_call(i+1,j+1,k+1,l+1)),index_call(i+1,k+1,j+1,l+1),values(index_call(i+1,k+1,j+1,l+1)));
				}
			}
		}
	}



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













