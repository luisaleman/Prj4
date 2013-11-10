#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "time.h"
using namespace std;
/*I declare prototype to avoid the compiler believing that the functions declared
at the end don't exist*/
void expliciteuler(double tstep, double step, int npoints, int timesteps);
void impliciteuler(double tstep, double step, int npoints, int timesteps);
void cranknicolson(double tstep, double step, int npoints, int timesteps);
void tridiag(double a, double b, int npoints, double *sourcevector, double *resultvector);
void output(string name, int npoints, double step, double tstep, int timesteps, double *v);

//main program
int main(int argc, char* argv[]){
	double  tstep, step, total_time;
	int timesteps,npoints;
	clock_t start, end;
	cout<< endl <<"Please choose the time step: ";
	cin>> tstep;
	cout<< "Please choose the number of time steps: ";
	cin>> timesteps;
	cout<< "Please choose the number of grid points: ";
	cin>> npoints;
	step=1/(double)(npoints+1);
	cout <<endl<< "Step : " << step;
	
	expliciteuler(tstep, step, npoints, timesteps);

	impliciteuler(tstep, step, npoints, timesteps);

	cranknicolson(tstep, step, npoints, timesteps);

	return 0;
}


void expliciteuler(double tstep, double step, int npoints, int timesteps)
{
	clock_t start, end;
	start=clock();
		double *v;
		v=new double[npoints+2];
		double *vnew;
		vnew = new double[npoints+2];
		double alpha= tstep/(pow(step,2));
		double x;
	//setting the initial condition
		v[0]=0.00;
		vnew[0]=0.00;
	for(int i=1; i<npoints+2;i++){
		x=i*step;
		v[i]= x-1.00;
		vnew[i]=0.00;
	}
	//implementing the computation
		for(int t=1; t<= timesteps; t++){
			//computing the new values for the next time step
			for(int i=1; i<=npoints; i++){
				vnew[i]=alpha*v[i-1]+(1-2*alpha)*v[i]+alpha*v[i+1];
			}
			//update the solution without changing the boundary points
			for(int i=1; i<=npoints; i++){
				v[i]=vnew[i];
			}
		}
	cout << endl << "Alpha: " << alpha;
	end=clock();
	double total_time;
	total_time=((end-start)/(double)CLOCKS_PER_SEC);
    cout << endl << "The elapsed time for explicit euler was " << total_time<< " seg"<< endl;

    output("explicit", npoints, step, tstep, timesteps,v);
}

void impliciteuler(double tstep, double step,int npoints, int timesteps){
	clock_t start, end;
	start=clock();
	//set the mtrix constants
	double alpha= tstep/(pow(step,2));
	double a= 1+2*alpha;
	double b= -alpha;
	//set the vectors
	double *v;
	v=new double[npoints+2];
	double *vpast;
	vpast= new double[npoints+2];
	double x;
	//initial and boundary conditions
	v[0]=vpast[0]=0.00;
	for(int i=1;i<npoints+2;i++){
		x=i*step;
		vpast[i]=x-1;
		v[i]=0.00;
	}
	//implementing the computation
	for(int t=1; t<=timesteps; t++){
		tridiag(a,b,npoints, vpast, v);
		//replace the previous time solution with the new one
		for(int i=0;i<npoints+2;i++){
			vpast[i]=v[i];
		}
	}
	end=clock();
	double total_time;
	total_time=((end-start)/(double)CLOCKS_PER_SEC);
    cout << endl << "The elapsed time for implicit euler was " << total_time<< " seg"<< endl;

    output("implicit", npoints, step, tstep, timesteps, v);
}

void cranknicolson(double tstep, double step, int npoints, int timesteps){
	clock_t start, end;
	start=clock();
	//set the mtrix constants
	double alpha= tstep/(pow(step,2));
	double a= 2+2*alpha;
	double b= -alpha;
	//set the vectors
	double *v;
	v=new double[npoints+2];
	double *r;
	r= new double[npoints+2];
	double x;
	//initial and boundary conditions
	v[0]=r[0]=0.00;
	for(int i=1;i<npoints+2;i++){
		x=i*step;
		r[i]=0.00;
		v[i]=x-1;
	}
	//implementing the computation
	for(int t=1; t<=timesteps;t++){
		//Explicit part: computing the values that will be used in the implicit part
		for(int i= 1; i<=npoints; i++){
			r[i]=alpha*v[i-1]+(2-2*alpha)*v[i]+alpha*v[i+1];
		}
		//Implicit part: trigiagonal matrix
		tridiag(a,b, npoints,r,v);
	}
	end=clock();
	double total_time;
	total_time=((end-start)/(double)CLOCKS_PER_SEC);
	cout << endl << "The elapsed time for cranknicolson euler was " << total_time<< " seg"<< endl;

	output("cranknicolson", npoints, step, tstep, timesteps,v);
}

void tridiag(double a, double b, int npoints, double *sourcevector, double *resultvector){
	//We delcare our 3 needed vectors for diagonalization
	vector<double> maindiagonal (npoints,a);
    vector<double> diagonal (npoints, b);
    double coef;

    //We start forward substitution
    for(int i=0;i<npoints-1;i++){
        coef= diagonal[i]/maindiagonal[i];                          
        maindiagonal[i+1]-= diagonal[i]*coef;                       
        sourcevector[i+2]-=sourcevector[i+1]*coef;              
     }
    //We start backward sbustitution, ending the "diagonalization" process
    for(int i=npoints-1;i>0; i--){
        coef=diagonal[i]/maindiagonal[i];                           
        sourcevector[i]-=sourcevector[i+1]*coef;                    
    }
    for(int i=0;i<npoints;i++){
        resultvector[i+1]=sourcevector[i+1]/maindiagonal[i];          
    } 
}

void output(string name, int npoints, double step, double tstep, int timesteps, double *v){
	string outputname;
    stringstream str2, str1, str3 ;
    str1 << npoints;
    str2 << tstep;
    str3 << timesteps*tstep;
    outputname = str1.str();
    outputname += "points";
    outputname += str2.str();
    outputname += "tstep";
	outputname += str3.str();
    outputname += "time";
    outputname += name;
    outputname += ".txt";

	ofstream myfile;
	myfile.open(outputname.c_str());
	for(int i=0; i<=npoints+1;i++){
		myfile << setprecision(8) << i*step<< setw(15);

		myfile << setprecision(8) << v[i] << endl;
	}
	myfile.close();
}