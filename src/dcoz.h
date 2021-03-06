// a class for filtering in the x direction
// member functions are in fx.cpp

#ifndef DCOZ_H
#define DCOZ_H

class dcoz {
public:
	dcoz( double *f, const int Nx, const int Ny, const int Nz, const int Nm, double dyy );
//	void setmn( int m, int n );
	void derivate_co( int m, int n );

private:
	double *const f;		// the variable need to Interpolated to or derivated at
//	int m;					// the last no. of array of f
//	int n;
	const int Nx;			// number of nodes - 1 in the desired direction of interpolation
	const int Ny;			// no. of total elements in the y direction of f
	const int Nz;			// no. of total elements in the z direction of f
	const int Nm;			// no. of total element in the last array of f

	int N;
	int Nm1;
	int Nm2;			// Np1 = N + 1


	double dz, dz7_9, dz_36, dz_2, dz3_4;

	double *b;	// variables before of the diagonal
	double *d;	// variables at the diagonal
	double *a;	// variables after of the diagonal
	double *c;	// variables at the right hand side

	void bda_coeff();						// calculating b, d, a once for entire of code
	void a_const();
	void bda_Thomas();					// evaluating dp for derivation and Interpolation
};

#endif