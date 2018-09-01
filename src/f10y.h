// a class for filtering in the x direction
// member functions are in fx.cpp

#ifndef F10Y_H
#define F10Y_H

class f10y {
public:
	f10y( double *f, const int Nx, const int Ny, const int Nz, const int Nm );
	void setmn( int m, int n );
	void setalpha( double alff );
	void filter();

private:
	double *const f;		// the variable need to Interpolated to or derivated at
	int m;					// the last no. of array of f
	int n;
	const int Nx;			// number of nodes - 1 in the desired direction of interpolation
	const int Ny;			// no. of total elements in the y direction of f
	const int Nz;
	const int Nm;			// no. of total element in the last array of f

	int N;
	int Nm1;
	int Nm2;			// Np1 = N + 1
	int Nm3;			// Np2 = N + 2
	int Nm4;			// Nm1 = N - 1
	int Nm5;			// Nm2 = N - 2
	int Nm6;
	int Nm7;
	int Nm8;
	int Nm9;
	int Nm10;	

	double alf;
	double a[ 12 ][ 6 ], a0, a1, a2, a3, a4, a5;

	double *bf;	// variables before of the diagonal
	double *df;	// variables at the diagonal
	double *af;	// variables after of the diagonal
	double *cf;	// variables at the right hand side

	void bda_coeff();						// calculating b, d, a once for entire of code
	void a_const();
	void bda_Thomas();					// evaluating dp for derivation and Interpolation
};

#endif