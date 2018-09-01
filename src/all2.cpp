#include "dcox.h"
#include <iostream>
using namespace std;
dcox::dcox( double *ff, const int Nxf, const int Nyf, const int Nzf, const int Nmf, double dxx )
: f( ff ), Nx( Nxf ), Ny( Nyf ), Nz( Nzf ), Nm( Nmf )
{
	N = Nx - 1;
	Nm1 = N - 1;
	Nm2 = N - 2;

	dx = dxx;
	a_const();

	bda_coeff();
	bda_Thomas();
}
/*
void dcox::setmn( int mm, int nn )
{
	m = mm;
	n = nn;
}
*/
void dcox::a_const()
{
	dx7_9 = 7.0 / ( dx * 9 );
	dx_36 = 1.0 / ( dx * 36 );
	dx_2 = 1.0 / ( dx * 2 );
	dx3_4 = 3.0 / ( dx * 4 );
}

void dcox::bda_coeff()
{
	int i;

	b = new double[ Nx ];	// variables before of the diagonal
	d = new double[ Nx ];	// variables at the diagonal
	a = new double[ Nx ];	// variables after of the diagonal
	c = new double[ Nx ];	// variables at the right hand side

	for ( i = 3; i < Nm1; i++ )
	{
		b[ i ] = 1.0 / 3.0;
		d[ i ] = 1;
		a[ i ] = 1.0 / 3.0;
	}

	d[ 1 ] = 1;
	a[ 1 ] = 2;

	b[ 2 ] = 0.25;
	d[ 2 ] = 1;
	a[ 2 ] = 0.25;

	d[ N ] = 1;
	b[ N ] = 2;

	b[ Nm1 ] = 0.25;
	d[ Nm1 ] = 1;
	a[ Nm1 ] = 0.25;
}

void dcox::bda_Thomas()
{
	int im1, i;

	for ( i = 2; i <= N; i++ )
	{
		im1 = i - 1;
		d[ i ] -= b[ i ] * a[ im1 ] / d[ im1 ];
	}
}

void dcox::derivate_co( int m, int n )
{
	int i, im1, j, k, cst1, cst2, cst3;

	cst1 = Ny * Nz * Nm;

	for ( j = 1; j < Ny; j++ )
	{
		cst3 = j * Nz * Nm;

		for ( k = 1; k < Nz; k++ )
		{
			cst2 = cst3 + k * Nm + n;

			for ( i = 3; i < Nm1; i++ )
			{
				c[ i ] = dx7_9 * ( f[ ( i + 1 ) * cst1 + cst2 ] - f[ ( i - 1 ) * cst1 + cst2 ] ) +
							dx_36 * ( f[ ( i + 2 ) * cst1 + cst2 ] - f[ ( i - 2 ) * cst1 + cst2 ] ) ;
			}

			c[ 1 ] = dx_2 * ( -5.0 * f[ cst1 + cst2 ] + 4.0 * f[ 2 * cst1 + cst2 ] + f[ 3 * cst1 + cst2 ] );
			c[ N ] = dx_2 * ( 5.0 * f[ N * cst1 + cst2 ] - 4.0 * f[ Nm1 * cst1 + cst2 ] - f[ Nm2 * cst1 + cst2 ] );

			c[ 2 ] = dx3_4 * ( f[ 3 * cst1 + cst2 ] - f[ cst1 + cst2 ] );
			c[ Nm1 ] = dx3_4 * ( f[ N * cst1 + cst2 ] - f[ Nm2 * cst1 + cst2 ] );

			for ( i = 2; i <= N; i++ )
			{
				im1 = i - 1;
				c[ i ] -= c[ im1 ] * b[ i ] / d[ im1 ];
			}

			cst2 = cst3 + k * Nm + m;
		
			f[ N * cst1 + cst2 ] = c[ N ] / d[ N ];

			for ( i = Nm1; i >= 1; i-- )
			{
				f[ i * cst1 + cst2 ] = ( c[ i ] - a[ i ] * f[ ( i + 1 ) * cst1 + cst2 ] ) / d[ i ];
			}
		}
	}
}

#include "dcoy.h"
#include <iostream>

dcoy::dcoy( double *ff, const int Nxf, const int Nyf, const int Nzf, const int Nmf, double dyy )
: f( ff ), Nx( Nxf ), Ny( Nyf ), Nz( Nzf ), Nm( Nmf )
{
	N = Ny - 1;
	Nm1 = N - 1;
	Nm2 = N - 2;

	dy = dyy;
	a_const();

	bda_coeff();
	bda_Thomas();
}
/*
void dcoy::setmn( int mm, int nn )
{
	m = mm;
	n = nn;
}
*/
void dcoy::a_const()
{
	dy7_9 = 7.0 / ( dy * 9 );
	dy_36 = 1.0 / ( dy * 36 );
	dy_2 = 1.0 / ( dy * 2 );
	dy3_4 = 3.0 / ( dy * 4 );
}

void dcoy::bda_coeff()
{
	int i;

	b = new double[ Ny ];	// variables before of the diagonal
	d = new double[ Ny ];	// variables at the diagonal
	a = new double[ Ny ];	// variables after of the diagonal
	c = new double[ Ny ];	// variables at the right hand side

	for ( i = 3; i < Nm1; i++ )
	{
		b[ i ] = 1.0 / 3.0;
		d[ i ] = 1;
		a[ i ] = 1.0 / 3.0;
	}

	d[ 1 ] = 1;
	a[ 1 ] = 2;

	b[ 2 ] = 0.25;
	d[ 2 ] = 1;
	a[ 2 ] = 0.25;

	d[ N ] = 1;
	b[ N ] = 2;

	b[ Nm1 ] = 0.25;
	d[ Nm1 ] = 1;
	a[ Nm1 ] = 0.25;
}

void dcoy::bda_Thomas()
{
	int im1, i;

	for ( i = 2; i <= N; i++ )
	{
		im1 = i - 1;
		d[ i ] -= b[ i ] * a[ im1 ] / d[ im1 ];
	}
}

void dcoy::derivate_co( int m, int n )
{
	int j, jm1, i, k, cst1, cst2, cst3;

	cst1 = Nz * Nm;

	for ( i = 1; i < Nx; i++ )
	{
		cst3 = i * Ny * Nz * Nm;

		for ( k = 1; k < Nz; k++ )
		{
			cst2 = cst3 + k * Nm + n;

			for ( j = 3; j < Nm1; j++ )
			{
				c[ j ] = dy7_9 * ( f[ ( j + 1 ) * cst1 + cst2 ] - f[ ( j - 1 ) * cst1 + cst2 ] ) +
							dy_36 * ( f[ ( j + 2 ) * cst1 + cst2 ] - f[ ( j - 2 ) * cst1 + cst2 ] ) ;
			}

			c[ 1 ] = dy_2 * ( -5.0 * f[ cst1 + cst2 ] + 4.0 * f[ 2 * cst1 + cst2 ] + f[ 3 * cst1 + cst2 ] );
			c[ N ] = dy_2 * ( 5.0 * f[ N * cst1 + cst2 ] - 4.0 * f[ Nm1 * cst1 + cst2 ] - f[ Nm2 * cst1 + cst2 ] );

			c[ 2 ] = dy3_4 * ( f[ 3 * cst1 + cst2 ] - f[ cst1 + cst2 ] );
			c[ Nm1 ] = dy3_4 * ( f[ N * cst1 + cst2 ] - f[ Nm2 * cst1 + cst2 ] );

			for ( j = 2; j <= N; j++ )
			{
				jm1 = j - 1;
				c[ j ] -= c[ jm1 ] * b[ j ] / d[ jm1 ];
			}

			cst2 = cst3 + k * Nm + m;
		
			f[ N * cst1 + cst2 ] = c[ N ] / d[ N ];

			for ( j = Nm1; j >= 1; j-- )
			{
				f[ j * cst1 + cst2 ] = ( c[ j ] - a[ j ] * f[ ( j + 1 ) * cst1 + cst2 ] ) / d[ j ];
			}
		}
	}
}

#include "dcoz.h"
#include <iostream>

dcoz::dcoz( double *ff, const int Nxf, const int Nyf, const int Nzf, const int Nmf, double dzz )
: f( ff ), Nx( Nxf ), Ny( Nyf ), Nz( Nzf ), Nm( Nmf )
{
	N = Nz - 1;
	Nm1 = N - 1;
	Nm2 = N - 2;

	dz = dzz;
	a_const();

	bda_coeff();
	bda_Thomas();
}
/*
void dcoz::setmn( int mm, int nn )
{
	m = mm;
	n = nn;
}
*/
void dcoz::a_const()
{
	dz7_9 = 7.0 / ( dz * 9 );
	dz_36 = 1.0 / ( dz * 36 );
	dz_2 = 1.0 / ( dz * 2 );
	dz3_4 = 3.0 / ( dz * 4 );
}

void dcoz::bda_coeff()
{
	int i;

	b = new double[ Nz ];	// variables before of the diagonal
	d = new double[ Nz ];	// variables at the diagonal
	a = new double[ Nz ];	// variables after of the diagonal
	c = new double[ Nz ];	// variables at the right hand side

	for ( i = 3; i < Nm1; i++ )
	{
		b[ i ] = 1.0 / 3.0;
		d[ i ] = 1;
		a[ i ] = 1.0 / 3.0;
	}

	d[ 1 ] = 1;
	a[ 1 ] = 2;

	b[ 2 ] = 0.25;
	d[ 2 ] = 1;
	a[ 2 ] = 0.25;

	d[ N ] = 1;
	b[ N ] = 2;

	b[ Nm1 ] = 0.25;
	d[ Nm1 ] = 1;
	a[ Nm1 ] = 0.25;
}

void dcoz::bda_Thomas()
{
	int im1, i;

	for ( i = 2; i <= N; i++ )
	{
		im1 = i - 1;
		d[ i ] -= b[ i ] * a[ im1 ] / d[ im1 ];
	}
}

void dcoz::derivate_co( int m, int n )
{
	int j, km1, i, k, cst1, cst2, cst3;

	cst1 = Nm;

	for ( i = 1; i < Nx; i++ )
	{
		cst3 = i * Ny * Nz * Nm;

		for ( j = 1; j < Ny; j++ )
		{
			cst2 = cst3 + j * Nz * Nm + n;

			for ( k = 3; k < Nm1; k++ )
			{
				c[ k ] = dz7_9 * ( f[ ( k + 1 ) * cst1 + cst2 ] - f[ ( k - 1 ) * cst1 + cst2 ] ) +
							dz_36 * ( f[ ( k + 2 ) * cst1 + cst2 ] - f[ ( k - 2 ) * cst1 + cst2 ] ) ;
			}

			c[ 1 ] = dz_2 * ( -5.0 * f[ cst1 + cst2 ] + 4.0 * f[ 2 * cst1 + cst2 ] + f[ 3 * cst1 + cst2 ] );
			c[ N ] = dz_2 * ( 5.0 * f[ N * cst1 + cst2 ] - 4.0 * f[ Nm1 * cst1 + cst2 ] - f[ Nm2 * cst1 + cst2 ] );

			c[ 2 ] = dz3_4 * ( f[ 3 * cst1 + cst2 ] - f[ cst1 + cst2 ] );
			c[ Nm1 ] = dz3_4 * ( f[ N * cst1 + cst2 ] - f[ Nm2 * cst1 + cst2 ] );

			for ( k = 2; k <= N; k++ )
			{
				km1 = k - 1;
				c[ k ] -= c[ km1 ] * b[ k ] / d[ km1 ];
			}

			cst2 = cst3 + j * Nz * Nm + m;
		
			f[ N * cst1 + cst2 ] = c[ N ] / d[ N ];

			for ( k = Nm1; k >= 1; k-- )
			{
				f[ k * cst1 + cst2 ] = ( c[ k ] - a[ k ] * f[ ( k + 1 ) * cst1 + cst2 ] ) / d[ k ];
			}
		}
	}
}
