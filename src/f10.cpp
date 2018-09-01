#include "f10x.h"
#include <iostream>
using namespace std;

f10x::f10x( double *ff, const int Nxf, const int Nyf, const int Nzf, const int Nmf )
: f( ff ), Nx( Nxf ), Ny( Nyf ), Nz( Nzf ), Nm( Nmf )
{	
	N = Nx - 1;
	Nm1 = N - 1;
	Nm2 = N - 2;
	Nm3 = N - 3;
	Nm4 = N - 4;
	Nm5 = N - 5;
	Nm6 = N - 6;
	Nm7 = N - 7;
	Nm8 = N - 8;
	Nm9 = N - 9;
	Nm10 = N - 10;
}

void f10x::setmn( int mm, int nn )
{	
	m = mm;
	n = nn;
}

void f10x::setalpha( double alff )
{
	alf = alff;

	a_const();
	bda_coeff();
	bda_Thomas();
}

void f10x::a_const()
{
	a0 = ( 193.0 + 126.0 * alf ) / 256;
	a1 = ( 105.0 + 302.0 * alf ) / 512;
	a2 = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a3 = 45.0 * ( 1.0 - 2.0 * alf ) / 1024;
	a4 = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a5 = ( 1 - 2.0 * alf ) / 1024;

	a[1][2] = ( 1.0 + 1022.0 * alf ) / 1024;
	a[2][2] = ( 507.0 + 10.0 * alf ) / 512;
	a[3][2] = ( 45.0 + 934.0 * alf ) / 1024;
	a[4][2] = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a[5][2] = 105.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[6][2] = 63.0 * ( -1.0 + 2.0 * alf ) / 256;
	a[7][2] = 105.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[8][2] = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a[9][2] = 45.0 * ( 1.0 - 2.0 * alf ) / 1024;
	a[10][2] = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[11][2] = ( 1.0 - 2.0 * alf ) / 1024;

	a[1][3] = ( -1.0 + 2.0 * alf ) / 1024;
	a[2][3] = ( 5.0 + 502.0 * alf ) / 512;
	a[3][3] = ( 979.0 + 90.0 * alf ) / 1024;
	a[4][3] = ( 15.0 + 98.0 * alf ) / 128;
	a[5][3] = 105.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[6][3] = 63.0 * ( 1.0 - 2.0 * alf ) / 256;
	a[7][3] = 105.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[8][3] = 15.0 * ( 1.0 - 2.0 * alf ) / 128;
	a[9][3] = 45.0 * ( -1.0 + 2.0 * alf ) / 1024;
	a[10][3] = 5.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[11][3] = ( -1.0 + 2.0 * alf ) / 1024;

	a[1][4] = ( 1.0 - 2.0 * alf ) / 1024;
	a[2][4] = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[3][4] = ( 45.0 + 934.0 * alf ) / 1024;
	a[4][4] = ( 113.0 + 30.0 * alf ) / 128;
	a[5][4] = ( 105.0 + 302.0 * alf ) / 512;
	a[6][4] = 63.0 * ( -1.0 + 2.0 * alf ) / 256;
	a[7][4] = 105.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[8][4] = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a[9][4] = 45.0 * ( 1.0 - 2.0 * alf ) / 1024;
	a[10][4] = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[11][4] = ( 1.0 - 2.0 * alf ) / 1024;

	a[1][5] = ( -1.0 + 2.0 * alf ) / 1024;
	a[2][5] = 5.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[3][5] = 45.0 * ( -1.0 + 2 * alf ) / 1024;
	a[4][5] = ( 15.0 + 98.0 * alf ) / 128;
	a[5][5] = ( 407.0 + 210.0 * alf ) / 512;
	a[6][5] = ( 63.0 + 130.0 * alf ) / 256;
	a[7][5] = 105.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[8][5] = 15.0 * ( 1.0 - 2.0 * alf ) / 128;
	a[9][5] = 45.0 * ( -1.0 + 2.0 * alf ) / 1024;
	a[10][5] = 5.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[11][5] = ( -1.0 + 2.0 * alf ) / 1024;
}

void f10x::bda_coeff()
{	
	int i;

	bf = new double[ Nx ];	// variables before of the diagonal
	df = new double[ Nx ];	// variables at the diagonal
	af = new double[ Nx ];	// variables after of the diagonal
	cf = new double[ Nx ];	// variables at the right hand side

	for ( i = 2; i < N; i++ )
	{
		bf[ i ] = alf;
		df[ i ] = 1;
		af[ i ] = alf;
	}
}

void f10x::bda_Thomas()
{
	int im1, i;

	for ( i = 3; i < N; i++ )
	{
		im1 = i - 1;
		df[ i ] -= bf[ i ] * af[ im1 ] / df[ im1 ];
	}
}

void f10x::filter()
{
	int i, im1, j, k, cst1, cst2m, cst2n, cst3;

	cst1 = Ny * Nz * Nm;

	for ( j = 2; j < Ny - 1; j++ )	
	{
		cst3 = j * Nz * Nm;

	for ( k = 2; k < Nz - 1; k++ )	
	{
		cst2m = cst3 + k * Nm + m;

		for ( i = 6; i <= Nm5; i++ )
		{
			cf[ i ] =  a0 * f[ i * cst1 + cst2m ] + 
						  a1 * ( f[ ( i + 1 ) * cst1 + cst2m ] + f[ ( i - 1 ) * cst1 + cst2m ] ) +
						  a2 * ( f[ ( i + 2 ) * cst1 + cst2m ] + f[ ( i - 2 ) * cst1 + cst2m ] ) +
						  a3 * ( f[ ( i + 3 ) * cst1 + cst2m ] + f[ ( i - 3 ) * cst1 + cst2m ] ) +
						  a4 * ( f[ ( i + 4 ) * cst1 + cst2m ] + f[ ( i - 4 ) * cst1 + cst2m ] ) +
						  a5 * ( f[ ( i + 5 ) * cst1 + cst2m ] + f[ ( i - 5 ) * cst1 + cst2m ] ) ;
		}

		cf[ 2 ] = a[1][2] * f[ cst1 + cst2m ] + a[2][2] * f[ 2 * cst1 + cst2m ] + a[3][2] * f[ 3 * cst1 + cst2m ] +
					 a[4][2] * f[ 4 * cst1 + cst2m ] + a[5][2] * f[ 5 * cst1 + cst2m ] + a[6][2] * f[ 6 * cst1 + cst2m ] +
					 a[7][2] * f[ 7 * cst1 + cst2m ] + a[8][2] * f[ 8 * cst1 + cst2m ] + a[9][2] * f[ 9 * cst1 + cst2m ] +
					 a[10][2] * f[ 10 * cst1 + cst2m ] + a[11][2] * f[ 11 * cst1 + cst2m ] - alf * f[ cst1 + cst2m ];

		cf[ 3 ] = a[1][3] * f[ cst1 + cst2m ] + a[2][3] * f[ 2 * cst1 + cst2m ] + a[3][3] * f[ 3 * cst1 + cst2m ] +
					 a[4][3] * f[ 4 * cst1 + cst2m ] + a[5][3] * f[ 5 * cst1 + cst2m ] + a[6][3] * f[ 6 * cst1 + cst2m ] +
					 a[7][3] * f[ 7 * cst1 + cst2m ] + a[8][3] * f[ 8 * cst1 + cst2m ] + a[9][3] * f[ 9 * cst1 + cst2m ] +
					 a[10][3] * f[ 10 * cst1 + cst2m ] + a[11][3] * f[ 11 * cst1 + cst2m ];

		cf[ 4 ] = a[1][4] * f[ cst1 + cst2m ] + a[2][4] * f[ 2 * cst1 + cst2m ] + a[3][4] * f[ 3 * cst1 + cst2m ] +
					 a[4][4] * f[ 4 * cst1 + cst2m ] + a[5][4] * f[ 5 * cst1 + cst2m ] + a[6][4] * f[ 6 * cst1 + cst2m ] +
					 a[7][4] * f[ 7 * cst1 + cst2m ] + a[8][4] * f[ 8 * cst1 + cst2m ] + a[9][4] * f[ 9 * cst1 + cst2m ] +
					 a[10][4] * f[ 10 * cst1 + cst2m ] + a[11][4] * f[ 11 * cst1 + cst2m ];

		cf[ 5 ] = a[1][5] * f[ cst1 + cst2m ] + a[2][5] * f[ 2 * cst1 + cst2m ] + a[3][5] * f[ 3 * cst1 + cst2m ] +
					 a[4][5] * f[ 4 * cst1 + cst2m ] + a[5][5] * f[ 5 * cst1 + cst2m ] + a[6][5] * f[ 6 * cst1 + cst2m ] +
					 a[7][5] * f[ 7 * cst1 + cst2m ] + a[8][5] * f[ 8 * cst1 + cst2m ] + a[9][5] * f[ 9 * cst1 + cst2m ] +
					 a[10][5] * f[ 10 * cst1 + cst2m ] + a[11][5] * f[ 11 * cst1 + cst2m ];

		cf[ Nm1 ] = a[1][2] * f[ N * cst1 + cst2m ] + a[2][2] * f[ Nm1 * cst1 + cst2m ] + a[3][2] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][2] * f[ Nm3 * cst1 + cst2m ] + a[5][2] * f[ Nm4 * cst1 + cst2m ] + a[6][2] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][2] * f[ Nm6 * cst1 + cst2m ] + a[8][2] * f[ Nm7 * cst1 + cst2m ] + a[9][2] * f[ Nm8 * cst1 + cst2m ] +
						a[10][2] * f[ Nm9 * cst1 + cst2m ] + a[11][2] * f[ Nm10 * cst1 + cst2m ] - alf * f[ N * cst1 + cst2m ];

		cf[ Nm2 ] = a[1][3] * f[ N * cst1 + cst2m ] + a[2][3] * f[ Nm1 * cst1 + cst2m ] + a[3][3] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][3] * f[ Nm3 * cst1 + cst2m ] + a[5][3] * f[ Nm4 * cst1 + cst2m ] + a[6][3] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][3] * f[ Nm6 * cst1 + cst2m ] + a[8][3] * f[ Nm7 * cst1 + cst2m ] + a[9][3] * f[ Nm8 * cst1 + cst2m ] +
						a[10][3] * f[ Nm9 * cst1 + cst2m ] + a[11][3] * f[ Nm10 * cst1 + cst2m ];

		cf[ Nm3 ] = a[1][4] * f[ N * cst1 + cst2m ] + a[2][4] * f[ Nm1 * cst1 + cst2m ] + a[3][4] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][4] * f[ Nm3 * cst1 + cst2m ] + a[5][4] * f[ Nm4 * cst1 + cst2m ] + a[6][4] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][4] * f[ Nm6 * cst1 + cst2m ] + a[8][4] * f[ Nm7 * cst1 + cst2m ] + a[9][4] * f[ Nm8 * cst1 + cst2m ] +
						a[10][4] * f[ Nm9 * cst1 + cst2m ] + a[11][4] * f[ Nm10 * cst1 + cst2m ];

		cf[ Nm4 ] = a[1][5] * f[ N * cst1 + cst2m ] + a[2][5] * f[ Nm1 * cst1 + cst2m ] + a[3][5] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][5] * f[ Nm3 * cst1 + cst2m ] + a[5][5] * f[ Nm4 * cst1 + cst2m ] + a[6][5] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][5] * f[ Nm6 * cst1 + cst2m ] + a[8][5] * f[ Nm7 * cst1 + cst2m ] + a[9][5] * f[ Nm8 * cst1 + cst2m ] +
						a[10][5] * f[ Nm9 * cst1 + cst2m ] + a[11][5] * f[ Nm10 * cst1 + cst2m ];


		for ( i = 3; i <= Nm1; i++ )
		{
			im1 = i - 1;
			cf[ i ] -= cf[ im1 ] * bf[ i ] / df[ im1 ];
		}

		cst2n = cst3 + k * Nm + n;
		
		f[ Nm1 * cst1 + cst2n ] = cf[ Nm1 ] / df[ Nm1 ];

		for ( i = Nm2; i >= 2; i-- )
		{
			f[ i * cst1 + cst2n ] = ( cf[ i ] - af[ i ] * f[ ( i + 1 ) * cst1 + cst2n ] ) / df[ i ];
		}
	}
	}
}

#include "f10y.h"
#include <iostream>

f10y::f10y( double *ff, const int Nxf, const int Nyf, const int Nzf, const int Nmf )
: f( ff ), Nx( Nxf ), Ny( Nyf ), Nz( Nzf ), Nm( Nmf )
{
	N = Ny - 1;
	Nm1 = N - 1;
	Nm2 = N - 2;
	Nm3 = N - 3;
	Nm4 = N - 4;
	Nm5 = N - 5;
	Nm6 = N - 6;
	Nm7 = N - 7;
	Nm8 = N - 8;
	Nm9 = N - 9;
	Nm10 = N - 10;
}

void f10y::setmn( int mm, int nn )
{	
	m = mm;
	n = nn;
}

void f10y::setalpha( double alff )
{
	alf = alff;

	a_const();
	bda_coeff();
	bda_Thomas();
}

void f10y::a_const()
{
	a0 = ( 193.0 + 126.0 * alf ) / 256;
	a1 = ( 105.0 + 302.0 * alf ) / 512;
	a2 = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a3 = 45.0 * ( 1.0 - 2.0 * alf ) / 1024;
	a4 = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a5 = ( 1 - 2.0 * alf ) / 1024;

	a[1][2] = ( 1.0 + 1022.0 * alf ) / 1024;
	a[2][2] = ( 507.0 + 10.0 * alf ) / 512;
	a[3][2] = ( 45.0 + 934.0 * alf ) / 1024;
	a[4][2] = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a[5][2] = 105.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[6][2] = 63.0 * ( -1.0 + 2.0 * alf ) / 256;
	a[7][2] = 105.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[8][2] = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a[9][2] = 45.0 * ( 1.0 - 2.0 * alf ) / 1024;
	a[10][2] = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[11][2] = ( 1.0 - 2.0 * alf ) / 1024;

	a[1][3] = ( -1.0 + 2.0 * alf ) / 1024;
	a[2][3] = ( 5.0 + 502.0 * alf ) / 512;
	a[3][3] = ( 979.0 + 90.0 * alf ) / 1024;
	a[4][3] = ( 15.0 + 98.0 * alf ) / 128;
	a[5][3] = 105.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[6][3] = 63.0 * ( 1.0 - 2.0 * alf ) / 256;
	a[7][3] = 105.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[8][3] = 15.0 * ( 1.0 - 2.0 * alf ) / 128;
	a[9][3] = 45.0 * ( -1.0 + 2.0 * alf ) / 1024;
	a[10][3] = 5.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[11][3] = ( -1.0 + 2.0 * alf ) / 1024;

	a[1][4] = ( 1.0 - 2.0 * alf ) / 1024;
	a[2][4] = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[3][4] = ( 45.0 + 934.0 * alf ) / 1024;
	a[4][4] = ( 113.0 + 30.0 * alf ) / 128;
	a[5][4] = ( 105.0 + 302.0 * alf ) / 512;
	a[6][4] = 63.0 * ( -1.0 + 2.0 * alf ) / 256;
	a[7][4] = 105.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[8][4] = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a[9][4] = 45.0 * ( 1.0 - 2.0 * alf ) / 1024;
	a[10][4] = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[11][4] = ( 1.0 - 2.0 * alf ) / 1024;

	a[1][5] = ( -1.0 + 2.0 * alf ) / 1024;
	a[2][5] = 5.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[3][5] = 45.0 * ( -1.0 + 2 * alf ) / 1024;
	a[4][5] = ( 15.0 + 98.0 * alf ) / 128;
	a[5][5] = ( 407.0 + 210.0 * alf ) / 512;
	a[6][5] = ( 63.0 + 130.0 * alf ) / 256;
	a[7][5] = 105.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[8][5] = 15.0 * ( 1.0 - 2.0 * alf ) / 128;
	a[9][5] = 45.0 * ( -1.0 + 2.0 * alf ) / 1024;
	a[10][5] = 5.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[11][5] = ( -1.0 + 2.0 * alf ) / 1024;
}

void f10y::bda_coeff()
{
	int i;

	bf = new double[ Ny ];	// variables before of the diagonal
	df = new double[ Ny ];	// variables at the diagonal
	af = new double[ Ny ];	// variables after of the diagonal
	cf = new double[ Ny ];	// variables at the right hand side

	for ( i = 2; i < N; i++ )
	{
		bf[ i ] = alf;
		df[ i ] = 1;
		af[ i ] = alf;
	}
}

void f10y::bda_Thomas()
{
	int im1, i;

	for ( i = 3; i < N; i++ )
	{
		im1 = i - 1;
		df[ i ] -= bf[ i ] * af[ im1 ] / df[ im1 ];
	}
}

void f10y::filter()
{
	int j, jm1, i, k, cst1, cst2m, cst2n, cst3;

	cst1 = Nz * Nm;

	for ( k = 2; k < Nz - 1; k++ )
	{
		cst3 = k * Nm;

	for ( i = 2; i < Nx - 1; i++ )
	{
		cst2m = i * Ny * Nz * Nm + cst3 + m;

		for ( j = 4; j <= Nm3; j++ )
		{
			cf[ j ] =  a0 * f[ j * cst1 + cst2m ] + 
						  a1 * ( f[ ( j + 1 ) * cst1 + cst2m ] + f[ ( j - 1 ) * cst1 + cst2m ] ) +
						  a2 * ( f[ ( j + 2 ) * cst1 + cst2m ] + f[ ( j - 2 ) * cst1 + cst2m ] ) +
						  a3 * ( f[ ( j + 3 ) * cst1 + cst2m ] + f[ ( j - 3 ) * cst1 + cst2m ] ) +
						  a4 * ( f[ ( j + 4 ) * cst1 + cst2m ] + f[ ( j - 4 ) * cst1 + cst2m ] ) +
						  a5 * ( f[ ( j + 5 ) * cst1 + cst2m ] + f[ ( j - 5 ) * cst1 + cst2m ] ) ;
		}

		cf[ 2 ] = a[1][2] * f[ cst1 + cst2m ] + a[2][2] * f[ 2 * cst1 + cst2m ] + a[3][2] * f[ 3 * cst1 + cst2m ] +
					 a[4][2] * f[ 4 * cst1 + cst2m ] + a[5][2] * f[ 5 * cst1 + cst2m ] + a[6][2] * f[ 6 * cst1 + cst2m ] +
					 a[7][2] * f[ 7 * cst1 + cst2m ] + a[8][2] * f[ 8 * cst1 + cst2m ] + a[9][2] * f[ 9 * cst1 + cst2m ] +
					 a[10][2] * f[ 10 * cst1 + cst2m ] + a[11][2] * f[ 11 * cst1 + cst2m ] - alf * f[ cst1 + cst2m ];

		cf[ 3 ] = a[1][3] * f[ cst1 + cst2m ] + a[2][3] * f[ 2 * cst1 + cst2m ] + a[3][3] * f[ 3 * cst1 + cst2m ] +
					 a[4][3] * f[ 4 * cst1 + cst2m ] + a[5][3] * f[ 5 * cst1 + cst2m ] + a[6][3] * f[ 6 * cst1 + cst2m ] +
					 a[7][3] * f[ 7 * cst1 + cst2m ] + a[8][3] * f[ 8 * cst1 + cst2m ] + a[9][3] * f[ 9 * cst1 + cst2m ] +
					 a[10][3] * f[ 10 * cst1 + cst2m ] + a[11][3] * f[ 11 * cst1 + cst2m ];

		cf[ 4 ] = a[1][4] * f[ cst1 + cst2m ] + a[2][4] * f[ 2 * cst1 + cst2m ] + a[3][4] * f[ 3 * cst1 + cst2m ] +
					 a[4][4] * f[ 4 * cst1 + cst2m ] + a[5][4] * f[ 5 * cst1 + cst2m ] + a[6][4] * f[ 6 * cst1 + cst2m ] +
					 a[7][4] * f[ 7 * cst1 + cst2m ] + a[8][4] * f[ 8 * cst1 + cst2m ] + a[9][4] * f[ 9 * cst1 + cst2m ] +
					 a[10][4] * f[ 10 * cst1 + cst2m ] + a[11][4] * f[ 11 * cst1 + cst2m ];

		cf[ 5 ] = a[1][5] * f[ cst1 + cst2m ] + a[2][5] * f[ 2 * cst1 + cst2m ] + a[3][5] * f[ 3 * cst1 + cst2m ] +
					 a[4][5] * f[ 4 * cst1 + cst2m ] + a[5][5] * f[ 5 * cst1 + cst2m ] + a[6][5] * f[ 6 * cst1 + cst2m ] +
					 a[7][5] * f[ 7 * cst1 + cst2m ] + a[8][5] * f[ 8 * cst1 + cst2m ] + a[9][5] * f[ 9 * cst1 + cst2m ] +
					 a[10][5] * f[ 10 * cst1 + cst2m ] + a[11][5] * f[ 11 * cst1 + cst2m ];

		cf[ Nm1 ] = a[1][2] * f[ N * cst1 + cst2m ] + a[2][2] * f[ Nm1 * cst1 + cst2m ] + a[3][2] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][2] * f[ Nm3 * cst1 + cst2m ] + a[5][2] * f[ Nm4 * cst1 + cst2m ] + a[6][2] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][2] * f[ Nm6 * cst1 + cst2m ] + a[8][2] * f[ Nm7 * cst1 + cst2m ] + a[9][2] * f[ Nm8 * cst1 + cst2m ] +
						a[10][2] * f[ Nm9 * cst1 + cst2m ] + a[11][2] * f[ Nm10 * cst1 + cst2m ] - alf * f[ N * cst1 + cst2m ];

		cf[ Nm2 ] = a[1][3] * f[ N * cst1 + cst2m ] + a[2][3] * f[ Nm1 * cst1 + cst2m ] + a[3][3] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][3] * f[ Nm3 * cst1 + cst2m ] + a[5][3] * f[ Nm4 * cst1 + cst2m ] + a[6][3] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][3] * f[ Nm6 * cst1 + cst2m ] + a[8][3] * f[ Nm7 * cst1 + cst2m ] + a[9][3] * f[ Nm8 * cst1 + cst2m ] +
						a[10][3] * f[ Nm9 * cst1 + cst2m ] + a[11][3] * f[ Nm10 * cst1 + cst2m ];

		cf[ Nm3 ] = a[1][4] * f[ N * cst1 + cst2m ] + a[2][4] * f[ Nm1 * cst1 + cst2m ] + a[3][4] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][4] * f[ Nm3 * cst1 + cst2m ] + a[5][4] * f[ Nm4 * cst1 + cst2m ] + a[6][4] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][4] * f[ Nm6 * cst1 + cst2m ] + a[8][4] * f[ Nm7 * cst1 + cst2m ] + a[9][4] * f[ Nm8 * cst1 + cst2m ] +
						a[10][4] * f[ Nm9 * cst1 + cst2m ] + a[11][4] * f[ Nm10 * cst1 + cst2m ];

		cf[ Nm4 ] = a[1][5] * f[ N * cst1 + cst2m ] + a[2][5] * f[ Nm1 * cst1 + cst2m ] + a[3][5] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][5] * f[ Nm3 * cst1 + cst2m ] + a[5][5] * f[ Nm4 * cst1 + cst2m ] + a[6][5] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][5] * f[ Nm6 * cst1 + cst2m ] + a[8][5] * f[ Nm7 * cst1 + cst2m ] + a[9][5] * f[ Nm8 * cst1 + cst2m ] +
						a[10][5] * f[ Nm9 * cst1 + cst2m ] + a[11][5] * f[ Nm10 * cst1 + cst2m ];


		for ( j = 3; j <= Nm1; j++ )
		{
			jm1 = j - 1;
			cf[ j ] -= cf[ jm1 ] * bf[ j ] / df[ jm1 ];
		}

		cst2n = i * Ny * Nz * Nm + cst3 + n;
		
		f[ Nm1 * cst1 + cst2n ] = cf[ Nm1 ] / df[ Nm1 ];

		for ( j = Nm2; j >= 2; j-- )
		{
			f[ j * cst1 + cst2n ] = ( cf[ j ] - af[ j ] * f[ ( j + 1 ) * cst1 + cst2n ] ) / df[ j ];
		}
	}
	}
}

#include "f10z.h"
#include <iostream>

f10z::f10z( double *ff, const int Nxf, const int Nyf, const int Nzf, const int Nmf )
: f( ff ), Nx( Nxf ), Ny( Nyf ), Nz( Nzf ), Nm( Nmf )
{
	N = Nz - 1;
	Nm1 = N - 1;
	Nm2 = N - 2;
	Nm3 = N - 3;
	Nm4 = N - 4;
	Nm5 = N - 5;
	Nm6 = N - 6;
	Nm7 = N - 7;
	Nm8 = N - 8;
	Nm9 = N - 9;
	Nm10 = N - 10;
}

void f10z::setmn( int mm, int nn )
{
	m = mm;
	n = nn;
}

void f10z::setalpha( double alff )
{
	alf = alff;

	a_const();
	bda_coeff();
	bda_Thomas();
}

void f10z::a_const()
{
	a0 = ( 193.0 + 126.0 * alf ) / 256;
	a1 = ( 105.0 + 302.0 * alf ) / 512;
	a2 = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a3 = 45.0 * ( 1.0 - 2.0 * alf ) / 1024;
	a4 = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a5 = ( 1 - 2.0 * alf ) / 1024;

	a[1][2] = ( 1.0 + 1022.0 * alf ) / 1024;
	a[2][2] = ( 507.0 + 10.0 * alf ) / 512;
	a[3][2] = ( 45.0 + 934.0 * alf ) / 1024;
	a[4][2] = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a[5][2] = 105.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[6][2] = 63.0 * ( -1.0 + 2.0 * alf ) / 256;
	a[7][2] = 105.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[8][2] = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a[9][2] = 45.0 * ( 1.0 - 2.0 * alf ) / 1024;
	a[10][2] = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[11][2] = ( 1.0 - 2.0 * alf ) / 1024;

	a[1][3] = ( -1.0 + 2.0 * alf ) / 1024;
	a[2][3] = ( 5.0 + 502.0 * alf ) / 512;
	a[3][3] = ( 979.0 + 90.0 * alf ) / 1024;
	a[4][3] = ( 15.0 + 98.0 * alf ) / 128;
	a[5][3] = 105.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[6][3] = 63.0 * ( 1.0 - 2.0 * alf ) / 256;
	a[7][3] = 105.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[8][3] = 15.0 * ( 1.0 - 2.0 * alf ) / 128;
	a[9][3] = 45.0 * ( -1.0 + 2.0 * alf ) / 1024;
	a[10][3] = 5.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[11][3] = ( -1.0 + 2.0 * alf ) / 1024;

	a[1][4] = ( 1.0 - 2.0 * alf ) / 1024;
	a[2][4] = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[3][4] = ( 45.0 + 934.0 * alf ) / 1024;
	a[4][4] = ( 113.0 + 30.0 * alf ) / 128;
	a[5][4] = ( 105.0 + 302.0 * alf ) / 512;
	a[6][4] = 63.0 * ( -1.0 + 2.0 * alf ) / 256;
	a[7][4] = 105.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[8][4] = 15.0 * ( -1.0 + 2.0 * alf ) / 128;
	a[9][4] = 45.0 * ( 1.0 - 2.0 * alf ) / 1024;
	a[10][4] = 5.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[11][4] = ( 1.0 - 2.0 * alf ) / 1024;

	a[1][5] = ( -1.0 + 2.0 * alf ) / 1024;
	a[2][5] = 5.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[3][5] = 45.0 * ( -1.0 + 2 * alf ) / 1024;
	a[4][5] = ( 15.0 + 98.0 * alf ) / 128;
	a[5][5] = ( 407.0 + 210.0 * alf ) / 512;
	a[6][5] = ( 63.0 + 130.0 * alf ) / 256;
	a[7][5] = 105.0 * ( -1.0 + 2.0 * alf ) / 512;
	a[8][5] = 15.0 * ( 1.0 - 2.0 * alf ) / 128;
	a[9][5] = 45.0 * ( -1.0 + 2.0 * alf ) / 1024;
	a[10][5] = 5.0 * ( 1.0 - 2.0 * alf ) / 512;
	a[11][5] = ( -1.0 + 2.0 * alf ) / 1024;
}

void f10z::bda_coeff()
{
	int i;

	bf = new double[ Nz ];	// variables before of the diagonal
	df = new double[ Nz ];	// variables at the diagonal
	af = new double[ Nz ];	// variables after of the diagonal
	cf = new double[ Nz ];	// variables at the right hand side

	for ( i = 2; i < N; i++ )
	{
		bf[ i ] = alf;
		df[ i ] = 1;
		af[ i ] = alf;
	}
}

void f10z::bda_Thomas()
{
	int im1, i;

	for ( i = 3; i < N; i++ )
	{
		im1 = i - 1;
		df[ i ] -= bf[ i ] * af[ im1 ] / df[ im1 ];
	}
}

void f10z::filter()
{
	int i, km1, j, k, cst1, cst2m, cst2n, cst3;

	cst1 = Nm;

	for ( i = 2; i < Nx - 1; i++ )	
	{
		cst3 = i * Ny * Nz * Nm;

	for ( j = 2; j <= Ny - 1; j++ )
	{
		cst2m = cst3 + j * Nz * Nm + m;

		for ( k = 6; k < Nm4; k++ )
		{
			cf[ k ] =  a0 * f[ k * cst1 + cst2m ] + 
						  a1 * ( f[ ( k + 1 ) * cst1 + cst2m ] + f[ ( k - 1 ) * cst1 + cst2m ] ) +
						  a2 * ( f[ ( k + 2 ) * cst1 + cst2m ] + f[ ( k - 2 ) * cst1 + cst2m ] ) +
						  a3 * ( f[ ( k + 3 ) * cst1 + cst2m ] + f[ ( k - 3 ) * cst1 + cst2m ] ) +
						  a4 * ( f[ ( k + 4 ) * cst1 + cst2m ] + f[ ( k - 4 ) * cst1 + cst2m ] ) +
						  a5 * ( f[ ( k + 5 ) * cst1 + cst2m ] + f[ ( k - 5 ) * cst1 + cst2m ] ) ;
		}

		cf[ 2 ] = a[1][2] * f[ cst1 + cst2m ] + a[2][2] * f[ 2 * cst1 + cst2m ] + a[3][2] * f[ 3 * cst1 + cst2m ] +
					 a[4][2] * f[ 4 * cst1 + cst2m ] + a[5][2] * f[ 5 * cst1 + cst2m ] + a[6][2] * f[ 6 * cst1 + cst2m ] +
					 a[7][2] * f[ 7 * cst1 + cst2m ] + a[8][2] * f[ 8 * cst1 + cst2m ] + a[9][2] * f[ 9 * cst1 + cst2m ] +
					 a[10][2] * f[ 10 * cst1 + cst2m ] + a[11][2] * f[ 11 * cst1 + cst2m ] - alf * f[ cst1 + cst2m ];

		cf[ 3 ] = a[1][3] * f[ cst1 + cst2m ] + a[2][3] * f[ 2 * cst1 + cst2m ] + a[3][3] * f[ 3 * cst1 + cst2m ] +
					 a[4][3] * f[ 4 * cst1 + cst2m ] + a[5][3] * f[ 5 * cst1 + cst2m ] + a[6][3] * f[ 6 * cst1 + cst2m ] +
					 a[7][3] * f[ 7 * cst1 + cst2m ] + a[8][3] * f[ 8 * cst1 + cst2m ] + a[9][3] * f[ 9 * cst1 + cst2m ] +
					 a[10][3] * f[ 10 * cst1 + cst2m ] + a[11][3] * f[ 11 * cst1 + cst2m ];

		cf[ 4 ] = a[1][4] * f[ cst1 + cst2m ] + a[2][4] * f[ 2 * cst1 + cst2m ] + a[3][4] * f[ 3 * cst1 + cst2m ] +
					 a[4][4] * f[ 4 * cst1 + cst2m ] + a[5][4] * f[ 5 * cst1 + cst2m ] + a[6][4] * f[ 6 * cst1 + cst2m ] +
					 a[7][4] * f[ 7 * cst1 + cst2m ] + a[8][4] * f[ 8 * cst1 + cst2m ] + a[9][4] * f[ 9 * cst1 + cst2m ] +
					 a[10][4] * f[ 10 * cst1 + cst2m ] + a[11][4] * f[ 11 * cst1 + cst2m ];

		cf[ 5 ] = a[1][5] * f[ cst1 + cst2m ] + a[2][5] * f[ 2 * cst1 + cst2m ] + a[3][5] * f[ 3 * cst1 + cst2m ] +
					 a[4][5] * f[ 4 * cst1 + cst2m ] + a[5][5] * f[ 5 * cst1 + cst2m ] + a[6][5] * f[ 6 * cst1 + cst2m ] +
					 a[7][5] * f[ 7 * cst1 + cst2m ] + a[8][5] * f[ 8 * cst1 + cst2m ] + a[9][5] * f[ 9 * cst1 + cst2m ] +
					 a[10][5] * f[ 10 * cst1 + cst2m ] + a[11][5] * f[ 11 * cst1 + cst2m ];

		cf[ Nm1 ] = a[1][2] * f[ N * cst1 + cst2m ] + a[2][2] * f[ Nm1 * cst1 + cst2m ] + a[3][2] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][2] * f[ Nm3 * cst1 + cst2m ] + a[5][2] * f[ Nm4 * cst1 + cst2m ] + a[6][2] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][2] * f[ Nm6 * cst1 + cst2m ] + a[8][2] * f[ Nm7 * cst1 + cst2m ] + a[9][2] * f[ Nm8 * cst1 + cst2m ] +
						a[10][2] * f[ Nm9 * cst1 + cst2m ] + a[11][2] * f[ Nm10 * cst1 + cst2m ] - alf * f[ N * cst1 + cst2m ];

		cf[ Nm2 ] = a[1][3] * f[ N * cst1 + cst2m ] + a[2][3] * f[ Nm1 * cst1 + cst2m ] + a[3][3] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][3] * f[ Nm3 * cst1 + cst2m ] + a[5][3] * f[ Nm4 * cst1 + cst2m ] + a[6][3] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][3] * f[ Nm6 * cst1 + cst2m ] + a[8][3] * f[ Nm7 * cst1 + cst2m ] + a[9][3] * f[ Nm8 * cst1 + cst2m ] +
						a[10][3] * f[ Nm9 * cst1 + cst2m ] + a[11][3] * f[ Nm10 * cst1 + cst2m ];

		cf[ Nm3 ] = a[1][4] * f[ N * cst1 + cst2m ] + a[2][4] * f[ Nm1 * cst1 + cst2m ] + a[3][4] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][4] * f[ Nm3 * cst1 + cst2m ] + a[5][4] * f[ Nm4 * cst1 + cst2m ] + a[6][4] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][4] * f[ Nm6 * cst1 + cst2m ] + a[8][4] * f[ Nm7 * cst1 + cst2m ] + a[9][4] * f[ Nm8 * cst1 + cst2m ] +
						a[10][4] * f[ Nm9 * cst1 + cst2m ] + a[11][4] * f[ Nm10 * cst1 + cst2m ];

		cf[ Nm4 ] = a[1][5] * f[ N * cst1 + cst2m ] + a[2][5] * f[ Nm1 * cst1 + cst2m ] + a[3][5] * f[ Nm2 * cst1 + cst2m ] +
					   a[4][5] * f[ Nm3 * cst1 + cst2m ] + a[5][5] * f[ Nm4 * cst1 + cst2m ] + a[6][5] * f[ Nm5 * cst1 + cst2m ] +
					   a[7][5] * f[ Nm6 * cst1 + cst2m ] + a[8][5] * f[ Nm7 * cst1 + cst2m ] + a[9][5] * f[ Nm8 * cst1 + cst2m ] +
						a[10][5] * f[ Nm9 * cst1 + cst2m ] + a[11][5] * f[ Nm10 * cst1 + cst2m ];


		for ( k = 3; k <= Nm1; k++ )
		{
			km1 = k - 1;
			cf[ k ] -= cf[ km1 ] * bf[ k ] / df[ km1 ];
		}

		cst2n = cst3 + j * Nz * Nm + n;
		
		f[ Nm1 * cst1 + cst2n ] = cf[ Nm1 ] / df[ Nm1 ];

		for ( k = Nm2; k >= 2; k-- )
		{
			f[ k * cst1 + cst2n ] = ( cf[ k ] - af[ k ] * f[ ( k + 1 ) * cst1 + cst2n ] ) / df[ k ];
		}
	}
	}
}
