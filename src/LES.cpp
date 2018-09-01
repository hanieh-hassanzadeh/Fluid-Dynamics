//testing

#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <stdlib.h>

#include "dcox.h"
#include "dcoy.h"
#include "dcoz.h"
#include "f10x.h"
#include "f10y.h"
#include "f10z.h"
using namespace std;

const int n1ph = 85,
			 n2ph = 32,
			 n3ph = 32,
			 n1BZ = 8,
			 n2BZ = 0,
			 n3BZ = 0,
			 n1 = n1ph + n1BZ,
			 n2 = n2ph + n2BZ,
			 n3 = n3ph + n3BZ,
			 Nx = n1 + 1,
			 Ny = n2 + 1,
			 Nz = n3 + 1,
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 Nm = 51; //35//35//35//35/////////////////////////////adddssssssssssssssssssssssssssssssss//////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
const double pi = 3.141592654;

double U[ Nx ][ Ny ][ Nz ][ Nm ],
		 x[ Nx ][ Ny ][ Nz ],
		 y[ Nx ][ Ny ][ Nz ],
		 z[ Nx ][ Ny ][ Nz ],
		 Ja[ Nx ][ Ny ][ Nz ][ 10 ],
		 l1ph = 45,
		 l2ph = 30,
		 l3ph = 30,
		 l1BZ = 15,
		 l2BZ = 0,
		 l3BZ = 0,
		 l1 = l1ph + l1BZ,
		 l2 = l2ph + l2BZ,
		 l3 = l3ph + l3BZ,
		 dx = l1 / ( n1 - 1 ),
		 dy = l2 / ( n2 - 1 ),
		 dz = l3 / ( n3 - 1 ),
		 dmin = 10000000,
		 gama = 1.4, gama_m1 = 0.4,
		 Mi = 0.9,
		 Mn = 0,
		 Mt = 0,
		 ReL = 3600 / Mi,
		 RePr = ReL * 0.7,
		 nts = 1,
		 r0 = 1, b = 3.125;

dcox Uxi( &U[0][0][0][0], Nx, Ny, Nz, Nm, dx );
dcoy Uet( &U[0][0][0][0], Nx, Ny, Nz, Nm, dy );
dcoz Uze( &U[0][0][0][0], Nx, Ny, Nz, Nm, dz );
void SV();

f10x Ux_f( &U[0][0][0][0], Nx, Ny, Nz, Nm );
f10y Uy_f( &U[0][0][0][0], Nx, Ny, Nz, Nm );
f10z Uz_f( &U[0][0][0][0], Nx, Ny, Nz, Nm );

void grid_generation();
void grid_f();
void initialize();
void NS();
void F( int mro, int mrou1, int mrou2, int mrou3, int mE );
void out_put_data( double tmax );
void BC_nRef( int mro, int mrou1, int mrou2, int mrou3, int mE );
void Umean();


class stopwatch
	{
		public:
		stopwatch() : start(clock()){} //start counting time
		~stopwatch();
		private:
		clock_t start;
	};

int main()
{
	grid_generation();

//	initialize();
//	BC_nRef( 0, 1, 2, 3, 4 );
//	out_put_data( 0 );
	NS();

	cout << "			END				" << endl;
	int stop;
	cin >> stop;
	return 0;
}

void NS()
{
	double dt = 0.1 * dmin,
			 tmax = 10 * dt,
			 b1 = 0.53333333333333,
			 b2 = 0.25,
			 b3 = 0.416666666666667,
			 b4 = 0.25,
			 b5 = 0.75,
			 RERO, REROV, REROU, REROW, REE, J,
			 X, sigmax = 1;

	int i, j, k, ts = 0;

	Ux_f.setalpha( 0.49 );
	Uy_f.setalpha( 0.49 );
	Uz_f.setalpha( 0.49 );

	initialize();

	while ( ts * dt < tmax )
	{
		Umean();

//		cout << "11111111111111111111111111111111111111111111111111111111111111111111111" << endl;
		F( 0, 1, 2, 3, 4 );
		
		for ( i = 2; i < Nx-1; i++ )
			for ( j = 2; j < Ny-1; j++ )
				for ( k = 2; k < Nz-1; k++ )
				{
					J = Ja[ i ][ j ][ k ][ 0 ];
					RERO = ( U[ i ][ j ][ k ][ 5 ] + U[ i ][ j ][ k ][ 10 ] + U[ i ][ j ][ k ][ 15 ] ) * J;
					REROU = ( U[ i ][ j ][ k ][ 6 ] + U[ i ][ j ][ k ][ 11 ] + U[ i ][ j ][ k ][ 16 ] ) * J;
					REROV = ( U[ i ][ j ][ k ][ 7 ] + U[ i ][ j ][ k ][ 12 ] + U[ i ][ j ][ k ][ 17 ] ) * J;
					REROW = ( U[ i ][ j ][ k ][ 8 ] + U[ i ][ j ][ k ][ 13 ] + U[ i ][ j ][ k ][ 18 ] ) * J;
					REE = ( U[ i ][ j ][ k ][ 9 ] + U[ i ][ j ][ k ][ 14 ] + U[ i ][ j ][ k ][ 19 ] ) * J;

					U[ i ][ j ][ k ][ 5 ] = U[ i ][ j ][ k ][ 0 ] - dt * b1 * RERO;
					U[ i ][ j ][ k ][ 6 ] = U[ i ][ j ][ k ][ 1 ] - dt * b1 * REROU;
					U[ i ][ j ][ k ][ 7 ] = U[ i ][ j ][ k ][ 2 ] - dt * b1 * REROV;
					U[ i ][ j ][ k ][ 8 ] = U[ i ][ j ][ k ][ 3 ] - dt * b1 * REROW;
					U[ i ][ j ][ k ][ 9 ] = U[ i ][ j ][ k ][ 4 ] - dt * b1 * REE;


					b4 = 1;


					U[ i ][ j ][ k ][ 0 ] -= dt * b4 * RERO;
					U[ i ][ j ][ k ][ 1 ] -= dt * b4 * REROU;
					U[ i ][ j ][ k ][ 2 ] -= dt * b4 * REROV;
					U[ i ][ j ][ k ][ 3 ] -= dt * b4 * REROW;
					U[ i ][ j ][ k ][ 4 ] -= dt * b4 * REE;
				}
				
				if ( i > n1ph )
				{					
					X = sigmax * pow( ( x[ i ][ j ][ k ] - l1ph ) / ( l1 - l1ph ), 3 );
					RERO += X * ( U[ i ][ j ][ k ][ 0 ] - U[ i ][ j ][ k ][ 45 ] / nts );
					REROU += X * ( U[ i ][ j ][ k ][ 1 ] - U[ i ][ j ][ k ][ 46 ] / nts );
					REROV += X * ( U[ i ][ j ][ k ][ 2 ] - U[ i ][ j ][ k ][ 47 ] / nts );
					REROW += X * ( U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 48 ] / nts );
					REE += X * ( U[ i ][ j ][ k ][ 4 ] - U[ i ][ j ][ k ][ 49 ] / nts );
				}

/*		BC_nRef( 5, 6, 7, 8, 9 );
		cout << "22222222222222222222222222222222222222222222222222222222222222222222222222222222" << endl;
		F( 5, 6, 7, 8, 9 );
		
		for ( i = 2; i < Nx-1; i++ )
			for ( j = 2; j < Ny-1; j++ )
				for ( k = 2; k < Nz-1; k++ )
				{
					J = Ja[ i ][ j ][ k ][ 0 ];
					RERO = ( U[ i ][ j ][ k ][ 5 ] + U[ i ][ j ][ k ][ 10 ] + U[ i ][ j ][ k ][ 15 ] ) * J;
					REROU = ( U[ i ][ j ][ k ][ 6 ] + U[ i ][ j ][ k ][ 11 ] + U[ i ][ j ][ k ][ 16 ] ) * J;
					REROV = ( U[ i ][ j ][ k ][ 7 ] + U[ i ][ j ][ k ][ 12 ] + U[ i ][ j ][ k ][ 17 ] ) * J;
					REROW = ( U[ i ][ j ][ k ][ 8 ] + U[ i ][ j ][ k ][ 13 ] + U[ i ][ j ][ k ][ 18 ] ) * J;
					REE = ( U[ i ][ j ][ k ][ 9 ] + U[ i ][ j ][ k ][ 14 ] + U[ i ][ j ][ k ][ 19 ] ) * J;

					U[ i ][ j ][ k ][ 5 ] = U[ i ][ j ][ k ][ 0 ] - dt * b3 * RERO;
					U[ i ][ j ][ k ][ 6 ] = U[ i ][ j ][ k ][ 1 ] - dt * b3 * REROU;
					U[ i ][ j ][ k ][ 7 ] = U[ i ][ j ][ k ][ 2 ] - dt * b3 * REROV;
					U[ i ][ j ][ k ][ 8 ] = U[ i ][ j ][ k ][ 3 ] - dt * b3 * REROW;
					U[ i ][ j ][ k ][ 9 ] = U[ i ][ j ][ k ][ 4 ] - dt * b3 * REE;


				}

		BC_nRef( 5, 6, 7, 8, 9 );
		cout << "33333333333333333333333333333333333333333333333333333333333333333333333333333333" << endl;
		F( 5, 6, 7, 8, 9 );

		for ( i = 2; i < Nx-1; i++ )		
			for ( j = 2; j < Ny-1; j++ )
				for ( k = 2; k < Nz-1; k++ )
				{
					J = Ja[ i ][ j ][ k ][ 0 ];
					RERO = ( U[ i ][ j ][ k ][ 5 ] + U[ i ][ j ][ k ][ 10 ] + U[ i ][ j ][ k ][ 15 ] ) * J;
					REROU = ( U[ i ][ j ][ k ][ 6 ] + U[ i ][ j ][ k ][ 11 ] + U[ i ][ j ][ k ][ 16 ] ) * J;
					REROV = ( U[ i ][ j ][ k ][ 7 ] + U[ i ][ j ][ k ][ 12 ] + U[ i ][ j ][ k ][ 17 ] ) * J;
					REROW = ( U[ i ][ j ][ k ][ 8 ] + U[ i ][ j ][ k ][ 13 ] + U[ i ][ j ][ k ][ 18 ] ) * J;
					REE = ( U[ i ][ j ][ k ][ 9 ] + U[ i ][ j ][ k ][ 14 ] + U[ i ][ j ][ k ][ 19 ] ) * J;

//					cout << U[ n1 / 2 ][ n2 / 2 ][ n3 / 2 ][ 1 ]<< endl;
//					cout << REROU << endl;

					U[ i ][ j ][ k ][ 0 ] -= dt * b5 * RERO;
					U[ i ][ j ][ k ][ 1 ] -= dt * b5 * REROU;
					U[ i ][ j ][ k ][ 2 ] -= dt * b5 * REROV;
					U[ i ][ j ][ k ][ 3 ] -= dt * b5 * REROW;
					U[ i ][ j ][ k ][ 4 ] -= dt * b5 * REE;
				}

*/
		BC_nRef( 0, 1, 2, 3, 4 );

		if ( ts % 3 != 0 )
		{
//			cout << "00000000000000000000" << endl;
			Ux_f.setmn( 0, 0 );
			Ux_f.filter();

			Ux_f.setmn( 1, 1 );
			Ux_f.filter();

			Ux_f.setmn( 2, 2 );
			Ux_f.filter();

			Ux_f.setmn( 3, 3 );
			Ux_f.filter();

			Ux_f.setmn( 4, 4 );
			Ux_f.filter();
		}
		else if ( ts % 3 != 1 )
		{
//			cout << "1111111111111111111" << endl;
			Uy_f.setmn( 0, 0 );
			Uy_f.filter();

			Uy_f.setmn( 1, 1 );
			Uy_f.filter();

			Uy_f.setmn( 2, 2 );
			Uy_f.filter();

			Uy_f.setmn( 3, 3 );
			Uy_f.filter();

			Uy_f.setmn( 4, 4 );
			Uy_f.filter();
		}
		else
		{
//			cout << "2222222222222222222" << endl;
			Uz_f.setmn( 0, 0 );
			Uz_f.filter();

			Uz_f.setmn( 1, 1 );
			Uz_f.filter();

			Uz_f.setmn( 2, 2 );
			Uz_f.filter();

			Uz_f.setmn( 3, 3 );
			Uz_f.filter();

			Uz_f.setmn( 4, 4 );
			Uz_f.filter();
		}

		ts++;
		cout << ts << "\t" << ts * dt << endl;

//		if ( ts % 100 == 0 )
		{			
			out_put_data( tmax );
			SV();
			cout << "output\n";

//			if ( ts % 10000 == 0 )
//				File();
		}
	}

	clock_t total = clock();// - start; //get elapsed time
	cout<<"total of ticks for this activity: "<<total<<endl;
	cout<<"in seconds: "<<double(total)/CLOCKS_PER_SEC<<endl;

	out_put_data( tmax );
}

void F( int mro, int mrou1, int mrou2, int mrou3, int mE )
{
	int i, j, k;
	double ro, E, u1, u2, u3, p, EPp, u, v, w,
			 kxh, kyh, kzh, exh, eyh, ezh, zxh, zyh, zzh, Uh, Vh, Wh, rou, rov, row,
			 ux, uy, uz, vx, vy, vz, wx, wy, wz, Tx, Ty, Tz, psixx, psiyy, psizz, psixy, psixz, psiyz, Qx, Qy, Qz,
			 qx, qy, qz, Delta, Sxx, Sxy, Sxz, Syy, Syz, Szz, SMSM, SDxx, SDyy, SDzz, ReL_2, CC, CC1, CC2,
			 tauxx, tauxy, tauxz, tauyy, tauyz, tauzz, Fv2, Fv3, Fv4, Fv5, Gv2, Gv3, Gv4, Gv5, Hv2, Hv3, Hv4, Hv5,
			 kesixh, kesiyh, kesizh, ethaxh, ethayh, ethazh, zethaxh, zethayh, zethazh, Jac,
			 ukesi, uetha, uzetha, vkesi, vetha, vzetha, wkesi, wetha, wzetha, Tkesi, Tetha, Tzetha,
//			 Csgs = 0, CI = 0, Prt = 0.7;
			 Csgs = 0.012, CI = 0.0066, Prt = 0.7;

	for ( i = 1; i < Nx; i++ )
		for ( j = 1; j < Ny; j++ )
			for ( k = 1; k < Nz; k++ )
			{
				ro = U[ i ][ j ][ k ][ mro ];
				E = U[ i ][ j ][ k ][ mE ];
				rou = U[ i ][ j ][ k ][ mrou1 ];
				rov = U[ i ][ j ][ k ][ mrou2 ];
				row = U[ i ][ j ][ k ][ mrou3 ];
				u1 = rou / ro;
				u2 = rov / ro;
				u3 = row / ro;
				p = ( E - 0.5 * ro * ( u1 * u1 + u2 * u2 + u3 * u3 ) ) * gama_m1;
				EPp = E + p;

				kxh = Ja[ i ][ j ][ k ][ 1 ];
				kyh = Ja[ i ][ j ][ k ][ 2 ];
				kzh = Ja[ i ][ j ][ k ][ 3 ];

				exh = Ja[ i ][ j ][ k ][ 4 ];
				eyh = Ja[ i ][ j ][ k ][ 5 ];
				ezh = Ja[ i ][ j ][ k ][ 6 ];

				zxh = Ja[ i ][ j ][ k ][ 7 ];
				zyh = Ja[ i ][ j ][ k ][ 8 ];
				zzh = Ja[ i ][ j ][ k ][ 9 ];

				Uh = u1 * kxh + u2 * kyh + u3 * kzh;
				Vh = u1 * exh + u2 * eyh + u3 * ezh;
				Wh = u1 * zxh + u2 * zyh + u3 * zzh;

				U[ i ][ j ][ k ][ 5 ] = ro * Uh;
				U[ i ][ j ][ k ][ 6 ] = rou * Uh + p * kxh;
				U[ i ][ j ][ k ][ 7 ] = rov * Uh + p * kyh;
				U[ i ][ j ][ k ][ 8 ] = row * Uh + p * kzh;
				U[ i ][ j ][ k ][ 9 ] = EPp * Uh;

				U[ i ][ j ][ k ][ 10 ] = ro * Vh;
				U[ i ][ j ][ k ][ 11 ] = rou * Vh + p * exh;
				U[ i ][ j ][ k ][ 12 ] = rov * Vh + p * eyh;
				U[ i ][ j ][ k ][ 13 ] = row * Vh + p * ezh;
				U[ i ][ j ][ k ][ 14 ] = EPp * Vh;

				U[ i ][ j ][ k ][ 15 ] = ro * Wh;
				U[ i ][ j ][ k ][ 16 ] = rou * Wh + p * zxh;
				U[ i ][ j ][ k ][ 17 ] = rov * Wh + p * zyh;
				U[ i ][ j ][ k ][ 18 ] = row * Wh + p * zzh;
				U[ i ][ j ][ k ][ 19 ] = EPp * Wh;

				U[ i ][ j ][ k ][ 20 ] = u1;
				U[ i ][ j ][ k ][ 23 ] = u2;
				U[ i ][ j ][ k ][ 26 ] = u3;
				U[ i ][ j ][ k ][ 29 ] = gama * p / ro;

//				cout << U[ n1 / 2 ][ n2 / 2 ][ n3 / 2 ][ 14 ] << endl;
//				cout << E << endl;
			}

	Uze.derivate_co( 22, 20 );
	Uze.derivate_co( 25, 23 );
	Uze.derivate_co( 28, 26 );
	Uze.derivate_co( 31, 29 );

	Uet.derivate_co( 21, 20 );
	Uet.derivate_co( 24, 23 );
	Uet.derivate_co( 27, 26 );
	Uet.derivate_co( 30, 29 );

	Uxi.derivate_co( 20, 20 );
	Uxi.derivate_co( 23, 23 );
	Uxi.derivate_co( 26, 26 );
	Uxi.derivate_co( 29, 29 );

	for ( i = 1; i < Nx; i++ )
		for ( j = 1; j < Ny; j++ )
			for ( k = 1; k < Nz; k++ )
			{
				ro = U[ i ][ j ][ k ][ mro ];
				u = U[ i ][ j ][ k ][ mrou1 ] / ro;
				v = U[ i ][ j ][ k ][ mrou2 ] / ro;
				w = U[ i ][ j ][ k ][ mrou3 ] / ro;

				ukesi = U[ i ][ j ][ k ][ 20 ];
				uetha = U[ i ][ j ][ k ][ 21 ];
				uzetha = U[ i ][ j ][ k ][ 22 ];
				vkesi = U[ i ][ j ][ k ][ 23 ];
				vetha = U[ i ][ j ][ k ][ 24 ];
				vzetha = U[ i ][ j ][ k ][ 25 ];
				wkesi = U[ i ][ j ][ k ][ 26 ];
				wetha = U[ i ][ j ][ k ][ 27 ];
				wzetha = U[ i ][ j ][ k ][ 28 ];
				Tkesi = U[ i ][ j ][ k ][ 29 ];
				Tetha = U[ i ][ j ][ k ][ 30 ];
				Tzetha = U[ i ][ j ][ k ][ 31 ];

				kesixh = Ja[ i ][ j ][ k ][ 1 ];
				kesiyh = Ja[ i ][ j ][ k ][ 2 ];
				kesizh = Ja[ i ][ j ][ k ][ 3 ];
				ethaxh = Ja[ i ][ j ][ k ][ 4 ];
				ethayh = Ja[ i ][ j ][ k ][ 5 ];
				ethazh = Ja[ i ][ j ][ k ][ 6 ];
				zethaxh = Ja[ i ][ j ][ k ][ 7 ];
				zethayh = Ja[ i ][ j ][ k ][ 8 ];
				zethazh = Ja[ i ][ j ][ k ][ 9 ];
				Jac = Ja[ i ][ j ][ k ][ 0 ];
				Delta = 1.0 / pow( Jac, 1.0 / 3 );

				ux = ( ukesi * kesixh + uetha * ethaxh + uzetha * zethaxh ) * Jac;
				uy = ( ukesi * kesiyh + uetha * ethayh + uzetha * zethayh ) * Jac;
				uz = ( ukesi * kesizh + uetha * ethazh + uzetha * zethazh ) * Jac;

				vx = ( vkesi * kesixh + vetha * ethaxh + vzetha * zethaxh ) * Jac;
				vy = ( vkesi * kesiyh + vetha * ethayh + vzetha * zethayh ) * Jac;
				vz = ( vkesi * kesizh + vetha * ethazh + vzetha * zethazh ) * Jac;

				wx = ( wkesi * kesixh + wetha * ethaxh + wzetha * zethaxh ) * Jac;
				wy = ( wkesi * kesiyh + wetha * ethayh + wzetha * zethayh ) * Jac;
				wz = ( wkesi * kesizh + wetha * ethazh + wzetha * zethazh ) * Jac;

				Tx = ( Tkesi * kesixh + Tetha * ethaxh + Tzetha * zethaxh ) * Jac;
				Ty = ( Tkesi * kesiyh + Tetha * ethayh + Tzetha * zethayh ) * Jac;
				Tz = ( Tkesi * kesizh + Tetha * ethazh + Tzetha * zethazh ) * Jac;

				Sxx = ux;
				Sxy = 0.5 * ( uy + vx );
				Sxz = 0.5 * ( uz + wx );
				Syy = uy;
				Syz = 0.5 * ( vz + wy );
				Szz = wz;

				SMSM = 2 * ( Sxx * Sxx + Syy * Syy + Szz * Szz + 2 * ( Sxy * Sxy + Sxz * Sxz + Syz * Syz ) );
				SDxx = ( 2 * Sxx - Syy - Szz ) / 3;
				SDyy = ( 2 * Syy - Sxx - Szz ) / 3;
				SDzz = ( 2 * Szz - Sxx - Syy ) / 3;

				ReL_2 = 2 / ReL;
				
				psixx = ReL_2 * SDxx;
				psixy = ReL_2 * Sxy;
				psixz = ReL_2 * Sxz;				
				psiyy = ReL_2 * SDyy;
				psiyz = ReL_2 * Syz;
				psizz = ReL_2 * SDzz;

				CC = 2 * ro * Delta * Delta;
				CC1 = -CC * Csgs * sqrt( SMSM );
				CC2 = CC / 3 * CI * SMSM;

				tauxx = CC1 * SDxx + CC2;
				tauyy = CC1 * SDyy + CC2;
				tauzz = CC1 * SDzz + CC2;
				tauxy = CC1 * Sxy;
				tauxz = CC1 * Sxz;
				tauyz = CC1 * Syz;

				qx = -Tx / ( gama_m1 * RePr );
				qy = -Ty / ( gama_m1 * RePr );
				qz = -Tz / ( gama_m1 * RePr );

				Qx = CC1 / ( 2 * Prt ) * Tx;
				Qy = CC1 / ( 2 * Prt ) * Ty;
				Qz = CC1 / ( 2 * Prt ) * Tz;

				Fv2 = kesixh * ( psixx - tauxx ) + kesiyh * ( psixy - tauxy ) + kesizh * ( psixz - tauxz );
				Fv3 = kesixh * ( psixy - tauxy ) + kesiyh * ( psiyy - tauyy ) + kesizh * ( psiyz - tauyz );
				Fv4 = kesixh * ( psixz - tauxz ) + kesiyh * ( psiyz - tauyz ) + kesizh * ( psizz - tauzz );
				Fv5 = u * Fv2 + v * Fv3 + w * Fv4 - 
						kesixh * ( qx + Qx ) - kesiyh * ( qy + Qy ) - kesizh * ( qz + Qz );

				Gv2 = ethaxh * ( psixx - tauxx ) + ethayh * ( psixy - tauxy ) + ethazh * ( psixz - tauxz );
				Gv3 = ethaxh * ( psixy - tauxy ) + ethayh * ( psiyy - tauyy ) + ethazh * ( psiyz - tauyz );
				Gv4 = ethaxh * ( psixz - tauxz ) + ethayh * ( psiyz - tauyz ) + ethazh * ( psizz - tauzz );
				Gv5 = u * Gv2 + v * Gv3 + w * Gv4 - 
						ethaxh * ( qx + Qx ) - ethayh * ( qy + Qy ) - ethazh * ( qz + Qz );

				Hv2 = zethaxh * ( psixx - tauxx ) + zethayh * ( psixy - tauxy ) + zethazh * ( psixz - tauxz );
				Hv3 = zethaxh * ( psixy - tauxy ) + zethayh * ( psiyy - tauyy ) + zethazh * ( psiyz - tauyz );
				Hv4 = zethaxh * ( psixz - tauxz ) + zethayh * ( psiyz - tauyz ) + zethazh * ( psizz - tauzz );
				Hv5 = u * Hv2 + v * Hv3 + w * Hv4 - 
						zethaxh * ( qx + Qx ) - zethayh * ( qy + Qy ) - zethazh * ( qz + Qz );

				U[ i ][ j ][ k ][ 6 ] -= Fv2;
				U[ i ][ j ][ k ][ 7 ] -= Fv3;
				U[ i ][ j ][ k ][ 8 ] -= Fv4;
				U[ i ][ j ][ k ][ 9 ] -= Fv5;

				U[ i ][ j ][ k ][ 11 ] -= Gv2;
				U[ i ][ j ][ k ][ 12 ] -= Gv3;
				U[ i ][ j ][ k ][ 13 ] -= Gv4;
				U[ i ][ j ][ k ][ 14 ] -= Gv5;

				U[ i ][ j ][ k ][ 16 ] -= Hv2;
				U[ i ][ j ][ k ][ 17 ] -= Hv3;
				U[ i ][ j ][ k ][ 18 ] -= Hv4;
				U[ i ][ j ][ k ][ 19 ] -= Hv5;
//				cout << Tkesi << endl;
			}	

	Uxi.derivate_co( 5, 5 );
	Uxi.derivate_co( 6, 6 );
	Uxi.derivate_co( 7, 7 );
	Uxi.derivate_co( 8, 8 );
	Uxi.derivate_co( 9, 9 );

	Uet.derivate_co( 10, 10 );
	Uet.derivate_co( 11, 11 );
	Uet.derivate_co( 12, 12 );
	Uet.derivate_co( 13, 13 );
	Uet.derivate_co( 14, 14 );

	Uze.derivate_co( 15, 15 );
	Uze.derivate_co( 16, 16 );
	Uze.derivate_co( 17, 17 );
	Uze.derivate_co( 18, 18 );
	Uze.derivate_co( 19, 19 );
}

void BC_nRef( int mro, int mrou1, int mrou2, int mrou3, int mE )
{
	int i, j, k, ip, jp;
	double dro, dp, du, dv, c1, c2, c3, c4, P=1/gama, RHO=1, u=0, v=Mn, w = Mt, rou=Mi, rov=Mn, row = Mt;

//	int i, j, ip, jp;
//	double P=1/gama, u, v, E, ro,
//			 du, dv, dp, dro, c1, c2, c3, c4,
//			 epsilon, alfa = 0.045, um;
/*
	i = 1;
	ip = i + 1;
//	ip = 1;

	for ( j = 1; j <= n2; j++ )
		for ( k = 1; k <= n3; k++ )
		{
			epsilon = -1 + 2 * rand() / static_cast<double>( RAND_MAX );

			E = U[ ip ][ j ][ k ][ mE ];
			ro = U[ ip ][ j ][ k ][ mro ];
			u = U[ ip ][ j ][ k ][ mrou1 ] / ro;
			v = U[ ip ][ j ][ k ][ mrou2 ] / ro;
			w = U[ ip ][ j ][ k ][ mrou3 ] / ro;
			P = ( E - 0.5 * ro * ( u * u + v * v ) ) * gama_m1;
//			um = Uc + Um * tanh( 2 * y[ ip ][ j ][ k ] / 1.0 );

			du = u - um;
			dp = P - 1 / gama;

			c4 = dp - du;
				
			ro = 1 + 0.5 * c4;
			u = um - 0.5 * c4;
//			v = epsilon * alfa * Uc * exp( -y[ ip ][ j ][ k ] * y[ ip ][ j ][ k ] / dmin / dmin );
			P = 1 / gama + 0.5 * c4;

			U[ i ][ j ][ k ][ mro ] = ro;
			U[ i ][ j ][ k ][ mE ] = P / gama_m1 + 0.5 * ro * ( u * u + v * v );
			U[ i ][ j ][ k ][ mrou1 ] = ro * u;
			U[ i ][ j ][ k ][ mrou2 ] = ro * v;
			U[ i ][ j ][ k ][ mrou2 ] = ro * w;
		}

	i = n1;
	ip = i - 1;
//	ip = i;

	for ( j = 1; j <= n2; j++ )
		for ( k = 1; k <= n3; k++ )
		{
			E = U[ ip ][ j ][ k ][ mE ];
			ro = U[ ip ][ j ][ k ][ mro ];
			u = U[ ip ][ j ][ k ][ mrou1 ] / ro;
			v = U[ ip ][ j ][ k ][ mrou2 ] / ro;
			P = ( E - 0.5 * ro * ( u * u + v * v ) ) * gama_m1;
			um = u;

			dro = ro - 1;
			du = u - um;
			dv = v;
			dp = P - 1 / gama;

			c1 = dp - dro;
			c2 = dv;
			c3 = dp + du;
			c4 = 0;
	
			ro = 1 - c1 + 0.5 * c3;
			u = um + 0.5 * c3;
			v = c2;
			P = 1 / gama + 0.5 * c3;

			U[ i ][ j ][ k ][ mro ] = ro;
			U[ i ][ j ][ k ][ mE ] = P / gama_m1 + 0.5 * ro * ( u * u + v * v );
			U[ i ][ j ][ k ][ mrou1 ] = ro * u;
			U[ i ][ j ][ k ][ mrou2 ] = ro * v;
		}

	j = 1;
	jp = j + 1;
//	jp = j;

	for ( i = 1; i <= n1; i++ )
		for ( k = 1; k <= n3; k++ )
		{
			E = U[ i ][ jp ][ k ][ mE ];
			ro = U[ i ][ jp ][ k ][ mro ];
			u = U[ i ][ jp ][ k ][ mrou1 ] / ro;
			v = U[ i ][ jp ][ k ][ mrou2 ] / ro;
			P = ( E - 0.5 * ro * ( u * u + v * v ) ) * gama_m1;

			dv = v;
			dp = P - 1 / gama;
			c4 = dp - dv;
		
			ro = 1 + 0.5 * c4;
//			u = M2;
			v = -0.5 * c4;
			P = 1 / gama + 0.5 * c4;

			U[ i ][ j ][ k ][ mro ] = ro;
			U[ i ][ j ][ k ][ mE ] = P / gama_m1 + 0.5 * ro * ( u * u + v * v );
			U[ i ][ j ][ k ][ mrou1 ] = ro * u;
			U[ i ][ j ][ k ][ mrou2 ] = ro * v;
			U[ i ][ j ][ k ][ mrou3 ] = ro * w;
		}

	j = n2;
	jp = j - 1;
//	jp = j;

	for ( i = 1; i <= n1; i++ )
		for ( k = 1; k <= n3; k++ )
		{
			E = U[ i ][ jp ][ k ][ mE ];
			ro = U[ i ][ jp ][ k ][ mro ];
			u = U[ i ][ jp ][ k ][ mrou1 ] / ro;
			v = U[ i ][ jp ][ k ][ mrou2 ] / ro;
			P = ( E - 0.5 * ro * ( u * u + v * v ) ) * gama_m1;

			dro = ro - 1;
//			du = u - M1;
			dv = v - 0;
			dp = P - 1 / gama;

			c1 = dp - dro;
			c2 = du;
			c3 = dp + dv;
			c4 = 0;

			ro = 1 - c1 + 0.5 * c3;
//			u = M1 + c2;
			v = 0.5 * c3;
			P = 1 / gama + 0.5 * c3;

			U[ i ][ j ][ k ][ mro ] = ro;
			U[ i ][ j ][ k ][ mE ] = P / gama_m1 + 0.5 * ro * ( u * u + v * v );
			U[ i ][ j ][ k ][ mrou1 ] = ro * u;
			U[ i ][ j ][ k ][ mrou2 ] = ro * v;
		}
*/

	double r, ur;

	i = 1;

	for ( j = 1; j <= n2; j++ )
	for ( k = 1; k <= n3; k++ )
	{
		r = sqrt( z[ i ][ j ][ k ] * z[ i ][ j ][ k ] + y[ i ][ j ][ k ] * y[ i ][ j ][ k ] );
		ur = 0.5 * Mi * ( 1 - tanh( b * ( r / r0 - r0 / r ) ) );
//		ur = 0.5;
//		cout << ur << endl;
		v = w = 0;

		U[ 1 ][ j ][ k ][ mro ] = RHO;				// ro at left wall
		U[ 1 ][ j ][ k ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ 1 ][ j ][ k ][ mrou1 ] = RHO * ur;				
		U[ 1 ][ j ][ k ][ mrou2 ] = RHO * v;
		U[ 1 ][ j ][ k ][ mrou3 ] = RHO * w;

		U[ n1 ][ j ][ k ][ mro ] = RHO;				// ro at left wall
		U[ n1 ][ j ][ k ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ n1 ][ j ][ k ][ mrou1 ] = RHO * ur;				
		U[ n1 ][ j ][ k ][ mrou2 ] = RHO * v;
		U[ n1 ][ j ][ k ][ mrou3 ] = RHO * w;
	}

	for ( i = 1; i <= n1; i++ )
	for ( k = 1; k <= n3; k++ )
	{
		u = v = w = 0;
		U[ i ][ 1 ][ k ][ mro ] = RHO;				// ro at left wall
		U[ i ][ 1 ][ k ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ i ][ 1 ][ k ][ mrou1 ] = RHO * u;				
		U[ i ][ 1 ][ k ][ mrou2 ] = RHO * v;
		U[ i ][ 1 ][ k ][ mrou3 ] = RHO * w;

		U[ i ][ n2 ][ k ][ mro ] = RHO;				// ro at left wall
		U[ i ][ n2 ][ k ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ i ][ n2 ][ k ][ mrou1 ] = RHO * u;				
		U[ i ][ n2 ][ k ][ mrou2 ] = RHO * v;
		U[ i ][ n2 ][ k ][ mrou3 ] = RHO * w;
	}

	for ( i = 1; i <= n1; i++ )
	for ( j = 1; j <= n2; j++ )
	{
		u = v = w = 0;
		U[ i ][ j ][ 1 ][ mro ] = RHO;				// ro at left wall
		U[ i ][ j ][ 1 ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ i ][ j ][ 1 ][ mrou1 ] = RHO * u;				
		U[ i ][ j ][ 1 ][ mrou2 ] = RHO * v;
		U[ i ][ j ][ 1 ][ mrou3 ] = RHO * w;

		U[ i ][ j ][ n3 ][ mro ] = RHO;				// ro at left wall
		U[ i ][ j ][ n3 ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ i ][ j ][ n3 ][ mrou1 ] = RHO * u;				
		U[ i ][ j ][ n3 ][ mrou2 ] = RHO * v;
		U[ i ][ j ][ n3 ][ mrou3 ] = RHO * w;
	}

/*
	for ( j = 1; j <= n2; j++ )
	{
		RHO = U[ 2 ][ j ][ mro ];
		rou = U[ 2 ][ j ][ mrou1 ];
//		rou = 5.0 / 16 * rou1[ 1 ][ j ][ mrou1 ] + 15.0 / 16 * rou1[ 2 ][ j ][ mrou1 ] - 5.0 / 16 * rou1[ 3 ][ j ][ mrou1 ]
//																+1.0 / 16 * rou1[ 4 ][ j ][ mrou1 ];
		rov = U[ 2 ][ j ][ mrou2 ];
		P = ( U[ 2 ][ j ][ mE ] - 0.5 * ( rou * rou + rov * rov ) / RHO ) * gama_m1;

		dro = RHO - 1;
		du = rou / RHO - Mi;
		dv = rov / RHO - Mn;
		dp = P - 1 / gama;

		c1 = 0;//dp - dro;
		c2 = 0;//dv;
		c3 = 0;//du + dp;
		c4 = dp - du;

		P = 0.5 * ( c3 + c4 ) + 1 / gama;
		RHO = 0.5 * ( c3 + c4 ) - c1 + 1;
		u = 0.5 * ( c3 - c4 ) + Mi;
		v = c2 + Mn;

		U[ 1 ][ j ][ mro ] = RHO;				// ro at left wall
		U[ 1 ][ j ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v );
		U[ 1 ][ j ][ mrou1 ] = RHO * u;
		U[ 1 ][ j ][ mrou2 ] = RHO * v;
	}

	for ( j = 1; j <= n2; j++ )
	{
		RHO = U[ n1-1 ][ j ][ mro ];
		rou = U[ n1-1 ][ j ][ mrou1 ];
//		rou = 5.0 / 16 * rou1[ n1u1 ][ j ][ mrou1 ] + 15.0 / 16 * rou1[ n1 ][ j ][ mrou1 ]
//			 - 5.0 / 16 * rou1[ n1 - 1 ][ j ][ mrou1 ] + 1.0 / 16 * rou1[ n1 - 2 ][ j ][ mrou1 ];
		rov = U[ n1-1 ][ j ][ mrou2 ];
		P = ( U[ n1-1 ][ j ][ mE ] - 0.5 * ( rou * rou + rov * rov ) / RHO ) * gama_m1;

		dro = RHO - 1;
		du = rou / RHO - Mi;
		dv = rov / RHO - Mn;
		dp = P - 1 / gama;

		c1 = dp - dro;
		c2 = dv;
		c3 = du + dp;
		c4 = 0;//dp - du;

		P = 0.5 * ( c3 + c4 ) + 1 / gama;
		RHO = 0.5 * ( c3 + c4 ) - c1 + 1;
		u = 0.5 * ( c3 - c4 ) + Mi;
		v = c2 + Mn;

		U[ n1 ][ j ][ mro ] = RHO;				// ro at left wall
		U[ n1 ][ j ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v );
		U[ n1 ][ j ][ mrou1 ] = RHO * u;
		U[ n1 ][ j ][ mrou2 ] = RHO * v;
	}

	for ( i = 1; i <= n1; i++ )
	{
		RHO = U[ i ][ 2 ][ mro ];
		rov = U[ i ][ 2 ][ mrou2 ];
//		rou = 5.0 / 16 * rou1[ 1 ][ j ][ mrou1 ] + 15.0 / 16 * rou1[ 2 ][ j ][ mrou1 ] - 5.0 / 16 * rou1[ 3 ][ j ][ mrou1 ]
//																+1.0 / 16 * rou1[ 4 ][ j ][ mrou1 ];
		rou = U[ i ][ 2 ][ mrou1 ];
		P = ( U[ i ][ 2 ][ mE ] - 0.5 * ( rou * rou + rov * rov ) / RHO ) * gama_m1;

		dro = RHO - 1;
		du = rou / RHO - Mi;
		dv = rov / RHO - Mn;
		dp = P - 1 / gama;

//		if ( rov < 0 )
//		{
			c1 = 0;//dp - dro;
			c2 = 0;//du;
			c3 = 0;//dp + dv;
			c4 =  dp - dv;
//		}
//		else
//		{
//			c1 = dp - dro;
//			c2 = du;
//			c3 = 0;//dp + dv;
//			c4 = dp - dv;
//		}

		P = 0.5 * ( c3 + c4 ) + 1 / gama;
		RHO = 0.5 * ( c3 + c4 ) - c1 + 1;
		v = 0.5 * ( c3 - c4 ) + Mn;
		u = c2 + Mi;

		U[ i ][ 1 ][ mro ] = RHO;				// ro at left wall
		U[ i ][ 1 ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v );
		U[ i ][ 1 ][ mrou1 ] = RHO * u;
		U[ i ][ 1 ][ mrou2 ] = RHO * v;
	}

	for ( i = 1; i <= n1; i++ )
	{
		RHO = U[ i ][ n2-1 ][ mro ];
		rov = U[ i ][ n2-1 ][ mrou2 ];
//		rou = 5.0 / 16 * rou1[ 1 ][ j ][ mrou1 ] + 15.0 / 16 * rou1[ 2 ][ j ][ mrou1 ] - 5.0 / 16 * rou1[ 3 ][ j ][ mrou1 ]
//																+1.0 / 16 * rou1[ 4 ][ j ][ mrou1 ];
		rou = U[ i ][ n2-1 ][ mrou1 ];
		P = ( U[ i ][ n2-1 ][ mE ] - 0.5 * ( rou * rou + rov * rov ) / RHO ) * gama_m1;

		dro = RHO - 1;
		du = rou / RHO - Mi;
		dv = rov / RHO - Mn;
		dp = P - 1 / gama;

//		if ( rov < 0 )
//		{
			c1 = dp - dro;
			c2 = du;
			c3 = dp + dv;
			c4 = 0;//dp - dv;
//		}
//		else
//		{
//			c1 = dp - dro;
//			c2 = du;
//			c3 = 0;//dp + dv;
//			c4 = dp - dv;
//		}

		P = 0.5 * ( c3 + c4 ) + 1 / gama;
		RHO = 0.5 * ( c3 + c4 ) - c1 + 1;
		v = 0.5 * ( c3 - c4 ) + Mn;
		u = c2 + Mi;

		U[ i ][ n2 ][ mro ] = RHO;				// ro at left wall
		U[ i ][ n2 ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v );
		U[ i ][ n2 ][ mrou1 ] = RHO * u;
		U[ i ][ n2 ][ mrou2 ] = RHO * v;
	}
*/
/*
	for ( j = 1; j <= n2; j++ )
	for ( k = 1; k <= n3; k++ )
	{
		U[ 1 ][ j ][ k ][ mro ] = RHO;				// ro at left wall
		U[ 1 ][ j ][ k ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ 1 ][ j ][ k ][ mrou1 ] = RHO * u;				
		U[ 1 ][ j ][ k ][ mrou2 ] = RHO * v;
		U[ 1 ][ j ][ k ][ mrou3 ] = RHO * w;

		U[ n1 ][ j ][ k ][ mro ] = RHO;				// ro at left wall
		U[ n1 ][ j ][ k ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ n1 ][ j ][ k ][ mrou1 ] = RHO * u;				
		U[ n1 ][ j ][ k ][ mrou2 ] = RHO * v;
		U[ n1 ][ j ][ k ][ mrou3 ] = RHO * w;
	}

	for ( i = 1; i <= n1; i++ )
	for ( k = 1; k <= n3; k++ )
	{
		U[ i ][ 1 ][ k ][ mro ] = RHO;				// ro at left wall
		U[ i ][ 1 ][ k ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ i ][ 1 ][ k ][ mrou1 ] = RHO * u;				
		U[ i ][ 1 ][ k ][ mrou2 ] = RHO * v;
		U[ i ][ 1 ][ k ][ mrou3 ] = RHO * w;

		U[ i ][ n2 ][ k ][ mro ] = RHO;				// ro at left wall
		U[ i ][ n2 ][ k ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ i ][ n2 ][ k ][ mrou1 ] = RHO * u;				
		U[ i ][ n2 ][ k ][ mrou2 ] = RHO * v;
		U[ i ][ n2 ][ k ][ mrou3 ] = RHO * w;
	}

	for ( i = 1; i <= n1; i++ )
	for ( j = 1; j <= n2; j++ )
	{
		U[ i ][ j ][ 1 ][ mro ] = RHO;				// ro at left wall
		U[ i ][ j ][ 1 ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ i ][ j ][ 1 ][ mrou1 ] = RHO * u;				
		U[ i ][ j ][ 1 ][ mrou2 ] = RHO * v;
		U[ i ][ j ][ 1 ][ mrou3 ] = RHO * w;

		U[ i ][ j ][ n3 ][ mro ] = RHO;				// ro at left wall
		U[ i ][ j ][ n3 ][ mE ] = P / gama_m1 + 0.5 * RHO * ( u * u + v * v + w * w );
		U[ i ][ j ][ n3 ][ mrou1 ] = RHO * u;				
		U[ i ][ j ][ n3 ][ mrou2 ] = RHO * v;
		U[ i ][ j ][ n3 ][ mrou3 ] = RHO * w;
	}
*/
}

void grid_generation()
{
	int i, j, k;

	double deltax = l1 / ( n1 - 1 ),
			 deltay = l2 / ( n2 - 1 ),
			 deltaz = l3 / ( n3 - 1 ),
			 sig, s, st, del = 0.00001,
			 ds = l1ph / ( n1ph - 1 ),
			 dxmax = 7 * deltax,
			 smax = ( n1 - 1 ) * ds,
			 xp = l1ph, xmax1 = l1, xt,
			 A, alf = 0.5, beta = 7,
			 xkesi, ykesi, zkesi, xetha, yetha, zetha,  xzetha, yzetha, zzetha;

	A = log( ( 1 + ( exp( beta ) - 1 ) * alf ) / ( 1 + ( exp( -beta ) - 1 ) * alf ) ) / ( 2 * beta );

	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
//				x[ i ][ j ] = deltax * ( i - 1 );
				y[ i ][ j ][ k ] = alf * (l2) * ( 1 + sinh( beta * ( deltay / l2 * ( j - 1 ) - A ) ) / sinh( beta * A ) ) - l2 / 2;
				z[ i ][ j ][ k ] = alf * (l3) * ( 1 + sinh( beta * ( deltaz / l3 * ( k - 1 ) - A ) ) / sinh( beta * A ) ) - l3 / 2;


				s = ds * ( i - 1 );
				st = ( smax * ( 1 + dxmax / ds ) - xmax1 ) / ( dxmax / ds );
				sig = log( dxmax / ( del * ds ) ) / ( st - xp );
				x[ i ][ j ][ k ] = s + dxmax / ( sig * ds ) * log( exp( sig * ( s - st ) ) + 1 );

				U[ i ][ j ][ k ][ 0 ] = x[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 1 ] = y[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 2 ] = z[ i ][ j ][ k ];
			}
		}
	}

/*	j = 1;

	for ( i = 1; i <= n1; i++ )
	{
		if ( i != n1 )
		{
			if ( fabs( ( x[ i + 1 ][ 1 ] - x[ i ][ 1 ] ) - 1.25 * ds ) < pow( 10, -2 ) )
			{				
				xt = x[ i ][ j ];
				cout << "xt = " << "\t" << xt << "at i = " << i << endl;
			}
		}
	}	
*/
/*
	for ( i = 1; i <= n1; i++ )
	{
		a[ i ] = 5.0 / 8.0 + 3.0 / 8.0 / ( 1 + exp( log10( del ) * ( x[ i ][ j ] - xt ) / ( xp - xt ) ) );
		b[ i ] = 2.0 / 3.0 * ( 1 - a[ i ] );
		c[ i ] = 1.0 / 6.0 * ( a[ i ] - 1 );

		if ( i > n1ph )
			cout << i << "\t" << a[ i ] << endl;
	}
*/
	double xmin = -10.0, xmax = 10.0, ymin = -10.0, ymax = 10.0,
			 Ax = 3, Ay = 3, etx = 3, ety = 3, wtao = 0.25, phi = 0 * pi / 2,
			 dxmin=100000, dymin=100000, dmin2=100000,
			 dxx = 1.0 / n1, dyy = 1.0 / n2;
			 
/*
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			x[ i ][ j ] = xmin + dxx * l1 * ( ( i - 1 ) + Ax * sin( 2 * pi * wtao ) * 
											sin( etx * pi * ( j - 1 ) * dyy * l2 / l2 ) );
			y[ i ][ j ] = ymin + dyy * l2 * ( ( j - 1 ) + Ay * sin( 2 * pi * wtao ) * 
											sin( ety * pi * ( i - 1 ) * dxx * l1 / l1 ) );

			U[ i ][ j ][ 0 ] = x[ i ][ j ];
			U[ i ][ j ][ 1 ] = y[ i ][ j ];
		}
	}
*/

	i = 1;
	j = 1;
	dmin = fabs( z[ i ][ 1 ][ 2 ] - z[ i ][ 1 ][ 1 ] );
	for ( k = 1; k <= n3; k++ )
	{
//		for ( j = 1; j <= n2; j++ )
//		{
			if ( k != 1 )
			{
				dxmin = fabs( z[ i ][ j ][ k ] - z[ i ][ j ][ k - 1 ] );
			}

//			if ( j != 1 )
//			{
//				dymin = fabs( y[ i ][ j ][ k ] - y[ i ][ j - 1 ][ k ] );
//			}

//			dmin = dxmin;
//			if ( dmin2 > dymin )
//				dmin2 = dymin;

			if ( dmin > dxmin )
				dmin = dxmin;
//		}
	}

//	dmin = 
	cout << "dmin = " << dmin << endl;

	grid_f();

///////////////////////  Jacobian   //////////////////////////////////////
	Uxi.derivate_co( 3, 0 );
	Uet.derivate_co( 4, 0 );
	Uze.derivate_co( 5, 0 );

	Uxi.derivate_co( 6, 1 );
	Uet.derivate_co( 7, 1 );
	Uze.derivate_co( 8, 1 );

	Uxi.derivate_co( 9, 2 );
	Uet.derivate_co( 10, 2 );
	Uze.derivate_co( 11, 2 );

	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				xkesi = U[ i ][ j ][ k ][ 3 ];
				xetha = U[ i ][ j ][ k ][ 4 ];
				xzetha = U[ i ][ j ][ k ][ 5 ];

				ykesi = U[ i ][ j ][ k ][ 6 ];
				yetha = U[ i ][ j ][ k ][ 7 ];
				yzetha = U[ i ][ j ][ k ][ 8 ];

				zkesi = U[ i ][ j ][ k ][ 9 ];
				zetha = U[ i ][ j ][ k ][ 10 ];
				zzetha = U[ i ][ j ][ k ][ 11 ];

				Ja[ i ][ j ][ k ][ 0 ] = 1.0 / ( xkesi * ( yetha * zzetha - yzetha * zetha ) - 
														   xetha * ( ykesi * zzetha - yzetha * zkesi ) +
															xzetha * ( ykesi * zetha - yetha * zkesi ) );
			}
		}
	}
////////////////////////////////////////////////////////////////////////////

///////////////////////  kesix   //////////////////////////////////////
	Uet.derivate_co( 3, 1 );
	Uze.derivate_co( 4, 1 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				U[ i ][ j ][ k ][ 3 ] *= z[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 4 ] *= z[ i ][ j ][ k ];
			}
		}
	}
	Uze.derivate_co( 3, 3 );
	Uet.derivate_co( 4, 4 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				Ja[ i ][ j ][ k ][ 1 ] = U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 4 ];
			}
		}
	}
////////////////////////////////////////////////////////////////////////////

///////////////////////  kesiy   //////////////////////////////////////
	Uet.derivate_co( 3, 2 );
	Uze.derivate_co( 4, 2 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				U[ i ][ j ][ k ][ 3 ] *= x[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 4 ] *= x[ i ][ j ][ k ];
			}
		}
	}
	Uze.derivate_co( 3, 3 );
	Uet.derivate_co( 4, 4 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				Ja[ i ][ j ][ k ][ 2 ] = U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 4 ];
			}
		}
	}
////////////////////////////////////////////////////////////////////////////

///////////////////////  kesiz   //////////////////////////////////////
	Uet.derivate_co( 3, 0 );
	Uze.derivate_co( 4, 0 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				U[ i ][ j ][ k ][ 3 ] *= y[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 4 ] *= y[ i ][ j ][ k ];
			}
		}
	}
	Uze.derivate_co( 3, 3 );
	Uet.derivate_co( 4, 4 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				Ja[ i ][ j ][ k ][ 3 ] = U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 4 ];
			}
		}
	}
////////////////////////////////////////////////////////////////////////////


///////////////////////  ethax   //////////////////////////////////////
	Uze.derivate_co( 3, 1 );
	Uxi.derivate_co( 4, 1 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				U[ i ][ j ][ k ][ 3 ] *= z[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 4 ] *= z[ i ][ j ][ k ];
			}
		}
	}
	Uxi.derivate_co( 3, 3 );
	Uze.derivate_co( 4, 4 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				Ja[ i ][ j ][ k ][ 4 ] = U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 4 ];
			}
		}
	}
////////////////////////////////////////////////////////////////////////////

///////////////////////  ethay   //////////////////////////////////////
	Uze.derivate_co( 3, 2 );
	Uxi.derivate_co( 4, 2 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				U[ i ][ j ][ k ][ 3 ] *= x[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 4 ] *= x[ i ][ j ][ k ];
			}
		}
	}
	Uxi.derivate_co( 3, 3 );
	Uze.derivate_co( 4, 4 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				Ja[ i ][ j ][ k ][ 5 ] = U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 4 ];
			}
		}
	}
////////////////////////////////////////////////////////////////////////////

///////////////////////  ethaz   //////////////////////////////////////
	Uze.derivate_co( 3, 0 );
	Uxi.derivate_co( 4, 0 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				U[ i ][ j ][ k ][ 3 ] *= y[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 4 ] *= y[ i ][ j ][ k ];
			}
		}
	}
	Uxi.derivate_co( 3, 3 );
	Uze.derivate_co( 4, 4 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				Ja[ i ][ j ][ k ][ 6 ] = U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 4 ];
			}
		}
	}
////////////////////////////////////////////////////////////////////////////


///////////////////////  zethax   //////////////////////////////////////
	Uxi.derivate_co( 3, 1 );
	Uet.derivate_co( 4, 1 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				U[ i ][ j ][ k ][ 3 ] *= z[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 4 ] *= z[ i ][ j ][ k ];
			}
		}
	}
	Uet.derivate_co( 3, 3 );
	Uxi.derivate_co( 4, 4 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				Ja[ i ][ j ][ k ][ 7 ] = U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 4 ];
			}
		}
	}
////////////////////////////////////////////////////////////////////////////

///////////////////////  zethay   //////////////////////////////////////
	Uxi.derivate_co( 3, 2 );
	Uet.derivate_co( 4, 2 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				U[ i ][ j ][ k ][ 3 ] *= x[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 4 ] *= x[ i ][ j ][ k ];
			}
		}
	}
	Uet.derivate_co( 3, 3 );
	Uxi.derivate_co( 4, 4 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				Ja[ i ][ j ][ k ][ 8 ] = U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 4 ];
			}
		}
	}
////////////////////////////////////////////////////////////////////////////

///////////////////////  zethaz   //////////////////////////////////////
	Uxi.derivate_co( 3, 0 );
	Uet.derivate_co( 4, 0 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				U[ i ][ j ][ k ][ 3 ] *= y[ i ][ j ][ k ];
				U[ i ][ j ][ k ][ 4 ] *= y[ i ][ j ][ k ];
			}
		}
	}
	Uet.derivate_co( 3, 3 );
	Uxi.derivate_co( 4, 4 );
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
				Ja[ i ][ j ][ k ][ 9 ] = U[ i ][ j ][ k ][ 3 ] - U[ i ][ j ][ k ][ 4 ];
			}
		}
	}
////////////////////////////////////////////////////////////////////////////
}

void initialize()
{
	int i, j, k;

	double r2, p, R = 1, alf1 = log( 2 ) / 1, u, v, w, dp = 0.001, RHO,
			 xc = l1 / 2,
			 yc = l2 / 2,
			 zc = l3 / 2;
	double r, ur;

	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
		{
			for ( k = 1; k <= n3; k++ )
			{
/*				r2 = ( ( x[ i ][ j ][ k ] - xc ) * ( x[ i ][ j ][ k ] - xc ) + ( y[ i ][ j ][ k ] - yc ) * ( y[ i ][ j ][ k ] - yc ) + 
						 ( z[ i ][ j ][ k ] - zc ) * ( z[ i ][ j ][ k ] - zc ) ) / ( R * R );
*/				p = 1 / gama;// + dp * exp( -r2 * alf1 );

				r = sqrt( y[ i ][ j ][ k ] * y[ i ][ j ][ k ] + z[ i ][ j ][ k ] * z[ i ][ j ][ k ] ) + pow( 10, -5 );
				ur = 0.5 * Mi * ( 1 - tanh( b * ( r / r0 - r0 / r ) ) );
//				ur = 0.5;
//				cout << ur << endl;

				RHO = 1;// + dp * exp( -r2 * alf1 );

				U[ i ][ j ][ k ][ 0 ] = RHO;
				u = ur;
				v = 0;
				w = 0;
				U[ i ][ j ][ k ][ 1 ] = U[ i ][ j ][ k ][ 0 ] * u;
				U[ i ][ j ][ k ][ 2 ] = U[ i ][ j ][ k ][ 0 ] * v;
				U[ i ][ j ][ k ][ 3 ] = U[ i ][ j ][ k ][ 0 ] * w;
				U[ i ][ j ][ k ][ 4 ] = p / ( gama_m1 ) + 0.5 * U[ i ][ j ][ k ][ 0 ] * ( u * u + v * v + w * w );


				U[ i ][ j ][ k ][ 39 ] = 0;
				U[ i ][ j ][ k ][ 40 ] = 0;
				U[ i ][ j ][ k ][ 41 ] = 0;

		
			U[ i ][ j ][ k ][ 42 ] = 0;
			U[ i ][ j ][ k ][ 43 ] = 0;
			U[ i ][ j ][ k ][ 44 ] = 0;

			U[ i ][ j ][ k ][ 45 ] = 0;
			U[ i ][ j ][ k ][ 46 ] = 0;
			U[ i ][ j ][ k ][ 47 ] = 0;
			U[ i ][ j ][ k ][ 48 ] = 0;
			U[ i ][ j ][ k ][ 49 ] = 0;
			}
		}
	}
}

void out_put_data( double tmax )
{
	ofstream fout( "data.plt" );

	if(!fout)
	{
		cout << "Can't open file.\n";
	}

	int i, j, k;

	double p;

	fout << "VARIABLES=x,y,z,ro,p,u,v,w,E" << endl;
	fout << "\nZONE I=" << Nx - 1 << " J=" << Ny - 1 << " K=" << Nz - 1 << endl;

	
	for ( k = 1; k < Nz; k++ )
	for ( j = 1; j < Ny; j++ )
		for ( i = 1; i < Nx; i++ )
	{
		
	
				p = ( U[ i ][ j ][ k ][ 4 ] - 0.5 * ( U[ i ][ j ][ k ][ 1 ] * U[ i ][ j ][ k ][ 1 ] + 
						U[ i ][ j ][ k ][ 2 ] * U[ i ][ j ][ k ][ 2 ] + U[ i ][ j ][ k ][ 3 ] * U[ i ][ j ][ k ][ 3 ] )
						  / U[ i ][ j ][ k ][ 0 ] ) * gama_m1;

			fout << setprecision( 20 ) << x[ i ][ j ][ k ] << "\t\t" << y[ i ][ j ][ k ] << "\t\t" << z[ i ][ j ][ k ]
				  << "\t\t" << U[ i ][ j ][ k ][ 0 ] << "\t\t" << ( U[ i ][ j ][ k ][ 4 ] - 0.5 * ( U[ i ][ j ][ k ][ 1 ] * U[ i ][ j ][ k ][ 1 ] + 
						U[ i ][ j ][ k ][ 2 ] * U[ i ][ j ][ k ][ 2 ] + U[ i ][ j ][ k ][ 3 ] * U[ i ][ j ][ k ][ 3 ] )
						  / U[ i ][ j ][ k ][ 0 ] ) * gama_m1
				  << "\t\t" << U[ i ][ j ][ k ][ 1 ] / U[ i ][ j ][ k ][ 0 ]
				  << "\t\t" << U[ i ][ j ][ k ][ 2 ] / U[ i ][ j ][ k ][ 0 ]
				  << "\t\t" << U[ i ][ j ][ k ][ 3 ] / U[ i ][ j ][ k ][ 0 ]
				  << "\t\t" << U[ i ][ j ][ k ][ 4 ] << endl;

//				fout << setprecision( 30 ) << x[ i ][ j ][ k ] << "\t\t\t" << y[ i ][ j ][ k ] << "\t\t\t" << z[ i ][ j ][ k ] 
//					  << "\t\t\t" <<  U[ i ][ j ][ k ][ 0 ] << "\t\t\t"	<< p << "\t\t\t" 
//					  << U[ i ][ j ][ k ][ 1 ] / U[ i ][ j ][ k ][ 0 ] << "\t\t\t" << U[ i ][ j ][ k ][ 2 ] / U[ i ][ j ][ k ][ 0 ]
//					  << "\t\t\t" << U[ i ][ j ][ k ][ 3 ] / U[ i ][ j ][ k ][ 0 ] << "\t\t\t" << U[ i ][ j ][ k ][ 4 ] << endl;
		
	}

	fout.close();
}

void grid_f()
{
	ofstream fout( "grid.plt" );

	if(!fout)
	{
		cout << "Can't open file.\n";
	}

	int i, j, k;

	fout << "VARIABLES=x,y,z,xx" << endl;
	fout << "\nZONE I=" << Nx - 1 << " J=" << Ny - 1 << " K=" << Nz - 1 << endl;

	for ( k = 1; k < Nz; k++ )
	for ( j = 1; j < Ny; j++ )
	{
		for ( i = 1; i < Nx; i++ )
		{
			fout << setprecision( 20 ) << x[ i ][ j ][ k ] << "\t\t" << y[ i ][ j ][ k ] 
				<< "\t\t" << z[ i ][ j ][ k ] << "\t\t" << x[i][j][k] << endl;
		}
	}

	fout.close();
}





void Umean()
{
	int i, j, k;
	double u, v, w, um, vm, wm;
	
	for ( i = 1; i <= n1; i++ )
	{
		for ( j = 1; j <= n2; j++ )
			for ( k = 1; k <= n3; k++ )
		{
			u = U[ i ][ j ][ k ][ 1 ] / U[ i ][ j ][ k ][ 0 ];
			v = U[ i ][ j ][ k ][ 2 ] / U[ i ][ j ][ k ][ 0 ];
			w = U[ i ][ j ][ k ][ 3 ] / U[ i ][ j ][ k ][ 0 ];

			U[ i ][ j ][ k ][ 39 ] += u;
			U[ i ][ j ][ k ][ 40 ] += v;
			U[ i ][ j ][ k ][ 41 ] += w;

			um = U[ i ][ j ][ k ][ 39 ] / nts;
			vm = U[ i ][ j ][ k ][ 40 ] / nts;
			wm = U[ i ][ j ][ k ][ 41 ] / nts;

			U[ i ][ j ][ k ][ 42 ] += ( u - um ) * ( u - um );
			U[ i ][ j ][ k ][ 43 ] += ( v - vm ) * ( v - vm );
			U[ i ][ j ][ k ][ 44 ] += ( u - um ) * ( v - vm );

			U[ i ][ j ][ k ][ 45 ] += U[ i ][ j ][ k ][ 0 ];
			U[ i ][ j ][ k ][ 46 ] += U[ i ][ j ][ k ][ 1 ];
			U[ i ][ j ][ k ][ 47 ] += U[ i ][ j ][ k ][ 2 ];
			U[ i ][ j ][ k ][ 48 ] += U[ i ][ j ][ k ][ 3 ];
			U[ i ][ j ][ k ][ 49 ] += U[ i ][ j ][ k ][ 4 ];
		}
	}	
}

void SV()
{
	int i, j, k;
	double umean, sigxx, sigyy, sigxy;

//	cout << "nts = " << nts << endl;
	
	ofstream fout( "SV.plt" );

	fout << "VARIABLES=x,y,z,umean,sigxx,sigyy,sigxy" << endl;
	fout << "\nZONE I=" << Nx - 1 << " J=" << Ny - 1 << " K=" << Nz - 1 << endl;

	for ( k = 1; k <= n3; k++ )
	for ( j = 1; j < Ny; j++ )
	{
		for ( i = 1; i < Nx; i++ )
		{
			umean = U[ i ][ j ][ k ][ 39 ] / nts;

			sigxx = U[ i ][ j ][ k ][ 42 ] / nts;
			sigyy = U[ i ][ j ][ k ][ 43 ] / nts;
			sigxy = U[ i ][ j ][ k ][ 44 ] / nts;

			fout << x[ i ][ j ][ k ] << "\t\t" << y[ i ][ j ][ k ] << "\t\t" << z[ i ][ j ][ k ] << "\t\t" << umean << "\t\t" << sigxx << "\t\t"
				  << sigyy << "\t\t" << sigxy << endl;
		}
	}

	fout.close();
}
