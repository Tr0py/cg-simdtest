#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>

/* size for A and B. */
#define SIZE 10

using namespace std;

const double NEARZERO = 1.0e-12;       // interpretation of "zero"

using vec    = vector<double>;         // vector
using matrix = vector<vec>;            // matrix (=collection of (row) vectors)

// Prototypes
void print( string title, const vec &V );
void print( string title, const matrix &A );
vec matrixTimesVector( const matrix &A, const vec &V );
vec vectorCombination( double a, const vec &U, double b, const vec &V );
double innerProduct( const vec &U, const vec &V );
double vectorNorm( const vec &V );
vec conjugateGradientSolver( const matrix &A, const vec &B );
matrix positiveDefiniteMatrix( const matrix &A);


//======================================================================


int main()
{
	int i, j;
	matrix A;
	vec B(SIZE), tmp(SIZE);
	clock_t startTime, endTime;

	// Fixed random seed.
	srand(0);

	// Initialize A and B
	for (i=0; i < SIZE; i++) {
		B[i] = rand() % 50;
	}
	for (i=0; i < SIZE; i++) {
		for (j=0; j < SIZE; j++) {
			tmp[j] = rand() % 5;
		}
		//print( "\ntmp:", tmp );
		A.push_back(tmp);
	}
	A = positiveDefiniteMatrix(A);


	// Start timer
	startTime = clock();
	vec X = conjugateGradientSolver( A, B );
	endTime = clock();

	cout << "Solves AX = B\n";
	print( "\nA:", A );
	print( "\nB:", B );
	print( "\nX:", X );

	print( "\nCheck AX:", matrixTimesVector( A, X ) );
	cout << "The run time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}


//======================================================================


void print( string title, const vec &V )
{
	cout << title << '\n';

	int n = V.size();
	for ( int i = 0; i < n; i++ )
	{
		double x = V[i];   if ( abs( x ) < NEARZERO ) x = 0.0;
		cout << x << '\t';
	}
	cout << '\n';
}


//======================================================================


void print( string title, const matrix &A )
{
	cout << title << '\n';

	int m = A.size(), n = A[0].size();                      // A is an m x n matrix
	for ( int i = 0; i < m; i++ )
	{
		for ( int j = 0; j < n; j++ )
		{
			double x = A[i][j];   if ( abs( x ) < NEARZERO ) x = 0.0;
			cout << x << '\t';
		}
		cout << '\n';
	}
}


//======================================================================


vec matrixTimesVector( const matrix &A, const vec &V )     // Matrix times vector
{
	int n = A.size();
	vec C( n );
	for ( int i = 0; i < n; i++ ) C[i] = innerProduct( A[i], V );
	return C;
}


//======================================================================


vec vectorCombination( double a, const vec &U, double b, const vec &V )        // Linear combination of vectors
{
	int n = U.size();
	vec W( n );
	for ( int j = 0; j < n; j++ ) W[j] = a * U[j] + b * V[j];
	return W;
}


//======================================================================


double innerProduct( const vec &U, const vec &V )          // Inner product of U and V
{
	return inner_product( U.begin(), U.end(), V.begin(), 0.0 );
}


//======================================================================


double vectorNorm( const vec &V )                          // Vector norm
{
	return sqrt( innerProduct( V, V ) );
}


//======================================================================


vec conjugateGradientSolver( const matrix &A, const vec &B )
{
	double TOLERANCE = 1.0e-12;

	int n = A.size();
	vec X( n, 0.0 );

	vec R = B;
	vec P = R;
	int k = 0;

	while ( 1 )
	{
		vec Rold = R;                                         // Store previous residual
		vec AP = matrixTimesVector( A, P );

		double alpha = innerProduct( R, R ) / innerProduct( P, AP );
		X = vectorCombination( 1.0, X, alpha, P );            // Next estimate of solution
		R = vectorCombination( 1.0, R, -alpha, AP );          // Residual

		if ( vectorNorm( R ) < TOLERANCE ) break;             // Convergence test

		double beta = innerProduct( R, R ) / innerProduct( Rold, Rold );
		P = vectorCombination( 1.0, R, beta, P );             // Next gradient
		k++;
	}

	return X;
}


//======================================================================
// ret = A * At
matrix positiveDefiniteMatrix( const matrix &A ) 
{
	int i, j, ii;
	matrix At(A), ret(A);
	

	// transpose A
	for (i=0; i < SIZE; ++i) {
		for (j=0; j < SIZE; ++j) {
			At[i][j] = A[j][i];
		}
	}
	// Initialize ret to zeros
	for (i=0; i < SIZE; ++i) {
		for (j=0; j < SIZE; ++j) {
			ret[i][j] = 0;
		}
	}
	// Matrix multiply
	for (i=0; i < SIZE; ++i) {
		for (j=0; j < SIZE; ++j) {
			for (ii=0; ii < SIZE; ++ii) {
				ret[i][j] += A[i][ii] * At[ii][j];
			}
		}
	}
	return ret;
}
