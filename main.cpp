//
//  main.cpp
//  simd1054
//
//  Created by Jane on 4/9/21.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <x86intrin.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <stdlib.h>
#include <xmmintrin.h>
#include <pmmintrin.h>

/* size for A and B. */
#ifndef SIZE
#define SIZE 256
#endif
#ifndef LOOPTIME
#define LOOPTIME 10
#endif
//#define SIMD_OPT

using namespace std; //用std

const float NEARZERO = 1.0e-12;       // interpretation of "zero"


//using vec    = vector<float>;         // vector using=typedef
typedef vector<float> vec;

//using matrix = vector<vec>;            // matrix (=collection of (row) vectors)
typedef vector<vec> matrix;


// Prototypes
void print( string title, const vec &V );
void print( string title, const matrix &A );
vec matrixTimesVector( const matrix &A, const vec &V );
vec vectorCombination( float a, const vec &U, float b, const vec &V );
float innerProduct( const vec &U, const vec &V );
float vectorNorm( const vec &V ); //模，平方相加开根号
vec conjugateGradientSolver( const matrix &A, const vec &B );
matrix positiveDefiniteMatrix( const matrix &A); //正定 A*A转置


//======================================================================


int main() //初始化A b，cg解x（计时）
{
	int i, j;
	matrix A;
	vec B(SIZE), tmp(SIZE); //tmp初始化
	clock_t startTime, endTime;
	vec X[LOOPTIME];
	// Fixed random seed.
	srand(time(NULL));

	// Initialize A and B
	for (i=0; i < SIZE; i++) {
		B[i] = rand() % 50;
	}
	for (i=0; i < SIZE; i++) {
		for (j=0; j < SIZE; j++) {
			tmp[j] = rand() % 5;
		}
		//print( "\ntmp:", tmp );
		A.push_back(tmp); //每次push一行，一开始A是0*0，push第一次就变成1*10，10次就是10*10
	}
	A = positiveDefiniteMatrix(A); //正定


	// Start timer
	startTime = clock();
	for (int i=0; i < LOOPTIME; i++) {
		X[i] = conjugateGradientSolver( A, B );
	}
	endTime = clock();

	//cout << "Solves AX = B\n";
	//print( "\nA:", A );
	print( "\nB:", B );
	//print( "\nX:", X );

	print( "\nCheck AX:", matrixTimesVector( A, X[0] ) ); //检查
	cout << "The run time is: " <<(float)(endTime - startTime) / CLOCKS_PER_SEC << " s" << endl;
}


//======================================================================


void print( string title, const vec &V )
{
	cout << title << '\n';

	int n = V.size();
	for ( int i = 0; i < n; i++ )
	{
		float x = V[i];   if ( abs( x ) < NEARZERO ) x = 0.0; //一个非常小的数，就记为0
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
			float x = A[i][j];   if ( abs( x ) < NEARZERO ) x = 0.0;
			cout << x << '\t';
		}
		cout << '\n';
	}
}


//======================================================================


vec matrixTimesVector( const matrix &A, const vec &V )     // Matrix times vector 矩阵*向量
{
	int n = A.size();
	vec C( n );
	for ( int i = 0; i < n; i++ ) C[i] = innerProduct( A[i], V );
	return C;
}


//======================================================================


#ifdef SIMD_OPT
vec vectorCombination(float a, vec &U, float b, vec &V) // Linear a*U+b*V
{
	int n = U.size();
	vec W(n), bV(n);
	__m128 m1, m2, m3;
	int nloop = n/4;
	float *p1 = V.data(), *p3 = bV.data(), *p2;

	for(int i = 0; i<nloop; i++){
		//bV[i] = b * V[i];
		m1 = _mm_loadu_ps(p1);
		m2 = _mm_set1_ps(b);
		m3 = _mm_mul_ps(m1, m2);
		_mm_storeu_ps(p3, m3);
		p1 += 4;
		p3 += 4;
	}

	p1 = U.data(), p2 = bV.data(), p3 = W.data();

	for(int i = 0; i<nloop; i++){
		m1 = _mm_load_ps(p1);
		m2 = _mm_load_ps(p2);
		m3 = _mm_add_ps(m1, m2);
		_mm_store_ps(p3, m3);
		p1 += 4;
		p2 += 4;
		p3 += 4;
	}
	return W;
}
#else
vec vectorCombination( float a, const vec &U, float b, const vec &V )        // Linear a*U+b*V combination of vectors
{
	int n = U.size();
	vec W( n );
	for ( int j = 0; j < n; j++ ) W[j] = a * U[j] + b * V[j];
	return W;
}
#endif


//======================================================================


#ifdef SIMD_OPT
float innerProduct( const vec &U, const vec &V )          // Inner product of U and V 内积
{
	// multiply UV
	int size = U.size();
	vec UV(size);
	__m128 m1, m2, m3;
	float *p1 = U.data(), *p2 = V.data(), *p3 = UV.data();
	int nloop = size / 4;
	float sum = 0.0;
	for (int i = 0; i < nloop; ++i)
	{
		m1 = _mm_loadu_ps(p1);
		m2 = _mm_loadu_ps(p2);
		m3 = _mm_mul_ps(m1, m2);
		_mm_storeu_ps(p3, m3);
		p1 += 4;
		p2 += 4;
		p3 += 4;
	}
	// sum UV
	
}

#else

float innerProduct( const vec &U, const vec &V )          // Inner product of U and V 内积
{
	return inner_product( U.begin(), U.end(), V.begin(), 0.0 );
}
#endif

//std库里内积的方法，U的开始到结束，V的开始到结束，这个是自带优化的，

//======================================================================


float vectorNorm( const vec &V )                          // Vector norm
{
	return sqrt( innerProduct( V, V ) );
}


//======================================================================


vec conjugateGradientSolver( const matrix &A, const vec &B )
{
	float TOLERANCE = 1.0e-12;

	int n = A.size();
	vec X( n, 0.0 );

	vec R = B;
	vec P = R;
	int k = 0;

	while ( 1 )
	{
		vec Rold = R;                                         // Store previous residual
		vec AP = matrixTimesVector( A, P );

		float alpha = innerProduct( R, R ) / innerProduct( P, AP );
		X = vectorCombination( 1.0, X, alpha, P );            // Next estimate of solution
		R = vectorCombination( 1.0, R, -alpha, AP );          // Residual

		if ( vectorNorm( R ) < TOLERANCE ) break;             // Convergence test

		float beta = innerProduct( R, R ) / innerProduct( Rold, Rold );
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

