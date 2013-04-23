#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/******************************************************************************************/
//           2 BY 2 MATRIX WITH ITS METHODS 
/******************************************************************************************/

typedef struct matrix *Matrix;

/* 2D matrix structure */
struct matrix
{
	int dim;
	double **m;
};

/* initialize a null nxn matrix */
Matrix initMatrix(int n)
{
	Matrix newMat = malloc(sizeof(struct matrix));
	newMat->dim = n;
	double **mat = malloc(n*sizeof(double *));
	int i, j;
	for(i = 0; i < n; i++)
	{
		mat[i] = malloc(n*sizeof(double));
		for(j = 0; j < n; j++)
			mat[i][j] = 0;
	}
	newMat->m = mat;
	return newMat;
}

/* free a given matrix */
void destroyMatrix(Matrix mat)
{
	int i,j;
	for(i = 0; i < mat->dim; i++)
	{
		free(mat->m[i]);
	}
	free(mat->m);
	free(mat);
}

/* print out the given matrix in 2D format*/
void printMatrix(Matrix mat)
{
	int i,j;
	for(i = 0; i < mat->dim; i++)
	{
		for(j = 0; j < mat->dim; j++)
		{
			//printf("%g ", mat->m[i][j]);
			printf("%f ", mat->m[i][j]);
		}
		printf("\n");
	}
	puts("");
}

/* create a clone of matrix provided */
Matrix copyMatrix(Matrix mat)
{
	Matrix copy = initMatrix(mat->dim);
	copy->dim = mat->dim;
	int i, j;
	for(i = 0; i < mat->dim; i++)
	{
		for(j = 0; j < mat->dim; j++)
		{
			copy->m[i][j] = mat->m[i][j];
		}
	}
	return copy;
}

/* create a 2D matrix from a 1D array */
Matrix createMatrix(double *mtx, int n)
{
	Matrix mat = initMatrix(n);
	int i,j;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++, mtx++)
			mat->m[i][j] = *mtx;
	}
	return mat;
}

/* transpose a given matrix in place */
void transposeMatrix(Matrix mat)
{
	int i, j;
	for(i = 0; i < mat->dim; i++)
	{
		for(j = i; j < mat->dim; j++)
		{
			double temp = mat->m[i][j];
			mat->m[i][j] = mat->m[j][i];
			mat->m[j][i] = temp;
		}
	}
}

/* flatten a 2D matrix back to a 1D double* */
void flattenMatrix(Matrix mat, double *m)
{
	int i, j;
	for(i = 0; i < mat->dim; i++)
	{
		for(j = 0; j < mat->dim; j++, m++)
		{
			*m = mat->m[i][j];
		}
	}
}

/* helper function for matrix multiplication
	returns row by column sum of specified row and col
	in matrices a and b respectively */
double rowByCol(int row, int col, Matrix a, Matrix b)
{
	double d = 0;
	int i;
	for(i = 0; i < a->dim; i++)
	{
		d += a->m[row][i] * b->m[i][col];
	}
	return d;
}

/* multiplies matrices : a*b = prod */
void matMul(Matrix a, Matrix b, Matrix prod)
{
	int i,j;
	for(i = 0; i < prod->dim; i++)
	{
		for(j = 0; j < prod->dim; j++)
			prod->m[i][j] = rowByCol(i,j,a,b);
	}
}

/* returns an nxn identity matrix */
Matrix idMat(int n)
{
	Matrix mat = initMatrix(n);
	int i;
	for(i = 0; i < n; i++)
	{
		mat->m[i][i] = 1;
	}
	return mat;
}

/******************************************************************************************/
//               FUNCTIONS FOR GIVENS ROTATIONS
/******************************************************************************************/

#define CLOCKWISE (2)
#define COUNTER_CLOCKWISE (3)
#define GIVEN_BY_MAT (4)
#define MAT_BY_GIVEN (5)

/* 	calculate sqrt(x^2 + y^2) using hypot() function 
	|x|sqrt(1 + (y/x)^2) given that |x| > |y|
	avoids overflows when x or y are too big 
*/
double hypot(double x, double y)
{
	double t;
	x = fabs(x);
	y = fabs(y);

	if(x > y)
	{
		t = y;
	}
	else
	{
		t = x;
		x = y;
	} 
	//t = min(x,y);
	//x = max(x,y);
	t = t/x;
	return x*sqrt(1+t*t);
}

/* returns values of a Givens matrix as follows:
	g[0] = c; g[1] = s;
*/
void givensMatrix(int i, int j, Matrix a, double *g)
{
	double x, y, r, c, s;
	// y is the element in A that goes to 0
	y = a->m[i][j];
	// x is the element in A above y
	x = a->m[i-1][j];
	//r = sqrt(x*x+y*y);
	//use hypot function instead of sqrt to avoid overflow
	r = hypot(x,y);
	// avoid division by zero...nothing to rotate if x = y = 0
	if(r == 0)
	{
		c = 0;
		s = 0;
	}
	else 
	{
		c = x/r;
		s = -y/r;
	}

	g[0] = c;
	g[1] = s;
}

/*	|c  s| |x| = |r|
	|-s c| |y|   |0|
	r = sqrt(x^2 + y^2)
	c = x/r
	s = -y/r

	Givens rotation matrix g:
	g(k,k) = 1 (for k != i,j)
	g(i,i) = c
	g(j,j) = c
	g(j,i) = -s (for i > j)
	g(i,j) = s (for i > j)
	g(a,b) = 0 (for all other values)
*/
void givensRotate(double *g, int i, Matrix a, int order, int dir)
{
	double c = g[0];
	double s = g[1];
	if(dir == CLOCKWISE) s = -s;

	int k;
	double ik, jk, idk, jdk;
	//pre-multiply matrix by Givens rotation matrix
	//only ith and (i-1)th ROWS of matrix affected
	if(order == GIVEN_BY_MAT)
	{
		for(k = 0; k < a->dim; k++)
		{
			ik = c*a->m[i-1][k] + s*a->m[i][k];
			jk = -s*a->m[i-1][k] + c*a->m[i][k];

			a->m[i-1][k] = ik;
			a->m[i][k] = jk;
		}
	}
	//post-multiply matrix by Givens rotation matrix
	//only ith and (i-1)th COLUMNS of matrix affected
	else if(order == MAT_BY_GIVEN)
	{
		for(k = 0; k < a->dim; k++)
		{
			ik = c*a->m[k][i-1] - s*a->m[k][i];
			jk = s*a->m[k][i-1] + c*a->m[k][i];

			a->m[k][i-1] = ik;
			a->m[k][i] = jk;
		}
	}
}

/* brings a given matrix A to upper Hessenberg form and returns it in B
	UAU^T = B */
void upperhes(double *a, int n, double *u, double *b)
{
	//stores composition of Givens transformations applied = G
	Matrix id = idMat(n);
	//create input matrix struct
	Matrix aMat = createMatrix(a,n);
	//stores c and s values for Givens rotation matrix
	double *g = malloc(2*sizeof(double));

	int i,j;
	//iterate over all elements to be zeroes (below the lower subdiagonal)
	for(j = 0; j < n - 2; j++)
	{
		for(i = n - 1; i >= j + 2; i--)
		{
			//create the Givens matrix, g(i-1,i) s.t. A[i][j] = 0
			givensMatrix(i,j,aMat,g);
			if(g[0] != 0 || g[1] != 0)
			{
				//Accumulate G = g(i-1,i) X G...composition of G's
				givensRotate(g,i,id,GIVEN_BY_MAT,CLOCKWISE);
				//Update A = (g(i-1,i))^T X A...rotate A clockwise
				givensRotate(g,i,aMat,GIVEN_BY_MAT,CLOCKWISE);
				//Update A = A x g(i-1,i)
				givensRotate(g,i,aMat,MAT_BY_GIVEN,COUNTER_CLOCKWISE);
			}
		}
	}

	//aMat is now B
	flattenMatrix(aMat,b);
	//id is now U
	flattenMatrix(id,u);

	//exit after freeing malloc'ed memory
	destroyMatrix(id);
	destroyMatrix(aMat);
	free(g);
}

int main(int argc, char *argv[])
{
	//double *a = malloc(dim*dim*sizeof(double));
	//double a[9] = {1,0,0,0,1,0,0,0,1};
	//double a[9] = {6,5,4,5,1,4,5,4,3};
	//double a[16] = {1,2,3,4,5,6,3,8,9,3,4,6,3,7,1,9};
	//double a[25] = {1,2,3,4,5,5,6,3,8,9,9,3,4,6,7,3,7,1,9,1,1,2,3,4,5};
	//double a[9] = {1,2,3,4,5,6,3,8,9};
	//double a[9] = {7.8102,4.4813,2.5607,0,-2.4327,3.0729,0,4,3};
	//double a[9] = {5.4,4,7.7,3.5,-.7,2.8,-3.2,5.1,.8};
	//double a[16] = {1,9,8,7,9,2,6,5,8,6,3,1,7,5,1,4};
	
	//choose dimension of matrix here!
	int n = 10;

	//create a test matrix A using rand()
	double *a = malloc(n*n*sizeof(double));
	double *atemp = a;
	int i;
	for(i = 0; i < n*n; i++, atemp++)
	{
		*atemp = ((double)rand()/(double)RAND_MAX);
	}
	
	//allocate space for output matrices B and U
	double *b = malloc(n*n*sizeof(double));
	double *u = malloc(n*n*sizeof(double));

	//run the main process for Hessenberg decomposition
	upperhes(a,n,u,b);

	//create matrix structures from input and output for testing
	Matrix aMat = createMatrix(a,n);
	Matrix bMat = createMatrix(b,n);
	Matrix uMat = createMatrix(u,n);
	Matrix uTran = createMatrix(u,n);
	transposeMatrix(uTran);
	
	//result2 = UAU^T
	Matrix result = initMatrix(n);
	Matrix result2 = initMatrix(n);
	matMul(uMat,aMat,result);
	matMul(result,uTran,result2);

	puts("A = ");
	printMatrix(aMat);
	puts("B = ");
	printMatrix(bMat);
	puts("UAU^T");
	printMatrix(result2);

	//free malloc'ed memory before exiting
	free(a);
	free(b);
	free(u);
	destroyMatrix(aMat);
	destroyMatrix(bMat);
	destroyMatrix(uMat);
	destroyMatrix(uTran);
	destroyMatrix(result);
	destroyMatrix(result2);

	return 0;
}