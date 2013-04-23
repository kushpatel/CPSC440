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

/* free the pointers of a backup matrix */
void destroyBackupMatrix(Matrix mat)
{
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

/* backs up the contents of a given matrix by pointers without
	making an explicit copy of the matrix */
Matrix backupMatrixPointers(Matrix mat)
{
	Matrix newMat = malloc(sizeof(struct matrix));
	newMat->dim = mat->dim;
	double **m = malloc(mat->dim*sizeof(double *));
	int i, j;
	for(i = 0; i < mat->dim; i++)
	{
		m[i] = mat->m[i];
	}
	newMat->m = m;
	return newMat;
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

/* multiplies matrices : a*b = prod 
	runs in O(n^3) therefore NOT used in QR algorithm*/
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

/* returns a submatrix of A constructed by deleting 1st row and 1st column of A
	assumes that the matrix to be copied is a tridiagonal */
/*Matrix createSubMatrix(Matrix mat)
{
	Matrix sub = initMatrix(mat->dim - 1);
	int i,j;
	for(i = 1; i < mat->dim; i++)
	{
		//only copy elements from the three diagonals
		for(j = i-1; j <= i+1; j++)
		{
			if(j < 1 || j >= mat->dim) continue;
			sub->m[i-1][j-1] = mat->m[i][j];
		}
	}
	return sub;
}*/

/* returns a submatrix of mat constructed by deleting 1st row and 1st column of mat
	reassigns pointer to do this! mat has to be backed up by caller to be freed later!!!
	*/
void createSubMatrix(Matrix mat)
{
	int i;
	for(i = 1; i < mat->dim; i++)
	{
		mat->m[i-1] = mat->m[i] + 1;
	}
	mat->dim--;
}

/* makes all the entries of a matrix zeros --- O(n^2) */
void nullifyMatrix(Matrix mat)
{
	int i, j;
	for(i = 0; i < mat->dim; i++)
	{
		for(j = 0; j < mat->dim; j++)
			mat->m[i][j] = 0;
	}
}

/* multiplies A and B that result into a tridiagonal matrix
	stores the result in A 
	efficient multiplication in O(n^2) */
void tridiagMatMul(Matrix a, Matrix b)
{
	//create copy of A in O(n^2)
	Matrix acopy = copyMatrix(a);
	//make all the entries in A zero in O(n^2)
	nullifyMatrix(a);
	int i, j;
	for(i = 0; i < a->dim; i++)
	{
		//compute elements on diagonal and sub-diagonal of A
		for(j = i-1; j <= i; j++)
		{
			if(j < 0) continue;
			a->m[i][j] = rowByCol(i,j,acopy,b);
		}
	}
	//elements on super-diagonal of A = elements on sub-diagonal of A
	for(i = 1; i < a->dim; i++)
		a->m[i-1][i] = a->m[i][i-1];

	destroyMatrix(acopy);
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
void givensMatrix(int i, int j, int h, Matrix a, double *g)
{
	double x, y, r, c, s;
	// y is the element in A that goes to 0
	y = a->m[i][j];
	// x is the element in A above/below y
	x = a->m[h][j];
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

	x = int h (param)
	y = int i (param)

	Givens rotation matrix g:
	g(k,k) = 1 (for k != i,j)
	g(i,i) = c
	g(j,j) = c
	g(j,i) = -s (for i > j)
	g(i,j) = s (for i > j)
	g(a,b) = 0 (for all other values)

	runs in O(n)
*/
void givensRotate(double *g, int i, int h, Matrix a, int order, int dir)
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
			ik = c*a->m[h][k] + s*a->m[i][k];
			jk = -s*a->m[h][k] + c*a->m[i][k];

			a->m[h][k] = ik;
			a->m[i][k] = jk;
		}
	}
	//post-multiply matrix by Givens rotation matrix
	//only ith and (i-1)th COLUMNS of matrix affected
	else if(order == MAT_BY_GIVEN)
	{
		for(k = 0; k < a->dim; k++)
		{
			ik = c*a->m[k][h] - s*a->m[k][i];
			jk = s*a->m[k][h] + c*a->m[k][i];

			a->m[k][h] = ik;
			a->m[k][i] = jk;
		}
	}
}



/******************************************************************************************/
//               FUNCTIONS FOR ZEROING ENTRIES IN A MATRIX BY GIVENS ROTATIONS
/******************************************************************************************/

/* modified upperhes procedure that takes aMat and turns it into
	the upper hessenberg form in place
	@param aMat = A that goes to its upper hessenberg form in place
	@param n = matrix dimension
	@param id = starts off with I and stores the composition of Givens rotations (U)
	runs in O(n^3)
*/
void upperhes(Matrix aMat, int n, Matrix id)
{
	double *g = malloc(2*sizeof(double));

	int i,j;
	//iterate over all elements to be zeroes (below the lower subdiagonal)
	for(j = 0; j < n - 2; j++)
	{
		for(i = n - 1; i >= j + 2; i--)
		{
			//create the Givens matrix, g(i-1,i) s.t. A[i][j] = 0
			givensMatrix(i,j,i-1,aMat,g);
			if(g[0] != 0 || g[1] != 0)
			{
				//Accumulate G = g(i-1,i) X G...composition of G's
				givensRotate(g,i,i-1,id,GIVEN_BY_MAT,CLOCKWISE);
				//Update A = (g(i-1,i))^T X A...rotate A clockwise
				givensRotate(g,i,i-1,aMat,GIVEN_BY_MAT,CLOCKWISE);
				//Update A = A x g(i-1,i)
				givensRotate(g,i,i-1,aMat,MAT_BY_GIVEN,COUNTER_CLOCKWISE);
			}
		}
	}

	free(g);
}

/* reduces a tridiagonal matrix (aMat) to lower triangular matrix in place
	by givens rotations
	@param aMat = A that is tridiagonal and goes to lower triangular in place
	@param n = dimension of A
	@param id = identity matrix that will accumulate givens transformations applied to A
	runs in O(n^2)
*/
void lowerTriangular(Matrix aMat, int n, Matrix id)
{
	double *g = malloc(2*sizeof(double));

	int i,j;
	//iterate over all elements to be zeroes (on the upper sub-diagonal)
	//because A is assumed to be tridiagonal
	for(i = n - 2 ; i >= 0; i--)
	{
		int j = i+1;
		//create the Givens matrix, g(i-1,i) s.t. A[i][j] = 0
		givensMatrix(i,j,i+1,aMat,g);
		if(g[0] != 0 || g[1] != 0)
		{
			//apply givens rotation to element A(i,j)
			givensRotate(g,i,i+1,aMat,GIVEN_BY_MAT,CLOCKWISE);
			//accumulate the Givens transformation in I
			givensRotate(g,i,i+1,id,GIVEN_BY_MAT,CLOCKWISE);
		}
	}

	free(g);
}

/******************************************************************************************/
//               FUNCTIONS FOR FINDING EIGENVALUES OF A SYMMETRIC MATRIX
/******************************************************************************************/

#define EPSILON (0.00000001)

/* finds the first eigenvalue of B (stored at B(0,0) at the end) in O(n^2) */
void findEigenvalue(Matrix bMat)
{
	int n = bMat->dim;

	//assuming that the #iterations << n => O(n^2)
	while(bMat->m[0][1] > EPSILON || bMat->m[0][1] < -EPSILON)
	{
		//step 2: B = B - uI (where u = B(0,0)) ---- O(n)
		double mu = bMat->m[0][0];
		int i;
		for(i = 0; i < n; i++)
			bMat->m[i][i] -= mu;

		//step 3: Bring B to lower triangular form ---- O(n^2)
		Matrix qMat = idMat(n);
		lowerTriangular(bMat,n,qMat);

		//step 4: Apply adjoints of Givens rotation to L ---- each O(n^2)
		transposeMatrix(qMat);
		tridiagMatMul(bMat,qMat);

		//step 5: Add the shift mu back to the diagonal ---- O(n)
		for(i = 0; i < n; i++)
			bMat->m[i][i] += mu;

		destroyMatrix(qMat);
	}		
}

/* computes a diagonal matrix containing eigenvalues of A
	runs in O(n^3) */
void qr_symmetric(double * a, int n, double * b)
{
	//given symmetric matrix A ---- each O(n^2)
	Matrix aMat = createMatrix(a,n);
	Matrix uMat = idMat(n);
	Matrix eigenVals = idMat(n);

	//step 1: feed into upperhes(A,n,U,B)...B will be symmetric tridiagonal ---- O(n^3)
	//upperhes was modified so that A (aMat) is modified to B
	upperhes(aMat,n,uMat);
	//so that aMat can be freed later!
	Matrix bMat = backupMatrixPointers(aMat);

	int i;
	//each iteration is O(n^2) to give total run-time O(n^3)
	for(i = 0; i < n - 1; i++)
	{
		//find each eigenvalue in O(n^2)
		findEigenvalue(bMat);
		eigenVals->m[i][i] = bMat->m[0][0];
		//pointer magic to do this in O(n)
		createSubMatrix(bMat);
	}
	eigenVals->m[i][i] = bMat->m[0][0];
		
	//fill the diagonal matrix in b as output
	flattenMatrix(eigenVals,b);
	
	//exit after freeing malloc'ed memory
	destroyMatrix(aMat);
	destroyBackupMatrix(bMat);
	destroyMatrix(uMat);
	destroyMatrix(eigenVals);
}

int main(int argc, char *argv[])
{
	//set dimensions of matrix here
	int n = 10;
	//double a[9] = {7,5,4,5,8,3,4,3,9};

	//create a test matrix
	Matrix aMat = initMatrix(n);

	//make A a symmetric matrix
	int i, j;
	for(i = 0; i < n; i++)
	{
		for(j = i; j < aMat->dim; j++)
		{
			double r = ((double)rand()/(double)RAND_MAX);
			aMat->m[i][j] = r;
			aMat->m[j][i] = r;
		}
	}
	
	//allocate space for input matrix A and output matrix B
	double *a = malloc(n*n*sizeof(double));
	double *b = malloc(n*n*sizeof(double));

	flattenMatrix(aMat,a);

	//run the main process for Hessenberg decomposition
	qr_symmetric(a,n,b);

	//create matrix structures from output for testing
	Matrix bMat = createMatrix(b,n);

	puts("A = ");
	printMatrix(aMat);
	puts("B = ");
	printMatrix(bMat);

	puts("Eigenvalues = ");
	for(i = 0; i < n; i++)
		printf("%f \n", bMat->m[i][i]);

	//free malloc'ed memory before exiting
	free(a);
	free(b);
	destroyMatrix(aMat);
	destroyMatrix(bMat);

	return 0;
}