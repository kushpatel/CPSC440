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
			printf("%f ", mat->m[i][j]);
		}
		printf("\n");
	}
	puts("");
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
//               VECTOR MANIPULATION FUNCTIONS FOR GRAM-SCHMIDT
/******************************************************************************************/

/* prints out contents of a double * vector
	handy in debugging */
void printVector(double *vector, int dim)
{
	int i;
	for(i = 0; i < dim; i++, vector++)
	{
		printf("%f ",*vector);
	}
	printf("\n");
}

/* returns the inner product of vector a and b */
double dotProd(double *a, double *b, int dim)
{
	double prod = 0;
	int i;
	for(i = 0; i < dim; i++)
	{
		prod += a[i] * b[i];
	}
	return prod;
}

/* returns the magnitude of the given vector */
double magnitude(double *vector, int dim)
{
	double mag = 0;
	int i;
	for(i = 0; i < dim; i++)
	{
		mag += vector[i] * vector[i];
	}
	return sqrt(mag);
}

/* returns the projection vector of u on v inside v
	same transformation is applied to v_id */
void proj(double *u, double *v, double *u_id, double *v_id, int dim)
{
	double u_v = dotProd(u,v,dim);
	double u_u = dotProd(u,u,dim);
	int i;
	for(i = 0; i < dim; i++)
	{
		v[i] -= (u[i] * u_v) / u_u;
		v_id[i] -= (u_id[i] * u_v) / u_u;
	}
}

/* normalizes the given vector v in place 
	same transformation is applied to v_id*/
void norm(double *v, double *v_id, int dim)
{
	double mag = magnitude(v,dim);
	int i;
	for(i = 0; i < dim; i++)
	{
		v[i] /= mag;
		v_id[i] /= mag;
	}
}

/* returns the row number of the vector with max magnitude */
int maxMagVector(Matrix mat, int pos)
{
	double maxMag = magnitude(mat->m[pos],mat->dim);
	int maxPos = pos;
	int i;
	for(i = pos + 1; i < mat->dim; i++)
	{
		double mag = magnitude(mat->m[i],mat->dim);
		if(mag > maxMag)
		{
			maxMag = mag;
			pos = i;
		}
	}
	return maxPos;
}

/******************************************************************************************/
//               GRAM-SCHMIDT METHOD FOR INVERTING A MATRIX
/******************************************************************************************/

/* orthogonalizes both a and id s.t A -> U and I -> G
	s.t AG = U */
void gramSchmidt(Matrix a, Matrix id)
{
	int i, j, k;
	for(i = 0; i < a->dim; i++)
	{
		//get the vector with maximum magnitude and switch it to be at pos = i
		int maxMag = maxMagVector(a,i);
		double *temp_a = a->m[i];
		double *temp_id = id->m[i];
		a->m[i] = a->m[maxMag];
		a->m[maxMag] = temp_a;
		id->m[i] = id->m[maxMag];
		id->m[maxMag] = temp_id;


		//orthogonalize vector at i with all the vectors preceding it
		for(k = 0; k < i; k++)
		{
			proj(a->m[k],a->m[i],id->m[k],id->m[i],a->dim);
		}

		//normalize the vector at i
		norm(a->m[i],id->m[i],a->dim);

		//orthogonalize all the vectors after i to i
		for(j = i+1; j < a->dim; j++)
		{
			proj(a->m[i],a->m[j],id->m[i],id->m[j],a->dim);
		}
	}
}

/* main method that inverts an nxn matrix A and returns U and A inverse in its parameters
	double *a = square invertible matrix A that is to be inverted
	int n = dimension of A
	double *u = the orthogonalized matrix U returned in this
	double *b = the inverse of A returned in this
*/
void inv_double_gs(double *a, int n, double *u, double *b)
{
	Matrix aMat = createMatrix(a,n);
	// an nxn identity matrix
	Matrix id = idMat(n);
	// transpose aMat so as to run GS on its columns
	transposeMatrix(aMat);
	// run GS ... orthogonalizes aMat by columns
	gramSchmidt(aMat,id);	// aMat is now U transpose
	// get G matrix
	transposeMatrix(id);
	//printMatrix(aMat);		//print U transpose
	//printMatrix(id);			//print G

	// A inverse = G * U transpose
	Matrix ainv = initMatrix(n);
	matMul(id,aMat,ainv);
	//printMatrix(ainv);		//print inverse of A

	//return U in double *u and A inverse in double *b params
	transposeMatrix(aMat);
	flattenMatrix(aMat,u);
	flattenMatrix(ainv,b);

	// clean up all matrices
	destroyMatrix(aMat);
	destroyMatrix(id);
	destroyMatrix(ainv);
}

int main(int argc, char *argv[])
{
	//choose dimension here
	int dim = 20;

	//allocate memory for A, U and B
	double *a = malloc(dim*dim*sizeof(double));
	double *u = malloc(dim*dim*sizeof(double));
	double *b = malloc(dim*dim*sizeof(double));
	
	//create a test matrix A using rand()
	double *atemp = a;
	int i;
	for(i = 0; i < dim*dim; i++, atemp++)
	{
		*atemp = ((double)rand()/(double)RAND_MAX);
	}

	//invert matrix by double GS procedure
	inv_double_gs(a,dim,u,b);

	//print out matrices A, U and B (A inverse)
	double *a2 = a;
	double *u2 = u;
	double *b2 = b;
	for(i = 0; i < dim*dim; i++, a2++)
	{
		if(!(i % dim))
			puts("");
		printf("%f ",*a2);
	}
	puts("");

	for(i = 0; i < dim*dim; i++, u2++)
	{
		if(!(i % dim))
			puts("");
		printf("%f ",*u2);
	}

	puts("");
	for(i = 0; i < dim*dim; i++, b2++)
	{
		if(!(i % dim))
			puts("");
		printf("%f ",*b2);
	}
	puts("");

	// check if A*B = I
	Matrix amat = createMatrix(a,dim);
	Matrix bmat = createMatrix(b,dim);
	Matrix cmat = initMatrix(dim);

	matMul(amat,bmat,cmat);
	puts("");

	printMatrix(cmat);

	destroyMatrix(amat);
	destroyMatrix(bmat);
	destroyMatrix(cmat);

	//free all memory
	free(a);
	free(u);
	free(b);

	return 0;
}