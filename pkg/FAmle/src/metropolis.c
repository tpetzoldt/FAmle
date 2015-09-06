#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

int CHOLESKY(double *mat, double *diag, int n, double *det)
{
	*det = 1;
	int i, j, k;
	double sum;
	int err = 0;
	for(i = 0; i < n; i++)
	{
		for(j = i; j < n; j++)
		{
			for(sum = mat[i + j*n],k = i - 1; k >= 0; k--) sum -= mat[i + k*n]*mat[j + k*n];
			if(i == j)
			{
				if(sum <= 0.0)
				{
					err = 1;
					*det = R_NaN;
					return err;
				}
				diag[i] = sqrt(sum);
			}
			else mat[j + i*n] = sum/diag[i];
		}
	}
	for(i = 0; i < n; i++) *det *= diag[i]*diag[i];
	return err;
}

void set_sqrt_sigma(double *sigma, int d, double *sqrt_sigma)
{
	int i, j;
	double mat[d*d], diag[d], det;
	for(i = 0; i < d*d; i++) mat[i] = *(sigma+i);
	CHOLESKY(&mat[0],&diag[0],d,&det);
	for(i = 0; i < d; i++)
	{
		for(j = 0; j < d; j++)
		{
			if(i > j) *(sqrt_sigma+i+j*d) = mat[i + j*d];
			else if(i == j) *(sqrt_sigma+i+j*d) = diag[i];
			else *(sqrt_sigma+i+j*d) = 0;
		}
	}
}

void rmvnorm(double *mu, double *sqrt_sigma, int d, double *out)
{
	int i, j;
	double z[d];
	for(i = 0; i < d; i++) z[i] = rnorm(0,1);
	for(i = 0; i < d; i++)
	{
		*(out+i) = 0;
		for(j = 0; j < d; j++) *(out + i) += *(sqrt_sigma+i+d*j)*z[j];
		*(out + i) += *(mu+i);
	}
}

void set_temp_vec(double *values, double *vec, int d)
{
	int i;
	for(i = 0; i < d; i++) *(vec+i) = *(values+i);
}

SEXP metropolis(SEXP iter, SEXP start, SEXP covariance, SEXP expr, SEXP rho)
{
	int i, d = length(start), n = INTEGER(iter)[0];
	double sqrt_sigma[d*d];
	set_sqrt_sigma(REAL(covariance),d,&sqrt_sigma[0]);
	SEXP a, b, sims, rate, list_out;
	PROTECT(a = allocVector(REALSXP,d));
	PROTECT(b = allocVector(REALSXP,d));
	PROTECT(sims = allocMatrix(REALSXP,d,n));
	PROTECT(rate = allocVector(REALSXP,n));
	PROTECT(list_out = allocVector(VECSXP,2));
	defineVar(install("a"),a,rho);
	defineVar(install("b"),b,rho);
	for(i = 0; i < d; i++) *(REAL(sims)+i) = *(REAL(start)+i);
	*(REAL(rate)+0) = 0;
	double ratio;
	set_temp_vec(REAL(sims),REAL(a),d);
	GetRNGstate();
	for(i = 1; i < n; i++)
	{
		rmvnorm(REAL(a),&sqrt_sigma[0],d,REAL(b));
		ratio = *REAL(eval(expr,rho));
		if(!ISNAN(ratio) && runif(0,1) < ratio)
		{
			set_temp_vec(REAL(b),REAL(a),d);
			*(REAL(rate)+i) = 1;
		}
		else *(REAL(rate)+i) = 0;
		set_temp_vec(REAL(a),&(*(REAL(sims)+i*d)),d);
	}
	PutRNGstate();
	SET_VECTOR_ELT(list_out,0,sims);
	SET_VECTOR_ELT(list_out,1,sims);
	UNPROTECT(5);
	return list_out;
}
