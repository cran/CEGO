#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
SEXP numericDistanceLevenshtein(SEXP Rs, SEXP Rt) {
	int *v0,*v1,*vtmp;
	int i,j,ins,del,sub;
	int ns,nt;	
	double *x,*y,*s,*t;  
	
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = REAL(Rs);
	t = REAL(Rt);

	// vector lengths
	ns = length(Rs)+1;
	nt = length(Rt)+1;
	
	v0 = (int *) R_alloc(sizeof(int), ns);
	v1 = (int *) R_alloc(sizeof(int), ns);

	for (i = 0; i < (ns); i++)
		v0[i] = i;

	for (y = t, i = 1; i < nt; y++, i++){
		v1[0] = i;
		for (x = s, j = 1; j < ns; x++, j++){
			ins = v0[j] + 1; 
			del = v1[j - 1] + 1;
			sub = v0[j - 1] + ((*x == *y) ? 0 : 1);

			v1[j] = ((ins < del) ? ins : del);
			v1[j] = ((v1[j] < sub) ? v1[j] : sub);
		}
		vtmp = v1;
		v1 = v0;
		v0 = vtmp;
	}
	
	INTEGER(Rval)[0] = v0[ns-1];
	UNPROTECT(1);
	return Rval;
}
		
SEXP numericDistanceLongestCommonSubstring(SEXP Rs, SEXP Rt) {
	int *v0,*v1,*vtmp;
	double *x,*y,*s,*t;  
	int i,j,ns,nt,max=0;	
		
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = REAL(Rs);
	t = REAL(Rt);

	// vector lengths
	ns = length(Rs);
	nt = length(Rt);
	
	v0 = (int *) R_alloc(sizeof(int), ns);
	v1 = (int *) R_alloc(sizeof(int), ns);

	for(i=0;i<ns;i++){
		v0[i]=0;
		v1[i]=0;
	}	
	
	for (y = t, i = 0; i < nt; y++, i++){
		for (x = s, j = 0; j < ns; x++, j++){
			if(*x == *y){
				if((j==0) || (i==0)){
					v1[j]=1;
				}else{
					v1[j]=1+v0[j-1];
				}
				if(v1[j]>max){
					max = v1[j];	
				}
			}else{
				v1[j]=0;
			}			
		}
		vtmp = v1;
		v1 = v0;
		v0 = vtmp;
	}
	
	max = ((ns > nt) ? ns-max : nt - max);
	INTEGER(Rval)[0] = max;
	UNPROTECT(1);
	return Rval;
}
