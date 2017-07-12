#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
SEXP stringDistanceHamming(SEXP Rs, SEXP Rt) {
	int i,n1,n2,n,dis=0;
	char *s,*t;  
	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	PROTECT(Rs = AS_CHARACTER(Rs)); 
	PROTECT(Rt = AS_CHARACTER(Rt)); 
	
	// allocate memory:
	s = R_alloc(strlen(CHAR(STRING_ELT(Rs, 0))), sizeof(char)); 
	t = R_alloc(strlen(CHAR(STRING_ELT(Rt, 0))), sizeof(char)); 

	// ... and copy Rs to s: 
	strcpy(s, CHAR(STRING_ELT(Rs, 0))); 
	strcpy(t, CHAR(STRING_ELT(Rt, 0))); 

	// string length and minimum length
	n1 = strlen(CHAR(STRING_ELT(Rs, 0)));
	n2 = strlen(CHAR(STRING_ELT(Rt, 0)));
	if(n1<=n2){
		n = n1;
	}else{
		n = n2;
	}			
	
	for (i = 0; i < n; i++){
		if(s[i]!=t[i]){
			dis++;
		}
	}
	
	//debug
	//Rprintf(" Debug Printed:");
	//Rprintf(" %s %s \n",s,t);
	//Rprintf(" %c \n",s[1]);	 
	//Rprintf(" %d \n",n1);	 
	//Rprintf(" %d \n",n2);	 
	//Rprintf(" %d \n",n);	 

	INTEGER(Rval)[0] = dis+abs(n1-n2);
	UNPROTECT(3);
	return Rval;
}
		
SEXP stringDistanceLevenshtein(SEXP Rs, SEXP Rt) {
	int *v0,*v1,*vtmp;
	int i,j,ins,del,sub;
	int ns,nt;	
	char *x,*y;		
	char *s,*t;  
	
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	PROTECT(Rs = AS_CHARACTER(Rs)); 
	PROTECT(Rt = AS_CHARACTER(Rt)); 
	
	// allocate memory:
	s = R_alloc(strlen(CHAR(STRING_ELT(Rs, 0))), sizeof(char)); 
	t = R_alloc(strlen(CHAR(STRING_ELT(Rt, 0))), sizeof(char)); 

	// ... and copy Rs to s: 
	strcpy(s, CHAR(STRING_ELT(Rs, 0))); 
	strcpy(t, CHAR(STRING_ELT(Rt, 0))); 
	
	// string lengths
	ns = strlen(CHAR(STRING_ELT(Rs, 0)))+1;
	nt = strlen(CHAR(STRING_ELT(Rt, 0)))+1;
	
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
	UNPROTECT(3);
	return Rval;
}
		
SEXP stringDistanceLongestCommonSubstring(SEXP Rs, SEXP Rt) {
	int *v0,*v1,*vtmp;
	char *x,*y;		
	char *s,*t;  
	int i,j,ns,nt,max=0;	
		
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	PROTECT(Rs = AS_CHARACTER(Rs)); 
	PROTECT(Rt = AS_CHARACTER(Rt)); 
	
	// allocate memory:
	s = R_alloc(strlen(CHAR(STRING_ELT(Rs, 0))), sizeof(char)); 
	t = R_alloc(strlen(CHAR(STRING_ELT(Rt, 0))), sizeof(char)); 

	// ... and copy Rs to s: 
	strcpy(s, CHAR(STRING_ELT(Rs, 0))); 
	strcpy(t, CHAR(STRING_ELT(Rt, 0))); 
	
	// string lengths
	ns = strlen(CHAR(STRING_ELT(Rs, 0)));
	nt = strlen(CHAR(STRING_ELT(Rt, 0)));
	
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
	UNPROTECT(3);
	return Rval;
}
