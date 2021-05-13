#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

SEXP permutationDistanceLongestCommonSubsequence(SEXP Rs, SEXP Rt){
	int i,n,*s,*t,*v0,*v1,*vtmp,j,*x,*y;
	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs)+1;
	
	v0 = (int *) R_alloc(sizeof(int), n);
	v1 = (int *) R_alloc(sizeof(int), n);
	
	for(i=0;i<n;i++){
		v0[i]=0;
		v1[i]=0;
	}

	for (y = t, i = 1; i < n; y++, i++){
		for (x = s, j = 1; j < n; x++, j++){
			if(*x == *y){
				v1[j]= v0[j-1]+1;
			}else{				
				v1[j]= ((v0[j] > v1[j-1]) ? v0[j] : v1[j-1]);
			}			
		}
		vtmp = v1;
		v1 = v0;
		v0 = vtmp;
	}
	
	INTEGER(Rval)[0] = v0[n-1];
	UNPROTECT(1);
	return Rval;
}


SEXP permutationDistanceLevenshtein(SEXP Rs, SEXP Rt){
	int *v0;
	int *v1;
	int *vtmp;
	int i,j,ins,del,sub;
	int *x;
	int *y;
	
	int *s,*t,n;	
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs)+1;
	
	v0 = (int *) R_alloc(sizeof(int), n);
	v1 = (int *) R_alloc(sizeof(int), n);

	for (i = 0; i < n; i++)
		v0[i] = i;

	for (y = t, i = 1; i < n; y++, i++){
		v1[0] = i;
		for (x = s, j = 1; j < n; x++, j++){
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
	
	INTEGER(Rval)[0] = v0[n - 1];
	UNPROTECT(1);
	return Rval;
}

SEXP permutationDistanceLongestCommonSubstring(SEXP Rs, SEXP Rt){
	int *v0;
	int *v1;
	int *vtmp;
	int i,j;
	int *x;
	int *y;
	int max=0;
	
	int *s,*t,n;	
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);
	
	v0 = (int *) R_alloc(sizeof(int), n);
	v1 = (int *) R_alloc(sizeof(int), n);

	for(i=0;i<n;i++){
		v0[i]=0;
		v1[i]=0;
	}	
	
	for (y = t, i = 0; i < n; y++, i++){
		for (x = s, j = 0; j < n; x++, j++){
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
	
	INTEGER(Rval)[0] = max;
	UNPROTECT(1);
	return Rval;
}

SEXP permutationDistanceSwapInv(SEXP Rs, SEXP Rt){
	int i,j;
  int tau=0;
	
	int *s,*t,n;	
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);
	
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if((s[i]<s[j]) && (t[i]>t[j]))
				tau++;			
		}
	}

	INTEGER(Rval)[0] = tau;
	UNPROTECT(1);
	return Rval;
}

SEXP permutationDistanceSwap(SEXP Rs, SEXP Rt){
	int i,j;
  int tau=0;
	
	int *s,*t,n;	
	
 	SEXP Rval;
 	SEXP Rvals;
 	SEXP Rvalt;
	
	PROTECT(Rval=allocVector(INTSXP,1));
	
	PROTECT(Rvals=Rf_lang1(Rs));
	PROTECT(Rvalt=Rf_lang1(Rt));	
	
	//s = INTEGER(Rs);
	//t = INTEGER(Rt);
	n = length(Rs);
	
	s = (int *) R_alloc(sizeof(int), n);
	t = (int *) R_alloc(sizeof(int), n);
	
	R_orderVector(s,n,Rvals,TRUE,FALSE);
	R_orderVector(t,n,Rvalt,TRUE,FALSE);
		
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if((s[i]<s[j]) && (t[i]>t[j]))
				tau++;			
		}
	}

	INTEGER(Rval)[0] = tau;
	UNPROTECT(3);
	return Rval;
}


SEXP permutationDistanceInterchange(SEXP Rs, SEXP Rt){
	int i,j;
	int res=0;
	int *piunchecked;
	
	int *s,*s2,*t,n;
 	SEXP Rval;
 	SEXP Rval2;
	PROTECT(Rval=allocVector(INTSXP,1));
	PROTECT(Rval2=Rf_lang1(Rs));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);

	// copy s
	s2= (int *) R_alloc(sizeof(int), n);
	for(i=0;i<n;i++){
		s2[i]= s[i];
	}
	
	//x<-y[order(x)]
	int *idx;
	idx= (int *) R_alloc(sizeof(int), n);
	R_orderVector(idx,n,Rval2,TRUE,FALSE);
	//R_orderVector(idx,n,Rf_lang1(Rs),TRUE,FALSE);
	for(i=0;i<n;i++){
		s2[i]= t[idx[i]];
	}
	
	piunchecked = (int *) R_alloc(sizeof(int), n);
	for(i=0;i<n;i++){
		piunchecked[i]=1;
	}	

	for (i = 0; i < n; i++){
		if(piunchecked[i]==1){
			res++;
			piunchecked[i]=0;
			j = s2[i]-1;
			while(j!=i){
				piunchecked[j] = 0;
				j = s2[j]-1;
			}	
		}		
	}
	
	INTEGER(Rval)[0] = res;
	UNPROTECT(2);
	return Rval;
}

SEXP permutationDistanceR(SEXP Rs, SEXP Rt){
	int i,j;
	int dis=0;

	int *s,*t,n;
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);

	for (i = 0; i < (n-1); i++){
		for (j = 0; j < n; j++){
			if(s[i]==t[j])
				break;			
		}
		if(j<(n-1)){
			dis = dis + ((s[i+1]==t[j+1])? 0 : 1);
		}else{
			dis++;
		}
	}
	INTEGER(Rval)[0] = dis;
	UNPROTECT(1);
	return Rval;
}

SEXP permutationDistanceAdjacency(SEXP Rs, SEXP Rt){
	int i,j;
	int dis=0;
	
	int *s,*t,n;
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);
	
	for (i = 0; i < (n-1); i++){
		for (j = 0; j < n; j++){
			if(j==0){
				if(s[i]==t[j]){
					if(s[i+1]==t[j+1])
						dis++;
				}		
			}else if(j==(n-1)){
				if(s[i]==t[j]){
					if(s[i+1]==t[j-1])
						dis++;
				}		
			}else{
				if(s[i]==t[j]){
					if((s[i+1]==t[j+1]) || (s[i+1]==t[j-1]))
						dis++;
				}
			}	
		}
	}
	INTEGER(Rval)[0] = dis;
	UNPROTECT(1);
	return Rval;
}


SEXP permutationDistancePosition(SEXP Rs, SEXP Rt){
	int i,j;
	int dis=0;
	
	int *s,*t,n;
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);
	
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if(s[i]==t[j]){
				dis=dis+abs(j-i);
			}
		}
	}
	INTEGER(Rval)[0] = dis;
	UNPROTECT(1);
	return Rval;
}

SEXP permutationDistancePosition2(SEXP Rs, SEXP Rt){
	int i,j;
	int dis=0;
	
	int *s,*t,n;
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);
	
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if(s[i]==t[j]){
				dis=dis+((j-i)*(j-i));
			}
		}
	}
	INTEGER(Rval)[0] = dis;
	UNPROTECT(1);
	return Rval;
}

SEXP permutationDistanceHamming(SEXP Rs, SEXP Rt){
	int i;
	int dis=0;
	
	int *s,*t,n;
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);
	
	for (i = 0; i < n; i++){
		if(s[i]!=t[i]){
			dis++;
		}
	}
	
	INTEGER(Rval)[0] = dis;
	UNPROTECT(1);
	return Rval;
}

SEXP permutationDistanceEuclidean(SEXP Rs, SEXP Rt){
	int i;
	int dis=0;
	
	int *s,*t,n;
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);
	
	for (i = 0; i < n; i++){
		dis= dis+(s[i]-t[i]) *(s[i]-t[i]);
	}
	
	INTEGER(Rval)[0] = dis;
	UNPROTECT(1);
	return Rval;
}

SEXP permutationDistanceLee(SEXP Rs, SEXP Rt){
	int i,tmp,ntmp;
	int dis=0;
	
	int *s,*t,n;
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
	s = INTEGER(Rs);
	t = INTEGER(Rt);
	n = length(Rs);
	
	for (i = 0; i < n; i++){
		tmp=abs(s[i]-t[i]);
		ntmp=n-tmp;
		dis=dis+ ((tmp < ntmp)? tmp : ntmp);
	}
	
	INTEGER(Rval)[0] = dis;
	UNPROTECT(1);
	return Rval;
}

// ULAM metric
// Longest increasing subsequence. (2014, December 20). In Wikipedia, The Free Encyclopedia. Retrieved 18:21, December 26, 2014
SEXP permutationDistanceInsert(SEXP Rs, SEXP Rt){
	int i, L, newL, lo, hi, mid;
	
	int *s,*s2,n;
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,1));
 	SEXP Rval2;
	PROTECT(Rval2=Rf_lang1(Rt));
	s = INTEGER(Rs);
	n = length(Rs);

	// copy s
	s2= (int *) R_alloc(sizeof(int), n);
	for(i=0;i<n;i++){
		s2[i]= s[i];
	}
	
	//x<-order(y)[x]
	int *idx;
	idx= (int *) R_alloc(sizeof(int), n);
	R_orderVector(idx,n,Rval2,TRUE,FALSE);
	//R_orderVector(idx,n,Rf_lang1(Rt),TRUE,FALSE);
	for(i=0;i<n;i++){
		s2[i]= idx[s[i]-1]+1;
	}

	// further memory allocation
	int *p = (int *) R_alloc(sizeof(int), n);
	int *m = (int *) R_alloc(sizeof(int), n+1);
	
	L = 0;
	for(i=0; i<n; i++){ 
		lo = 1;
		hi = L;
		while(lo <= hi){
			mid = (lo+hi)/2;
			if(s2[m[mid]] < s2[i]){
				lo = mid+1;
			}else{
				hi = mid-1;
			}
		}
		newL = lo;

		p[i] = m[newL-1];
		m[newL] = i;
		
		if(newL > L){
			L = newL;
		}
	}		
	
	INTEGER(Rval)[0] = L;
	UNPROTECT(2);
	return Rval;
}

SEXP lexPermOrder(SEXP Rs){
	int i,j,sum;
	
	int *s,n,*val;
	s = INTEGER(Rs);
	n = length(Rs);
 	SEXP Rval;
	PROTECT(Rval=allocVector(INTSXP,n));
  val = INTEGER(Rval);
	
	
	for(i=0; i<n; i++){ 
		val[i] = s[i]-1;
		if(i>0){
			sum=0;
			for(j=0; j<i; j++){ 
				if(s[j]<s[i])
					sum++;
			}
			val[i] = val[i] -sum;
		}
	}
	
	UNPROTECT(1);
	return Rval;
}
