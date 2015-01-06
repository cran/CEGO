#include <R.h>
#include <stdio.h>
#include <stdlib.h>

void permutationDistanceLongestCommonSubsequence(int *s, int *t, int *n, int *ret){
	int *v0;
	int *v1;
	int *vtmp;
	int i,j;
	int *x;
	int *y;
    
	v0 = (int *) R_alloc(sizeof(int), *n +1);
	v1 = (int *) R_alloc(sizeof(int), *n +1);
	
	for(i=0;i<(*n+1);i++){
		v0[i]=0;
		v1[i]=0;
	}

	for (y = t, i = 1; i < (*n+1); y++, i++){
		for (x = s, j = 1; j < (*n+1); x++, j++){
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
	ret[0]= v0[*n];
}

void permutationDistanceLevenshtein(int *s, int *t, int *n, int *ret){
	int *v0;
	int *v1;
	int *vtmp;
	int i,j,ins,del,sub;
	int *x;
	int *y;
    
	v0 = (int *) R_alloc(sizeof(int), *n);
	v1 = (int *) R_alloc(sizeof(int), *n);

	for (i = 0; i < *n; i++)
		v0[i] = i;

	for (y = t, i = 1; i < *n; y++, i++){
		v1[0] = i;
		for (x = s, j = 1; j < *n; x++, j++){
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
	ret[0]= v0[*n - 1];
}

void permutationDistanceLongestCommonSubstring(int *s, int *t, int *n, int *ret){
	int *v0;
	int *v1;
	int *vtmp;
	int i,j;
	int *x;
	int *y;
    int max=0;
	
	v0 = (int *) R_alloc(sizeof(int), *n);
	v1 = (int *) R_alloc(sizeof(int), *n);

	for(i=0;i<*n;i++){
		v0[i]=0;
		v1[i]=0;
	}	
	
	for (y = t, i = 0; i < *n; y++, i++){
		for (x = s, j = 0; j < *n; x++, j++){
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
	ret[0]= max;
}

void permutationDistanceSwap(int *s, int *t, int *n, int *ret){
	int i,j;
    int tau=0;
	
	for (i = 0; i < *n; i++){
		for (j = 0; j < *n; j++){
			if((s[i]<s[j]) && (t[i]>t[j]))
				tau++;			
		}
	}
	ret[0]= tau;
}




void permutationDistanceInterchange(int *s, int *n, int *ret){
	int i,j;
    int cc=0;
	int *piunchecked;
	piunchecked = (int *) R_alloc(sizeof(int), *n);
	for(i=0;i<*n;i++){
		piunchecked[i]=1;
	}	

	for (i = 0; i < *n; i++){
		if(piunchecked[i]==1){
			cc++;
			piunchecked[i]=0;
			j = s[i]-1;
			while(j!=i){
				piunchecked[j] = 0;
				j = s[j]-1;
			}	
		}		
	}
	ret[0]= cc;
}




void permutationDistanceR(int *s, int *t, int *n, int *ret){
	int i,j;
	int dis=0;
	
	for (i = 0; i < (*n-1); i++){
		for (j = 0; j < *n; j++){
			if(s[i]==t[j])
				break;			
		}
		if(j<(*n-1)){
			dis = dis + ((s[i+1]==t[j+1])? 0 : 1);
		}else{
			dis++;
		}
	}
	ret[0]= dis;
}

void permutationDistanceAdjacency(int *s, int *t, int *n, int *ret){
	int i,j;
	int dis=0;
	
	for (i = 0; i < (*n-1); i++){
		for (j = 0; j < *n; j++){
			if(j==0){
				if(s[i]==t[j]){
					if(s[i+1]==t[j+1])
						dis++;
				}		
			}else if(j==(*n-1)){
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
	ret[0]= dis;
}

void permutationDistancePosition(int *s, int *t, int *n, int *ret){
	int i,j;
	int dis=0;
	
	for (i = 0; i < *n; i++){
		for (j = 0; j < *n; j++){
			if(s[i]==t[j]){
				dis=dis+abs(j-i);
			}
		}
	}
	ret[0]= dis;
}

void permutationDistancePosition2(int *s, int *t, int *n, int *ret){
	int i,j;
	int dis=0;
	
	for (i = 0; i < *n; i++){
		for (j = 0; j < *n; j++){
			if(s[i]==t[j]){
				dis=dis+((j-i)*(j-i));
			}
		}
	}
	ret[0]= dis;
}

void permutationDistanceHamming(int *s, int *t, int *n, int *ret){
	int i;
	int dis=0;
	
	for (i = 0; i < *n; i++){
		if(s[i]!=t[i]){
			dis++;
		}
	}
	ret[0]= dis;
}

void permutationDistanceEuclidean(int *s, int *t, int *n, int *ret){
	int i;
	int dis=0;
	
	for (i = 0; i < *n; i++){
		dis= dis+(s[i]-t[i]) *(s[i]-t[i]);
	}
	ret[0]= dis;
}

void permutationDistanceLee(int *s, int *t, int *n, int *ret){
	int i,tmp,ntmp;
	int dis=0;
	
	for (i = 0; i < *n; i++){
		tmp=abs(s[i]-t[i]);
		ntmp=*n-tmp;
		dis=dis+ ((tmp < ntmp)? tmp : ntmp);
	}
	ret[0]= dis;
}

//not a distance, but used for ULAM metric
//based on:
// Longest increasing subsequence. (2014, December 20). In Wikipedia, The Free Encyclopedia. Retrieved 18:21, December 26, 2014, from http://en.wikipedia.org/w/index.php?title=Longest_increasing_subsequence&oldid=638943901
void longestIncreasingSubsequence(int *s, int *n, int *ret) {
	int i, L, newL, lo, hi, mid;
	int *p = (int *) R_alloc(sizeof(int), *n);
	int *m = (int *) R_alloc(sizeof(int), *n+1);
	
	L = 0;
	for(i=0; i<*n; i++){ 
		lo = 1;
		hi = L;
		while(lo <= hi){
			mid = (lo+hi)/2;
			if(s[m[mid]] < s[i]){
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
	ret[0]=L;
}

void lexPermOrder(int *s, int *n, int *ret){
	int i,j,sum;

	for(i=0; i<*n; i++){ 
		ret[i] = s[i]-1;
		if(i>0){
			sum=0;
			for(j=0; j<i; j++){ 
				if(s[j]<s[i])
					sum++;
			}
			ret[i] = ret[i] -sum;
		}
	}
	
}
