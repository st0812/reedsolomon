#include "RS.h"
#include "GF.h"
#include "Poly.h"
#include "Matrix.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

static unsigned int s_N;
static unsigned int s_K;
static unsigned int s_t;

static Poly calcI(unsigned int* data, unsigned int len);
static Poly calcG();
static Poly calcGFactor(unsigned int i);
static void calcSyndrome(Poly Y, unsigned int * syndrome);
static int calcNumberOfErrors(unsigned int* syndrome);
static int calcPositionOfErrors(unsigned int * syndrome, int number_of_errors, unsigned int* positions);
static int calcValueOfErrors(unsigned int * syndrome, unsigned int * positions, int number_of_errors, unsigned int * values);
static Matrix calcSyndromeMatrix(unsigned int* syndrome, int len);
static Poly calcE(unsigned int * positions, unsigned int * values, unsigned int number_of_errors);
static Matrix calcErrorEquationsMatrix(unsigned int * syndrome, unsigned int * positions, int number_of_errors);

void RS_Init(unsigned int g, unsigned int N, unsigned int K){
	GF_Init(g);
	s_N = N;
	s_K = K;
	s_t = (N-K)/2;
}

int RS_Encode(unsigned int* plain, unsigned int len, unsigned int* encoded){
	Poly I=NULL,G=NULL,x_N_K=NULL,Ix_N_K=NULL, P=NULL, C=NULL;
	I = calcI(plain,len);
	if(!I)goto cleanup;
	G = calcG();
	if(!G)goto cleanup;
	x_N_K = Poly_GetMono(1,s_N-s_K);
	if(!x_N_K)goto cleanup;
	Ix_N_K = Poly_Mul(x_N_K,I);
	if(!Ix_N_K)goto cleanup;
	P = Poly_Mod(Ix_N_K,G);
	if(!P)goto cleanup;
	C = Poly_Add(Ix_N_K,P);
	if(!C)goto cleanup;
	
	for(unsigned int i=0;i<s_N;i++){
		encoded[i]=Poly_GetCoef(C,s_N-1-i);
	}
	return 0;
cleanup:
	Poly_Free(I);
	Poly_Free(G);
	Poly_Free(x_N_K);
	Poly_Free(Ix_N_K);
	Poly_Free(P);
	Poly_Free(C);
	return -1;
}

int RS_Decode(unsigned int* encoded, unsigned int len, unsigned int* plain){
	Poly Y=NULL, E=NULL, C=NULL;
	unsigned int * syndrome=NULL;
	unsigned int * positions = NULL;
	unsigned int * values=NULL;

	Y=calcI(encoded, len);
	if(!Y)goto cleanup;
	syndrome = (unsigned int *) malloc(sizeof(unsigned int )*2*s_t);
	if(!syndrome)goto cleanup;
	calcSyndrome(Y,syndrome);
	int k = calcNumberOfErrors(syndrome);
	if(k<0)goto cleanup;
	if(k==0){
		memcpy(plain, encoded, sizeof(unsigned int)*len);
		return 0;
	}
	positions = (unsigned int *) malloc(sizeof(unsigned int )*k);
	if(!positions)goto cleanup;
	if(calcPositionOfErrors(syndrome,k,positions))goto cleanup;
	values = (unsigned int *) malloc(sizeof(unsigned int )*k);
	if(!values)goto cleanup;
	if(calcValueOfErrors(syndrome, positions, k, values))goto cleanup;
	E = calcE(positions, values, k);
	if(!E)goto cleanup;
	C = Poly_Add(Y,E);
	if(!C)goto cleanup;
	for(unsigned int i=0;i<s_N;i++){
		plain[i]=Poly_GetCoef(C,s_N-i-1);
	}

	free(positions);
	free(syndrome);
	free(values);
	Poly_Free(Y);
	Poly_Free(E);
	Poly_Free(C);
	return 0 ;
cleanup:
	Poly_Free(C);
	Poly_Free(E);
	Poly_Free(Y);
	free(values);
	free(positions);
	free(syndrome);
	return -1;
}

Poly calcI(unsigned int* data, unsigned int len){
	Poly I=NULL, p=NULL, tmp=NULL;
	I = Poly_GetMono(0,0);
	if(!I)goto cleanup;
	for(unsigned int i = len;i-->0;){
		p = Poly_GetMono(data[i],len-1-i);
		if(!p)goto cleanup;
		tmp = Poly_Add(I,p);
		Poly_Free(I);
		Poly_Free(p);
		I = tmp;
	}
	return I;
cleanup:
	Poly_Free(I);
	Poly_Free(p);
	Poly_Free(tmp);
	return NULL;
}

Poly calcG(){
	Poly G=NULL, factor=NULL, tmp=NULL;
	G = Poly_GetMono(1,0);
	if(!G)goto cleanup;
	for(unsigned int i =0;i<2*s_t;i++){
		factor = calcGFactor(i);
		if(!factor)goto cleanup;
		tmp = Poly_Mul(G,factor);
		if(!tmp)goto cleanup;
		Poly_Free(G);
		Poly_Free(factor);
		G=tmp;
	}
	return G;
cleanup:
	Poly_Free(G);
	Poly_Free(factor);
	Poly_Free(tmp);
	return NULL;
}

Poly calcGFactor(unsigned int i){
	Poly x=NULL, c=NULL, result=NULL;
	x = Poly_GetMono(1,1);
	if(!x)goto cleanup;
	c = Poly_GetMono(GF_Pow(i),0);
	if(!c)goto cleanup;
	result = Poly_Sub(x,c);
	if(!result)goto cleanup;
	return result;
cleanup:
	Poly_Free(x);
	Poly_Free(c);
	Poly_Free(result);
	return NULL;
}


void calcSyndrome(Poly Y, unsigned int * syndrome){
	for(unsigned int i=0;i<2*s_t;i++){
		syndrome[i] = Poly_Substitute(Y,GF_Pow(i));
	}
}

static int calcNumberOfErrors(unsigned int* syndrome){
	Matrix m=NULL, tmp=NULL;
	for(unsigned int k=s_t;k>0;k--){
		tmp=calcSyndromeMatrix(syndrome,k);
		if(!tmp)return -1;
		m= Matrix_Inv(tmp);
		Matrix_Free(tmp);
		if (m!=NULL){
			Matrix_Free(m);
			return k;
		}
	}
	return 0;
}
static Matrix calcSyndromeMatrix(unsigned int* syndrome, int k){
	unsigned int ** elems =NULL;
	Matrix result = NULL;
	elems = (unsigned int **)malloc(sizeof(unsigned int *)*k);
	if(!elems)return NULL;
	for ( int i= 0; i<k;i++){
		elems[i]= (unsigned int *)malloc(sizeof(unsigned int)*k);
		if(!elems[i]){
			for(int j=0;j<i;j++){
				free(elems[j]);
			}
			free(elems);
			return NULL;
		}
		for(int j=0;j<k;j++){
			elems[i][j]=syndrome[k-1+i-j];
		}
	}
	result = Matrix_Init(elems,k,k);
	for(int i=0;i<k;i++){
		free(elems[i]);
	}
	free(elems);
	return result;
}


static int calcPositionOfErrors(unsigned int * syndrome, int number_of_errors, unsigned int* positions){
	Matrix S = NULL, S_inv=NULL;
	Poly p=NULL;
	unsigned int * sigmas=NULL;
	S = calcSyndromeMatrix(syndrome,number_of_errors);
	if(!S)goto cleanup;
	S_inv = Matrix_Inv(S);
	if(!S_inv)goto cleanup;
	sigmas = malloc(sizeof(unsigned int ) * (number_of_errors+1) );
	if(!sigmas)goto cleanup;
	sigmas[0]=1;
	int count = 0;
	Matrix_Mul(S_inv,syndrome+number_of_errors,number_of_errors,sigmas+1);
	p = Poly_Init(sigmas, number_of_errors +1);
	if(!p)goto cleanup;
	for(int i=0;i<GF_Order()-1;i++){
		if(Poly_Substitute(p,GF_Pow(i))==0){
			positions[count++]=GF_Order()-1-i;
		}
	}
	assert(number_of_errors== count);

  	Matrix_Free(S);
	Matrix_Free(S_inv);
	Poly_Free(p);
	free(sigmas);

	return 0;
cleanup:
	Matrix_Free(S);
	Matrix_Free(S_inv);
	Poly_Free(p);
	free(sigmas);
	return -1;
}



static int calcValueOfErrors(unsigned int * syndrome, unsigned int * positions, int number_of_errors, unsigned int * values){
	Matrix mat =calcErrorEquationsMatrix(syndrome,positions,number_of_errors);
	if(!mat)return -1;
	Matrix imat = Matrix_Inv(mat);
	Matrix_Free(mat);
	if(!imat)return -1;
	Matrix_Mul(imat,syndrome,number_of_errors,values);
	Matrix_Free(imat);
	return 0;
}

static Matrix calcErrorEquationsMatrix(unsigned int * syndrome, unsigned int * positions, int number_of_errors){
	unsigned int ** elems = malloc(sizeof(unsigned int *)*number_of_errors);
	if(!elems)return NULL;
	for (unsigned int i = 0; i< number_of_errors;i++){
		elems[i]= malloc(sizeof(unsigned int )* number_of_errors);
		if(!elems[i]){
			for (unsigned int j =0;j<i;j++){
				free(elems[j]);
			}
			free(elems);
			return NULL;
		}
		for (unsigned int j =0;j<number_of_errors;j++){
			elems[i][j]= GF_Pow(positions[j]*i);
		}
	}
	Matrix mat = Matrix_Init(elems,number_of_errors,number_of_errors);
	for (unsigned int i =0;i<number_of_errors;i++){
		free(elems[i]);
	}
	free(elems);
	return mat;
}


static Poly calcE(unsigned int * positions, unsigned int * values, unsigned int number_of_errors){
	Poly E=NULL, tmp=NULL, tmp2=NULL;
	E  = Poly_GetMono(0,0);
	if(!E)goto cleanup;
	for(unsigned int i=0;i<number_of_errors;i++){
		tmp  = Poly_GetMono(values[i],positions[i]);
		if(!tmp)goto cleanup;
		tmp2 = Poly_Add(E, tmp);
		if(!tmp2)goto cleanup;
		Poly_Free(E);
		Poly_Free(tmp);
		E=tmp2;
	}
	return E;
cleanup:
	Poly_Free(E);
	Poly_Free(tmp);
	Poly_Free(tmp2);
}
