#include "RS.h"
#include "GF.h"
#include "Poly.h"
#include "Matrix.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

struct RS_ContextStruct{
	unsigned int N;
	unsigned int K;
	unsigned int t;
};

static Poly correctErrors(RS_Context ctx, Poly y, const unsigned int* syndrome);
static Poly buildMessagePolynomial(const RS_Message* message);
static Poly buildGeneratorPolynomial(RS_Context ctx );
static Poly buildGeneratorFactor(unsigned int i);
static unsigned int* computeSyndromes(RS_Context ctx, Poly Y);
static int estimateErrorCount(RS_Context ctx, const unsigned int* syndrome);
static unsigned int* findErrorLocations(const unsigned int * syndrome, int number_of_errors);
static unsigned int* computeErrorValues(const unsigned int * syndrome, const unsigned int * positions, int number_of_errors);
static Matrix buildSyndromeMatrix(const unsigned int* syndrome, int len);
static Poly buildErrorPolynomial(const unsigned int * positions, const unsigned int * values, unsigned int number_of_errors);
static Matrix buildErrorEquationsMatrix(const unsigned int * syndrome, const unsigned int * positions, int number_of_errors);
static Poly buildMonomial(unsigned int coef, unsigned int degree);
static void writeMessage(RS_Context ctx, Poly p, RS_Message* out);

RS_Context RS_Init(unsigned int g, unsigned int N, unsigned int K){
	RS_Context ctx = (RS_Context)malloc(sizeof(struct RS_ContextStruct));
	if(!ctx)return NULL;
	GF_Init(g);
	ctx->N= N;
	ctx->K = K;
	ctx->t = (N-K)/2;
	return ctx;
}

int RS_Encode(RS_Context ctx, const RS_Message* message, RS_Message* encoded_message){
	Poly I=NULL,G=NULL,x_N_K=NULL, P=NULL, C=NULL, tmp=NULL;
	//情報多項式Iの生成
	I = buildMessagePolynomial(message);
	if(!I)goto cleanup;

	// 生成多項式Gの生成
	G = buildGeneratorPolynomial(ctx);
	if(!G)goto cleanup;

	// 符号語に対応する多項式Cの生成
	x_N_K = buildMonomial(1,ctx->N-ctx->K);
	if(!x_N_K)goto cleanup;
	tmp = Poly_Mul(x_N_K,I);
	if(!tmp)goto cleanup;
	P = Poly_Mod(tmp,G);
	if(!P)goto cleanup;
	C = Poly_Add(tmp,P);
	if(!C)goto cleanup;
	
	// 多項式Cから符号語を書きだし
	writeMessage(ctx, C,encoded_message);
	return 0;
cleanup:
	Poly_Free(I);
	Poly_Free(G);
	Poly_Free(x_N_K);
	Poly_Free(tmp);
	Poly_Free(P);
	Poly_Free(C);
	return -1;
}

int RS_Decode(RS_Context ctx, const RS_Message* message, RS_Message* decoded_message){
	Poly Y=NULL,  C=NULL;
	unsigned int * syndrome=NULL;

	// 受信したデータから多項式Yを生成
	Y=buildMessagePolynomial(message);
	if(!Y)goto cleanup;

	// シンドロームの算出
	syndrome=computeSyndromes(ctx, Y);
	if(!syndrome)goto cleanup;

	// 誤り訂正を行い、復元した符号語に対応する多項式Cを生成
	C = correctErrors(ctx, Y,syndrome);
	if(!C)goto cleanup;

	// 多項式Cから復元した符号語を書きだし
	writeMessage(ctx,C,decoded_message);

	free(syndrome);
	Poly_Free(Y);
	Poly_Free(C);
	return 0 ;
cleanup:
	Poly_Free(C);
	Poly_Free(Y);
	free(syndrome);
	return -1;
}

void RS_Free(RS_Context ctx){
		free(ctx);
} 

static Poly correctErrors(RS_Context ctx, Poly y, const unsigned int* syndrome){
	Poly E=NULL, C=NULL;
	unsigned int * positions = NULL;
	unsigned int * values=NULL;

	// 誤りを含むシンボルの数を検出
	int k = estimateErrorCount(ctx, syndrome);
	if(k<0)goto cleanup;
	if(k==0){
		return Poly_Copy(y);
	}
	
	// 誤りを含むシンボルの位置を検出
	positions=findErrorLocations(syndrome,k);
	if(!positions)goto cleanup;

	// 誤りの値を検出
	values=computeErrorValues(syndrome, positions, k);
	if(!values)goto cleanup;

	// 誤りの訂正
	E = buildErrorPolynomial(positions, values, k);
	if(!E)goto cleanup;
	C = Poly_Add(y,E);
	if(!C)goto cleanup;

	Poly_Free(E);
	free(positions);
	free(values);
	return C;

cleanup:
	Poly_Free(E);
	free(positions);
	free(values);
	return NULL;
}

static void writeMessage(RS_Context ctx, Poly p, RS_Message* out){
	for(unsigned int i=0;i<ctx->N;i++){
		out->data[i]=Poly_Data(p)[ctx->N-1-i];
	}
}


static Poly buildMessagePolynomial(const RS_Message* message){
	Poly I=NULL, p=NULL, tmp=NULL;
	I = buildMonomial(0,0);
	if(!I)goto cleanup;
	for(size_t i = message->length;i-->0;){
		p = buildMonomial(message->data[i],message->length-1-i);
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

static Poly buildGeneratorPolynomial(RS_Context ctx){
	Poly G=NULL, factor=NULL, tmp=NULL;
	G = buildMonomial(1,0);
	if(!G)goto cleanup;
	for(unsigned int i =0;i<2*ctx->t;i++){
		factor = buildGeneratorFactor(i);
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

static Poly buildGeneratorFactor(unsigned int i){
	Poly x=NULL, c=NULL, result=NULL;
	x = buildMonomial(1,1);
	if(!x)goto cleanup;
	c = buildMonomial(GF_Pow(i),0);
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

static unsigned int * computeSyndromes(RS_Context ctx, Poly Y){
	unsigned int * syndrome = (unsigned int *) malloc(sizeof(unsigned int )*2*ctx->t);
	if(!syndrome)return NULL;
	for(unsigned int i=0;i<2*ctx->t;i++){
		syndrome[i] = Poly_Eval(Y,GF_Pow(i));
	}
	return syndrome;
}

static int estimateErrorCount(RS_Context ctx, const unsigned int* syndrome){
	Matrix m=NULL, tmp=NULL;
	for(unsigned int k=ctx->t;k>0;k--){
		tmp=buildSyndromeMatrix(syndrome,k);
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

static Matrix buildSyndromeMatrix(const unsigned int* syndrome, int k){
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

static unsigned int* findErrorLocations(const unsigned int * syndrome, int number_of_errors){
	Matrix S = NULL, S_inv=NULL;
	Poly p=NULL;
	unsigned int * sigmas=NULL;
	unsigned int * positions = (unsigned int *) malloc(sizeof(unsigned int )*number_of_errors);
	if(!positions)return NULL;
	S = buildSyndromeMatrix(syndrome,number_of_errors);
	if(!S)goto cleanup;
	S_inv = Matrix_Inv(S);
	if(!S_inv)goto cleanup;
	sigmas = malloc(sizeof(unsigned int ) * (number_of_errors+1) );
	if(!sigmas)goto cleanup;
	sigmas[0]=1;
	int count = 0;
	unsigned int * v=Matrix_Mul(S_inv,syndrome+number_of_errors,number_of_errors);
	if(!v)goto cleanup;
	memcpy(sigmas+1,v,sizeof(unsigned int)*number_of_errors);
	free(v);
	p = Poly_Init(sigmas, number_of_errors +1);
	if(!p)goto cleanup;
	for(int i=0;i<GF_Order()-1;i++){
		if(Poly_Eval(p,GF_Pow(i))==0){
			positions[count++]=GF_Order()-1-i;
		}
	}
	assert(number_of_errors== count);

  	Matrix_Free(S);
	Matrix_Free(S_inv);
	Poly_Free(p);
	free(sigmas);

	return positions;
cleanup:
	Matrix_Free(S);
	Matrix_Free(S_inv);
	Poly_Free(p);
	free(sigmas);
	free(positions);
	return NULL;
}

static unsigned int* computeErrorValues(const unsigned int * syndrome, const unsigned int * positions, int number_of_errors){
	Matrix mat =buildErrorEquationsMatrix(syndrome,positions,number_of_errors);
	if(!mat)return NULL;
	Matrix imat = Matrix_Inv(mat);
	Matrix_Free(mat);
	if(!imat)return NULL;
	unsigned int * result=Matrix_Mul(imat,syndrome,number_of_errors);
	Matrix_Free(imat);
	return result;
}

static Matrix buildErrorEquationsMatrix(const unsigned int * syndrome, const unsigned int * positions, int number_of_errors){
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

static Poly buildErrorPolynomial(const unsigned int * positions, const unsigned int * values, unsigned int number_of_errors){
	Poly E=NULL, tmp=NULL, tmp2=NULL;
	E  = buildMonomial(0,0);
	if(!E)goto cleanup;
	for(unsigned int i=0;i<number_of_errors;i++){
		tmp  = buildMonomial(values[i],positions[i]);
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
	return NULL;
}

static Poly buildMonomial(unsigned int coef, unsigned int degree){
        Poly p = NULL;
        unsigned int* coefs = (unsigned int *)malloc(sizeof(unsigned int)*(degree+1));
        if(!coefs)return NULL;
        memset(coefs,0,sizeof(unsigned int)*(degree+1));
        coefs[degree]=coef;
        p= Poly_Init(coefs, degree+1);
        free(coefs);
        if(!p)return NULL;
        return p;
}
