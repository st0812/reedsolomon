#ifndef RS_H
#define RS_H
typedef struct RS_ContextStruct* RS_Context;
RS_Context RS_Init(unsigned int g ,unsigned int N, unsigned int K);
int RS_Encode(RS_Context ctx, const unsigned int* plain, unsigned int len, unsigned int* encrypted);
int RS_Decode(RS_Context ctx, const unsigned int* encrypted, unsigned int len, unsigned int* decoded);
void RS_Free(RS_Context ctx);
#endif
