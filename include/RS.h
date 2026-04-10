#ifndef RS_H
#define RS_H
void RS_Init(unsigned int g ,unsigned int N, unsigned int K);
int RS_Encode(unsigned int* plain, unsigned int len, unsigned int* encrypted);
int RS_Decode(unsigned int* encrypted, unsigned int len, unsigned int* decoded);
#endif
