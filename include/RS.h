#ifndef RS_H
#define RS_H
#include <string.h>
typedef struct{
    unsigned int* data;
    size_t length;
} RS_Message;
typedef struct RS_ContextStruct* RS_Context;
RS_Context RS_Init(unsigned int g ,unsigned int N, unsigned int K);
int RS_Encode(RS_Context ctx, const RS_Message* message, RS_Message* encoded_message);
int RS_Decode(RS_Context ctx, const RS_Message* message, RS_Message* decoded_message);
void RS_Free(RS_Context ctx);
#endif
