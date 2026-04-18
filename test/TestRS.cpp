#include "CppUTest/TestHarness.h"
extern "C" {
 #include "RS.h"
 #include "Matrix.h"
 #include "GF.h"
}


TEST_GROUP(RS_GF2_4)
{
	void setup(){
		
	}
	void teardown(){

	}
};

TEST(RS_GF2_4, Encode)
{
	RS_Context ctx=RS_Init(0x13,6,2);
	unsigned int plain[2];
	plain[0]=GF_Pow(1);
	plain[1]=GF_Pow(2);
	RS_Message message = {plain,2};
	unsigned int encoded[6];
	RS_Message encoded_message = {encoded, 6};
	RS_Encode(ctx,&message,&encoded_message);
	

	LONGS_EQUAL(GF_Pow(1),encoded[0]);
	LONGS_EQUAL(GF_Pow(2),encoded[1]);
	LONGS_EQUAL(GF_Pow(3),encoded[2]);
	LONGS_EQUAL(GF_Pow(9),encoded[3]);
	LONGS_EQUAL(GF_Pow(1),encoded[4]);
	LONGS_EQUAL(GF_Pow(5),encoded[5]);

};

TEST(RS_GF2_4, Decode)
{
	RS_Context ctx=RS_Init(0x13,6,2);
	unsigned int plain[2];
	plain[0]=GF_Pow(1);
	plain[1]=GF_Pow(2);
	RS_Message message = {plain,2};
	unsigned int encoded[6];
	RS_Message encoded_message = {encoded, 6};
	RS_Encode(ctx,&message,&encoded_message);
	

	LONGS_EQUAL(GF_Pow(1),encoded[0]);
	LONGS_EQUAL(GF_Pow(2),encoded[1]);
	LONGS_EQUAL(GF_Pow(3),encoded[2]);
	LONGS_EQUAL(GF_Pow(9),encoded[3]);
	LONGS_EQUAL(GF_Pow(1),encoded[4]);
	LONGS_EQUAL(GF_Pow(5),encoded[5]);

	encoded[0]=GF_Pow(2);
	encoded[1]=GF_Pow(3);
	
	unsigned int decoded[6];
	RS_Message decoded_message = {decoded, 6};
	int result = RS_Decode(ctx, &encoded_message,&decoded_message);

	LONGS_EQUAL(GF_Pow(1),decoded[0]);
	LONGS_EQUAL(GF_Pow(2),decoded[1]);
	LONGS_EQUAL(GF_Pow(3),decoded[2]);
	LONGS_EQUAL(GF_Pow(9),decoded[3]);
	LONGS_EQUAL(GF_Pow(1),decoded[4]);
	LONGS_EQUAL(GF_Pow(5),decoded[5]);



};


TEST_GROUP(RS_GF2_8)
{
	void setup(){

	}
	void teardown(){

	}
};

TEST(RS_GF2_8, Decode)
{
	RS_Context ctx=RS_Init(0x11D,9,5);
	// Helloのメッセージを作成
	unsigned int plain[5];
	plain[0]=GF_Pow('H');
	plain[1]=GF_Pow('e');
	plain[2]=GF_Pow('l');
	plain[3]=GF_Pow('l');
	plain[4]=GF_Pow('o');
	RS_Message message={plain, 5};

	//RS符号化メッセージを作成
	unsigned int encoded[9];
	RS_Message encoded_message={encoded,9};
	RS_Encode(ctx, &message,&encoded_message);

	//通信路でメッセージが破損
	encoded_message.data[0]=1;
	encoded_message.data[2]=2;

	//RS復号を実施
	unsigned int decoded[9];
	RS_Message decoded_message={decoded,9};
	RS_Decode(ctx, &encoded_message,&decoded_message);

	// 復号されているか確認
	LONGS_EQUAL(GF_Pow('H'),decoded[0]);
	LONGS_EQUAL(GF_Pow('e'),decoded[1]);
	LONGS_EQUAL(GF_Pow('l'),decoded[2]);
	LONGS_EQUAL(GF_Pow('l'),decoded[3]);
	LONGS_EQUAL(GF_Pow('o'),decoded[4]);

};

