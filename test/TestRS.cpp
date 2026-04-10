#include "CppUTest/TestHarness.h"
extern "C" {
 #include "RS.h"
 #include "Matrix.h"
 #include "GF.h"
}


TEST_GROUP(RS)
{
	void setup(){
		RS_Init(0x13,6,2);
	}
	void teardown(){

	}
};

TEST(RS, Encode)
{
	unsigned int plain[2];
	plain[0]=GF_Pow(1);
	plain[1]=GF_Pow(2);
	unsigned int encoded[6];
	RS_Encode(plain,2,encoded);
	

	LONGS_EQUAL(GF_Pow(1),encoded[0]);
	LONGS_EQUAL(GF_Pow(2),encoded[1]);
	LONGS_EQUAL(GF_Pow(3),encoded[2]);
	LONGS_EQUAL(GF_Pow(9),encoded[3]);
	LONGS_EQUAL(GF_Pow(1),encoded[4]);
	LONGS_EQUAL(GF_Pow(5),encoded[5]);

};

TEST(RS, Decode)
{
	unsigned int plain[2];
	plain[0]=GF_Pow(1);
	plain[1]=GF_Pow(2);
	unsigned int encoded[6];
	RS_Encode(plain,2,encoded);
	

	LONGS_EQUAL(GF_Pow(1),encoded[0]);
	LONGS_EQUAL(GF_Pow(2),encoded[1]);
	LONGS_EQUAL(GF_Pow(3),encoded[2]);
	LONGS_EQUAL(GF_Pow(9),encoded[3]);
	LONGS_EQUAL(GF_Pow(1),encoded[4]);
	LONGS_EQUAL(GF_Pow(5),encoded[5]);

	encoded[0]=GF_Pow(2);
	encoded[1]=GF_Pow(3);
	
	unsigned int decoded[6];
	int result = RS_Decode(encoded, 6,decoded);

	LONGS_EQUAL(GF_Pow(1),decoded[0]);
	LONGS_EQUAL(GF_Pow(2),decoded[1]);
	LONGS_EQUAL(GF_Pow(3),decoded[2]);
	LONGS_EQUAL(GF_Pow(9),decoded[3]);
	LONGS_EQUAL(GF_Pow(1),decoded[4]);
	LONGS_EQUAL(GF_Pow(5),decoded[5]);



};

