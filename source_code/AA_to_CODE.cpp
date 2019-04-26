#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


//======================= I/O related ==========================//
//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}

//----------------------//
//-- Getline_Ending ---//
//--------------------//
void getline_end(string &input,char kill)
{
	int len=(int)input.length();
	if(input[len-1]==kill)input=input.substr(0,len-1);
}

//---- kill space ----//
void Kill_Invalid(string &in,string &out)
{
	int i;
	int len=(int)in.length();
	out="";
	for(i=0;i<len;i++)
	{
		if(in[i]!=' ' && in[i]!='-' )out=out+in[i];
	}
}

//---------- from amino acid to 0-20 code -----------//
void AA_to_CODE(string &input_fasta, string &output_code)
{
// maxlength
	const int BufLen = 100000;
	string sequence_ = "";
	string sequence = "";
	char buf[BufLen];
	const int AA2SUB[26]={0,20,1,2,3,4,5,6,7,20,8,9,10,11,20,12,13,14,15,16,20,17,18,20,19,20};

// generate SEQ file
	string seqfile0=input_fasta;
	ifstream seq0_in(seqfile0.c_str());
	if(!seq0_in.is_open())
	{
		cerr << "Cannot open input_fasta_file " << seqfile0 << endl;
		exit(1);
	}
	int first=1;
	while(seq0_in.getline(buf,BufLen))
	{
		string wstmp=buf;
		if(wstmp=="")continue;
		getline_end(wstmp,0x0D);
		if(buf[0]!='>') sequence_ = sequence_ + wstmp;
		else
		{
			if(first==1)first=0;
			else break;
		}
	}
	Kill_Invalid(sequence_,sequence);
	int fullLen = sequence.length();
	seq0_in.close();
	
// fix unknown amino acids other than ARNDCQEGHILKMFPSTWYVX
	for(int i=0;i<fullLen;i++)
	{
		char c=sequence[i];
		if(c<'A' || c>'Z')
		{
			sequence[i]='X';
			continue;
		}
		if(AA2SUB[c-'A']==20)
		{
			sequence[i]='X';
			continue;
		}
	}
	
	//===== output 0-20 code ======//
	FILE *fpseq=fopen(output_code.c_str(),"wb");
	for(int i=0;i<fullLen;i++)
	{
		int code=AA2SUB[sequence[i]-'A'];
		if(code==20)code=0;
		else code++;
		fprintf(fpseq,"%2d  0 \n",code);
	}
	fprintf(fpseq,"\n");
	fclose(fpseq);
}


//----------- main ------------//
int main(int argc, char** argv)
{
	//------- AA_to_CODE -----//
	{
		if(argc<3)
		{
			fprintf(stderr,"AA_to_CODE <input_fasta> <output_code> \n");
			exit(-1);
		}
		string input_fasta=argv[1];
		string output_code=argv[2];
		AA_to_CODE(input_fasta,output_code);
		//exit
		exit(0);
	}
}
