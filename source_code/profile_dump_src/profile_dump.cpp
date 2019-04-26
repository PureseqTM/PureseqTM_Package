#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <string>
#include "seq.h"
using namespace std;


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

//------------ Load TGT, output profile -------//
void profile_dump(string &tgt_file)
{
	//-- read tgt ---//
	string tgt_name,tgt_root;
	getBaseName(tgt_file,tgt_name,'/','.');
	getRootName(tgt_file,tgt_root,'/');
	SEQUENCE *s=new SEQUENCE(tgt_name,tgt_root,1,1);

	//======== process feature =======//
	vector <string> output;
	output.clear();
	for(int k=0;k<s->length;k++)
	{
		//----- feature -----//
		vector <double> features;
		features.clear();
		//-- profile realted --//
		//emission score
		for(int i=0;i<20;i++)
		{
			float template_amino_Score = 1.0*s->EmissionProb[k][i];
			features.push_back(template_amino_Score);
		}
		//PSM score
		for(int i=0;i<20;i++)
		{
			double template_amino_Prob = 1.0*s->PSM[k][i];
			template_amino_Prob=1.0/(1.0+exp(-1.0*template_amino_Prob));
			features.push_back(template_amino_Prob);
		}
		//------ output ------//
		int featdim=(int)features.size();
		stringstream oss;
		for(int i=0;i<featdim;i++)
		{
			int wsiii=(int)features[i];
			if(wsiii!=features[i])oss << features[i] << " ";
			else oss << wsiii << " ";
		}
		string wsbuf=oss.str();
		output.push_back(wsbuf);
	}

	//====== output =======//
	for(int k=0;k<s->length;k++)printf("%s\n",output[k].c_str());

	//---- delete ----//
	delete s;
}

//----------- main -------------//
int main(int argc,char **argv)
{
	//------- Dump 20+20 profile features from TGT file --------// 
	{
		if(argc<2)
		{
			fprintf(stderr,"profile_dump <tgt_file> \n");
			exit(-1);
		}
		string tgt_file=argv[1];
		//process
		profile_dump(tgt_file);
		exit(0);
	}
}
