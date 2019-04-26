#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <vector> 
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <cmath>
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

//--------- Parse_Str_Str --------//
int Parse_Str_Str(string &in,vector <string> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		if(buf=="")continue;
		out.push_back(buf);
		count++;
	}
	return count;
}
//--------- Parse_Str_Str (automatic) --------//
int Parse_Str_Str(string &in,vector <string> &out)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(! (www>>buf) )break;
		if(buf=="")continue;
		out.push_back(buf);
		count++;
	}
	return count;
}


//-------- load FASTA file -------//
int Read_FASTA_SEQRES(string &seqfile,string &seqres,int skip=1) //->from .fasta file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(seqfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",seqfile.c_str());
		exit(-1);
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file bad! %s \n",seqfile.c_str());
			exit(-1);
		}
	}
	//process
	temp="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		temp+=buf;
	}
	seqres=temp;
	//return
	return (int)seqres.length();
}


//----------- load Phobius PLP file -------------//
//-> example
/*
# Q8TCT8|PSL2_HUMAN
#         i         o         O         M         n         h         c         C
M      0.00029   0.00004   0.00005   0.00000   0.99962   0.00000   0.00000   0.00000
G      0.00029   0.00004   0.00005   0.00000   0.99962   0.00000   0.00000   0.00000
P      0.00029   0.00004   0.00005   0.00000   0.99962   0.00000   0.00000   0.00000
Q      0.00029   0.00004   0.00005   0.00000   0.99962   0.00000   0.00000   0.00000
R      0.00029   0.00004   0.00005   0.00000   0.99962   0.00000   0.00000   0.00000
R      0.00029   0.00004   0.00005   0.00000   0.99925   0.00037   0.00000   0.00000
...
*/
int Load_Phobius_PLP(string &fn, vector <int> &sp_label)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(fn.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",fn.c_str());
		exit(-1);
	}
	//load
	int col_size=9;
	int count=0;
	int first=1;
	sp_label.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		vector <string> tmp_str;
		int retv=Parse_Str_Str(buf,tmp_str);
		if(retv!=col_size)
		{
			fprintf(stderr,"retv %d not equal to col_size %d \n",
				retv,col_size);
			exit(-1);
		}
		//assign SP label
		if(first==1)
		{
			//-> get SP probability
			double nt_prob=atof(tmp_str[5].c_str());
			double sp_prob=atof(tmp_str[6].c_str())+atof(tmp_str[7].c_str());
			double cut_prob=atof(tmp_str[8].c_str());
			double tot_prob=nt_prob+sp_prob+cut_prob+cut_prob;
			//-> check cut
			if(cut_prob>0.5 || tot_prob<0.5)
			{
				first=0;
				sp_label.push_back(0);
			}
			//-> check SP
			else
			{
				if(sp_prob>0.5)
					sp_label.push_back(1);
				else
					sp_label.push_back(0);
			}
			count++;
		}
		else
		{
			sp_label.push_back(0);
			count++;
		}
	}
	//return
	return count;
}


//=========== merge SP and TM together =============//
//[note]: by default, tm_skip_lines shall be 1 or 2
void Merge_SP_TM(string &sp_file, string &tm_file, int tm_skip_lines)
{
	//-> load SP_file
	vector <int> sp_label;
	int sp_len=Load_Phobius_PLP(sp_file, sp_label);
	//-> load TM_label
	string tm_string;
	int tm_len=Read_FASTA_SEQRES(tm_file,tm_string,tm_skip_lines);
	//-> check length
	if(sp_len!=tm_len)
	{
		fprintf(stderr,"sp_len %d not equal to tm_len %d\n",
			sp_len,tm_len);
		exit(-1);
	}
	//-> merge SP and TM
	vector <int> tm_label(tm_len,0);
	for(int i=0;i<tm_len;i++)tm_label[i]=tm_string[i]-'0';
	for(int i=0;i<sp_len;i++)if(sp_label[i]==1 && tm_label[i]==0)tm_label[i]=2;
	//-> output
	for(int i=0;i<tm_len;i++)printf("%d",tm_label[i]);
	printf("\n");
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//---- Merge_SP_TM ----//
	{
		if(argc<4)
		{
			fprintf(stderr,"Merge_SP_TM <sp_file> <tm_label> <tm_skip_lines> \n");
			fprintf(stderr,"[note]: tm_skip_lines shall be set to 1 or 2 \n");
			exit(-1);
		}
		string sp_file=argv[1];
		string tm_label=argv[2];
		int tm_skip_lines=atoi(argv[3]);
		//process
		Merge_SP_TM(sp_file,tm_label,tm_skip_lines);
		//exit
		exit(0);
	}
}

