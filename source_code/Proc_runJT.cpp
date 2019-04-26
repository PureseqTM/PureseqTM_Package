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


//------------- load runJT_test.out --------------//
//-> file format
/*
Segment 0, after CE, log(prob(evidence)) = -1478.685031, per frame =-2.843625, per numUFrams = -2.843625
--------
Partition 0 (P), Clique 1: Printing Clique with 1 variables, 3 entries, H=3.111350e-02
0: 9.96987843e-01 topo10(0)=9
1: 3.45358649e-04 topo10(0)=8
2: 2.66679872e-03 topo10(0)=6
--------
Partition 0 (P), Clique 2: Printing Clique with 1 variables, 4 entries, H=3.111350e-02
0: 7.15843773e-11 topo10(1)=7
1: 2.66679864e-03 topo10(1)=1
2: 9.96987843e-01 topo10(1)=4
3: 3.45358649e-04 topo10(1)=3
--------
Partition 1 (C), Clique 1: Printing Clique with 1 variables, 5 entries, H=3.113115e-02
0: 1.43698618e-06 topo10(2)=7
1: 2.66536842e-03 topo10(2)=1
2: 7.15843773e-11 topo10(2)=2
3: 3.45351891e-04 topo10(2)=3
4: 9.96987843e-01 topo10(2)=4
--------
Partition 2 (C), Clique 1: Printing Clique with 1 variables, 5 entries, H=3.116164e-02
0: 2.73295395e-06 topo10(3)=7
1: 2.66283599e-03 topo10(3)=1
2: 1.43705776e-06 topo10(3)=2
3: 3.45151360e-04 topo10(3)=3
4: 9.96987843e-01 topo10(3)=4
.....

Partition 11 (C), Clique 1: Printing Clique with 1 variables, 5 entries, H=3.300772e-02
0: 1.58953501e-04 topo10(12)=7
1: 2.34795890e-04 topo10(12)=1
2: 2.30849280e-03 topo10(12)=2
3: 3.09915177e-04 topo10(12)=3
4: 9.96987843e-01 topo10(12)=4
--------
Partition 12 (C), Clique 1: Printing Clique with 1 variables, 6 entries, H=3.239088e-02
0: 3.19497854e-11 topo10(13)=8
1: 1.33237712e-04 topo10(13)=1
2: 2.46744630e-03 topo10(13)=2
3: 3.08601338e-04 topo10(13)=3
4: 9.96987843e-01 topo10(13)=4
5: 1.02872017e-04 topo10(13)=7
--------
Partition 13 (C), Clique 1: Printing Clique with 1 variables, 6 entries, H=3.187784e-02
0: 7.92366127e-11 topo10(14)=8
1: 1.00013883e-04 topo10(14)=1
2: 2.57031832e-03 topo10(14)=2
3: 3.07963549e-04 topo10(14)=3
4: 9.96987843e-01 topo10(14)=4
5: 3.38616491e-05 topo10(14)=7
--------
Partition 14 (C), Clique 1: Printing Clique with 1 variables, 6 entries, H=3.167804e-02
0: 3.65897066e-10 topo10(15)=8
1: 8.56812642e-05 topo10(15)=1
2: 2.60417996e-03 topo10(15)=2
3: 3.07605413e-04 topo10(15)=3
4: 9.96987842e-01 topo10(15)=4
5: 1.46908345e-05 topo10(15)=7

.....

Partition 518 (E), Clique 1: Printing Clique with 1 variables, 4 entries, H=2.776011e-02
0: 3.72785395e-06 topo10(519)=8
1: 9.97207069e-01 topo10(519)=1
2: 8.90806984e-09 topo10(519)=6
3: 2.78919459e-03 topo10(519)=3
--------
Partition 518 (E), Clique 2: Printing Clique with 1 variables, 4 entries, H=2.956585e-02
0: 3.01215736e-03 pType(519)=2
1: 9.96987843e-01 pType(519)=3
2: 1.28018460e-22 pType(519)=0
3: 5.59126064e-19 pType(519)=1
____ PROGRAM ENDED SUCCESSFULLY WITH STATUS 0 AT Wednesday April 03 2019, 18:33:10 +03 ____
*/


//--- index convert ----//
int index_convert(int &in)
{
	switch(in)
	{
		case 1: return 0;
		case 2: return 1;
		case 3: return 2;
		case 4: return 3;
		case 6: return 4;
		case 7: return 5;
		case 8: return 6;
		case 9: return 7;    //-> this is for START
		case 0: return 8;    //-> this is for TERMI
		default: return -1;
	}
}

//--- process runJT_test each part ---//
int Proc_runJT_single(vector <vector <string> > &in, vector <float> &out)
{
	int i;
	int size=(int)in.size();
	int colsize=7;
	int normal=1;
	out.resize(colsize,0);
	for(i=0;i<size;i++)
	{
		//position
		int pos=atoi(in[i][2].substr(in[i][2].length()-1,1).c_str());
		int index=index_convert(pos);
		if(index>=colsize)
		{
			if(index==8)normal=0;
			continue;
		}
		if(pos<0)
		{
			fprintf(stderr,"pos %d over-range \n",pos);
			exit(-1);
		}
		//value
		float val=atof(in[i][1].c_str());
		//assign
		out[index]=val;
	}
	return normal;
}

//--- the output shall be L*7 ----//
int Proc_runJT(string &fn, vector <vector <float> > &out)
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
	int col_size=3;
	int count=0;
	int first=1;
	out.clear();
	vector <vector <string> > tmp_rec_all;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		vector <string> tmp_str;
		int retv=Parse_Str_Str(buf,tmp_str);
		if(retv!=col_size)
		{
			if(first==0)
			{
				//-> process
				vector <float> val;
				int retv=Proc_runJT_single(tmp_rec_all,val);
				if(retv==0)break;
				//-> record
				out.push_back(val);
				count++;
				//-> clear
				tmp_rec_all.clear();
			}
			first=1;
		}
		else
		{
			tmp_rec_all.push_back(tmp_str);
			first=0;
		}
	}
	//return
	return count;
}

//------------ main -------------//
int main(int argc, char** argv)
{
	//---- Proc_runJT ----//
	{
		if(argc<3)
		{
			fprintf(stderr,"Proc_runJT <input_runJT> <output_file> \n");
			exit(-1);
		}
		string input_runJT=argv[1];
		string output_file=argv[2];
		//proc
		vector <vector <float> > out;
		int num=Proc_runJT(input_runJT, out);
		if(num<=0)exit(-1);
		//output
		int colsize=7;
		FILE *fp=fopen(output_file.c_str(),"wb");
		for(int i=0;i<num;i++)
		{
			for(int j=0;j<colsize;j++)
			{
				fprintf(fp,"%e ",out[i][j]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		//exit
		exit(0);
	}
}

