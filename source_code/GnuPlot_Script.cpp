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


//-------- select segment --------//
void Select_Segment(vector <int> &input, int label,vector < pair<int,int> > & output, int thres=0)
{
	output.clear();
	int first=1;
	int start,end;
	for(int i=0;i<(int)input.size();i++)
	{
		if(input[i]==label)
		{
			if(first==1)
			{
				first=0;
				start=i;
			}
		}
		else
		{
			if(first==0)
			{
				first=1;
				end=i-1;
				if(end-start+1>=thres)
				{
					output.push_back(pair<int,int>(start,end));
				}
			}
		}
	}
	if(first==0)
	{
		first=1;
		end=(int)input.size()-1;
		if(end-start+1>=thres)
		{
			output.push_back(pair<int,int>(start,end));
		}
	}
}


//-------- load predicted results -------------//
/*
# PureseqTM: 2-state transmembrane topology prediction results
# Transmembrane (TM) residues are marked with '1' above threshold 0.500
# Non-transmembrane (non-TM) residues are marked with '0' below threshold 0.500
# If ground-truth is given, the TM (non-TM) residues are marked with '1' ('0') at the right-most column
   1 M 0 0.003 0
   2 G 0 0.000 0
   3 S 0 0.000 0
   4 S 0 0.000 0
   5 H 0 0.000 0
   6 H 0 0.000 0
...
*/
int Load_Predicted_Results(string &input_file, 
	vector <double> &out_reso,
	vector <int> &pred_label,
	vector <int> &truth_label)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",input_file.c_str());
		exit(-1);
	}
	//process
	out_reso.clear();
	pred_label.clear();
	truth_label.clear();
	int count=0;
	int first=1;
	int first_num;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		//load
		vector <string> tmp_str;
		int retv=Parse_Str_Str(buf,tmp_str);
		if(first==1)
		{
			first_num=retv;
			first=0;
		}
		else
		{
			if(retv!=first_num)
			{
				fprintf(stderr,"retv %d not equal to first_num %d \n",
					retv,first_num);
				exit(-1);
			}
		}
		//parse
		double score=atof(tmp_str[3].c_str());
		out_reso.push_back(score);
		int pred_lab=atoi(tmp_str[2].c_str());
		pred_label.push_back(pred_lab);
		if(first_num==5)
		{
			int truth_lab=atoi(tmp_str[4].c_str());
			truth_label.push_back(truth_lab);
		}
		count++;
	}
	//return
	return count;
}

//------------ output to GNUPLOT script --------------//
//-> example
/*
set arrow from  1,-0.05 to 11,-0.05 nohead lt 3 lw 10
set arrow from 12,-0.05 to 29,-0.05 nohead lt 1 lw 10
set arrow from 30,-0.05 to 40,-0.05 nohead lt 3 lw 10
set arrow from 41,-0.05 to 66,-0.05 nohead lt 1 lw 10
set arrow from 67,-0.05 to 71,-0.05 nohead lt 3 lw 10
set title "Phobius posterior probabilities for UNNAMED "
set key below
set yrange [-0.1:1.0]
set ylabel "Posterior label probability"
set xrange [1:71]
set terminal pngcairo
set output 'aaaa.png'
plot "plp" using 1:5 title "transmembrane" with line lt 1 lw 2, \
"" using 1:3 title "non-transmembrane" with line lt 3 lw 2
exit
*/

void GnuPlot_Script(string &infile, double thres, int truth, string &pngfile, string &gff_file)
{
	//-> get name
	string name;
	getBaseName(infile,name,'/','.');
	//-> load file
	vector <double> out_reso;
	vector <int> pred_label;
	vector <int> truth_label;
	int length=Load_Predicted_Results(infile,out_reso,pred_label,truth_label);
	//-> select segment
	vector < pair<int,int> > tm_segment;
	vector < pair<int,int> > nontm_segment;
	if(truth==1 && truth_label.size()==length)
	{
		Select_Segment(truth_label,1,tm_segment);
		Select_Segment(truth_label,0,nontm_segment);
	}
	else
	{
		Select_Segment(pred_label,1,tm_segment);
		Select_Segment(pred_label,0,nontm_segment);
	}
	//--------- write GNUPLOT script ---------//
	//-> set arrow
	//--| nontm_segment
	for(int k=0;k<(int)nontm_segment.size();k++)
	{
		int begin=nontm_segment[k].first;
		int termi=nontm_segment[k].second;
		for(int i=-10;i<=10;i++)
		{
			float posi=-0.05+i*0.001;
			printf("set arrow from  %d,%5.3f to %d,%5.3f nohead lt rgb \"grey\" lw 1 \n",
				begin+1,posi,termi+1,posi);
		}
		
	}
	//--| tm_segment
	for(int k=0;k<(int)tm_segment.size();k++)
	{
		int begin=tm_segment[k].first;
		int termi=tm_segment[k].second;
		for(int i=-10;i<=10;i++)
		{
			float posi=-0.05+i*0.001;
			printf("set arrow from  %d,%5.3f to %d,%5.3f nohead lt rgb \"red\" lw 1 \n",
				begin+1,posi,termi+1,posi);
		}
	}
	//-> draw threshold
	printf("set arrow from 1,%5.3f to %d,%5.3f nohead lt rgb \"green\" lw 1 \n",
		thres,length,thres);
	//-> set title
	printf("set title \"PureseqTM posterior probabilities for %s \" \n",
		name.c_str());
	//-> others
	printf("set key below \n");
	printf("set yrange [-0.1:1.0] \n");
	printf("set ylabel \"Posterior label probability\" \n");
	printf("set xrange [1:%d] \n",length);
	printf("set terminal pngcairo \n");
	printf("set output \"%s\" \n",pngfile.c_str());
	//-> plot command
	printf("plot \"%s\" using 1:4 title \"transmembrane\" with line lt rgb \"red\" lw 1, \"\" title \"non-transmembrane\" with line lt rgb \"grey\" lw 1 \n",
		infile.c_str());
	//-> end
	printf("exit \n");

	//============ write to GFF file ===========//
	FILE *fp=fopen(gff_file.c_str(),"wb");
	fprintf(fp,"ID\t%s\n",name.c_str());
	//--| nontm_segment
	for(int k=0;k<(int)nontm_segment.size();k++)
	{
		int begin=nontm_segment[k].first;
		int termi=nontm_segment[k].second;
		fprintf(fp,"FT\tTOPO_DOM\t%d\t%d\n",begin+1,termi+1);
	}
	//--| tm_segment
	for(int k=0;k<(int)tm_segment.size();k++)
	{
		int begin=tm_segment[k].first;
		int termi=tm_segment[k].second;
		fprintf(fp,"FT\tTRANSMEM\t%d\t%d\n",begin+1,termi+1);
	}
	fclose(fp);
}



//------------ main -------------//
int main(int argc, char** argv)
{
	//---- GnuPlot_Script ----//
	{
		if(argc<6)
		{
			fprintf(stderr,"GnuPlot_Script <pred_file> <thres> <truth> <png_file> <gff_file> \n");
			fprintf(stderr,"[note]: if truth is 1, then output the 'truth_label' \n");
			exit(-1);
		}
		string pred_file=argv[1];
		double thres=atof(argv[2]);
		int truth=atoi(argv[3]);
		string png_file=argv[4];
		string gff_file=argv[5];
		//process
		GnuPlot_Script(pred_file,thres,truth,png_file,gff_file);
		//exit
		exit(0);
	}
}

