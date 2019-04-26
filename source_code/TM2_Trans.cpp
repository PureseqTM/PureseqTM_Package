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

//-------- merge segment --------//
void Merge_Segment(vector < pair<int,int> > &input,vector < pair<int,int> > &output, int thres=0)
{
	//init judge
	output.clear();
	if(input.size()==0)return;
	//prev assign
	int start_prev=input[0].first;
	int end_prev=input[0].second;
	output.push_back(input[0]);
	int count=0;
	//proc
	for(int i=1;i<(int)input.size();i++)
	{
		//curr assign
		int start_curr=input[i].first;
		int end_curr=input[i].second;
		//judge
		if(start_curr-end_prev<thres) //-> merge
		{
			output[count].second=end_curr;
		}
		else                            //-> new assign
		{
			start_prev=input[i].first;
			end_prev=input[i].second;
			output.push_back(input[i]);
			count++;
		}
	}
}

//-------- find min and max ----//
//-> merge_thres could be set to 0.001 to 0.02
int Find_MinMax(vector <double> &score, vector <int> &min_pos, double merge_thres)
{
	int i;
	int size=(int)score.size();
	//-> merge
	vector <double> score_merge;
	vector <int> posi_merge;
	int count=0;
	score_merge.push_back(score[0]);
	posi_merge.push_back(0);
	count++;
	for(i=1;i<size;i++)
	{
		if(fabs(score[i]-score_merge[count-1])>merge_thres)
		{
			score_merge.push_back(score[i]);
			posi_merge.push_back(i);
			count++;
		}
	}
	//-> diff
	vector <double> diff(count-1,0);
	for(i=1;i<count;i++)diff[i-1]=score_merge[i]-score_merge[i-1];
	//-> minmax
	min_pos.clear();
	int min_count=0;
	for(i=5;i<(int)diff.size()-5;i++)
	{
		if(diff[i]*diff[i+1]<0)
		{
			if(diff[i]<0 && diff[i+1]>0) //min found
			{
				min_pos.push_back(posi_merge[i+1]);
				min_count++;
			}
		}
	}
	//return
	return min_count;
}


//-------- cut segment --------//
void Cut_Segment(vector < pair<int,int> > &input, vector <double> &score, int cutoff)
{
	//proc
	for(int i=0;i<(int)input.size();i++)
	{
		//curr assign
		int start=input[i].first;
		int end=input[i].second;
		int len=end-start+1;
		//judge
		if(len>=cutoff)
		{
			//-> get score
			vector <double> score_rec;
			for(int k=start;k<=end;k++)score_rec.push_back(score[k]);
			//-> find min-point
			vector <int> min_pos;
			double merge_thres=0.005;
			int min_count=Find_MinMax(score_rec,min_pos,merge_thres);
			//-> reduce the probablity
			if(min_count==1)score[start+min_pos[0]]=0;
			else if(min_count==2)
			{
				int cur_start=min_pos[0];
				int cur_end=min_pos[1];
				if(cur_end-cur_start>=10)
				{
					score[start+min_pos[0]]=0;
					score[start+min_pos[1]]=0;
				}
				else
				{
					for(int k=cur_start;k<=cur_end;k++)score[start+k]=0; 
				}
			}
		}
		//kill small-length fragment
		if(i>0 && i<(int)input.size()-1)
		{
			int prev_end=input[i-1].second;
			int next_start=input[i+1].first;
			if(abs(start-prev_end)>10 && abs(next_start-end)>10)
			{
				if(len<=4)
				{
					//-> reduce the probablity
					for(int k=start;k<=end;k++)score[k]=0;
				}
			}
		}
		if(len<=2)
		{
			//-> reduce the probablity
			for(int k=start;k<=end;k++)score[k]=0;
		}
	}
}



//------- assign label -------//
int Assign_Label(vector < pair<int,int> > &input, int label, vector <int> &output, int &frag, int thres=0)
{
	//proc
	int totnum=0;
	frag=0;
	for(int i=0;i<(int)input.size();i++)
	{
		int start=input[i].first;
		int end=input[i].second;
		if(end-start+1<thres)continue;
		//assign
		for(int j=start;j<=end;j++)output[j]=label;
		totnum+=(end-start+1);
		frag++;
	}
	//return
	return totnum;
}


//-------- load predicted results -------------//
/*
#-> 123
-1 -> 0.909218 0.090782 -> 0.909218  0
-1 -> 0.958281 0.041719 -> 0.958281  0
-1 -> 0.970875 0.029125 -> 0.970875  0
-1 -> 0.974167 0.025833 -> 0.974167  0
...
*/
int Parse_String(string &in,vector <double> &out)
{
	string buf;
	istringstream www(in);
	out.clear();
	for(;;)
	{
		if(!getline(www,buf,' '))break;
		if(buf=="")continue;
		if(buf=="-")break;
		double value=atof(buf.c_str());
		out.push_back(value);
	}
	return (int)out.size();
}
int Load_Predicted_Results(string &input_file,vector <vector <double> > &out_reso,vector <int> &truth_label)
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
	truth_label.clear();
	int first=1;
	int first_num;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		istringstream www(buf);
		getline(www,temp,'>');
		int label=atoi(temp.c_str());
		truth_label.push_back(label);
		getline(www,temp,'>');
		//--- parse ----//
		vector <double> out;
		int retv=Parse_String(temp,out);
		if(first==1)
		{
			first_num=retv;
			first=0;
		}
		else
		{
			if(retv!=first_num)
			{
				fprintf(stderr,"file %s format bad at %s \n",input_file.c_str(),buf.c_str());
				exit(-1);
			}
		}
		//--- push_back ---//
		out_reso.push_back(out);
	}
	//return
	return (int)out_reso.size();
}


//---------- Membrane Protein Labels ----------//
char TM2_To_Char(int in)
{
	switch(in)
	{
		case 0: return '0'; // non-TM
		case 1: return '1'; // TM
		default: return '0';
	}
}
int Return_Label(vector <double> &in,double &maxval)
{
	int i;
	int size=(int)in.size();
	int label=0;
	double wsmax=0;
	for(i=0;i<size;i++)
	{
		if(in[i]>wsmax)
		{
			wsmax=in[i];
			label=i;
		}
	}
	maxval=wsmax;
	return label;
}

//-------- output to TopoPred format ---------//
/*
#TopoPred: transmembrane topology prediction results
#Transmembrane residues are marked with asterisks (*) above threshold 0.500
# Non-transmembrane residues are marked with dots (.) below threshold 0.500

   1 M * 0.988
   2 S * 0.977
   3 Q * 0.970
...
*/
void Output_TopoPred_Format(string &seq_file,string &reso_file,double thres,
	int truth,int cut_len,int estimate_num, int segment_len)
{
	//-> load seq_file
	string sequence;
	int seq_len=Read_FASTA_SEQRES(seq_file,sequence);
	//-> load reso_file
	vector <vector <double> > out_reso;
	vector <int> truth_label;
	int reso_len=Load_Predicted_Results(reso_file,out_reso,truth_label);
	//check length
	if(seq_len!=reso_len)
	{
		fprintf(stderr,"seq_len %d not equal to reso_len %d \n",seq_len,reso_len);
		exit(-1);
	}
	//-> get label
	vector <int> pred_label;
	vector <double> pred_score;
	for(int i=0;i<seq_len;i++)
	{
		//return maximal value
		double maxval;
		int label=Return_Label(out_reso[i],maxval);
		pred_label.push_back(label);
		pred_score.push_back(out_reso[i][1]);
	}
	//-> filter 
	vector < pair<int,int> > first_seg;
	Select_Segment(pred_label,1,first_seg);
	vector < pair<int,int> > second_seg=first_seg;
	for(int i=0;i<seq_len;i++)pred_label[i]=0;
	int fragnum;
	int totnum=Assign_Label(second_seg,1,pred_label,fragnum,cut_len);
	int remove=0;
	if(fragnum<=1 && estimate_num==0)
	{
		remove=1;
		for(int i=0;i<seq_len;i++)pred_label[i]=0;
	}
	//-> cut segment
	Cut_Segment(first_seg,pred_score,segment_len);

	//output
	{
		printf("# PureseqTM: 2-state transmembrane topology prediction results \n");
		printf("# Transmembrane (TM) residues are marked with '1' above threshold %5.3f \n",thres);
		printf("# Non-transmembrane (non-TM) residues are marked with '0' below threshold %5.3f \n",thres);
		if(truth==1)printf("# If ground-truth is given, the TM (non-TM) residues are marked with '1' ('0') at the right-most column \n");
		for(int i=0;i<seq_len;i++)
		{
			int label;
			if(totnum<cut_len || remove==1)
			{
				label=pred_label[i];
			}
			else
			{
				if(pred_score[i]>=thres)label=1;
				else label=0;
			}
			char lab_c=TM2_To_Char(label);
			//output
			if(truth==0)
				printf("%4d %c %c %5.3f \n",i+1,sequence[i],lab_c,out_reso[i][1]);
			else
				printf("%4d %c %c %5.3f %d \n",i+1,sequence[i],lab_c,out_reso[i][1],truth_label[i]);
		}
	}
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//---- TM2_Trans ----//
	{
		if(argc<8)
		{
			fprintf(stderr,"TM2_Trans <seq_file> <reso_file> <threshold> <truth> \n");
			fprintf(stderr,"          <cut_len> <estimate_num> <segment_len> \n");
			fprintf(stderr,"[note]: if truth is 1, then output the 'truth_label' \n");
			fprintf(stderr,"        cut_len shall be set to 10, estimate_num is from phobius \n");
			fprintf(stderr,"        segment_len shall be set to 30-35 \n");
			exit(-1);
		}
		string seq_file=argv[1];
		string reso_file=argv[2];
		double thres=atof(argv[3]);
		int truth=atoi(argv[4]);
		int cut_len=atoi(argv[5]);
		int estimate_num=atoi(argv[6]);
		int segment_len=atoi(argv[7]);
		//process
		Output_TopoPred_Format(seq_file,reso_file,thres,truth,cut_len,estimate_num,segment_len);
		//exit
		exit(0);
	}
}

