#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <iomanip>
#include <time.h>
#include <algorithm>
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

//=================== upper and lower case ====================//
//----------upper_case-----------//
void toUpperCase(char *buffer) 
{  
	for(int i=0;i<(int)strlen(buffer);i++) 
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
void toUpperCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
//----------lower_case-----------//
void toLowerCase(char *buffer)
{  
	for(int i=0;i<(int)strlen(buffer);i++) 
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}
void toLowerCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}


//========== Part 0: Data Preparation ==========//
//-> input examle
/*
 1 0.842039
 1 0.775209
 1 0.755237
 1 0.661877
 1 0.515135
 1 0.394574
 1 0.357573
 1 0.277307
 1 0.222007
 1 0.168835
 1 0.149776
 1 0.136482
 1 0.129423
 1 0.115352
 1 0.115104
 1 0.121574
 1 0.143086
 1 0.161823
 1 0.184361
 1 0.196237
 1 0.203620
 1 0.212567
 0 0.209723
 0 0.185400
 0 0.181789
 0 0.179994
*/
int Load_Data(string &in_file,vector <int> &label, vector <double> &value)
{
	ifstream fin;
	string buf,temp;
	fin.open(in_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"%s not found!\n",in_file.c_str());
		exit(-1);
	}
	//load
	label.clear();
	value.clear();
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		istringstream www(buf);
		int lab;
		double val;
		www>>lab>>val;
		label.push_back(lab);
		value.push_back(val);
		count++;
	}
	return count;
}

//-------- Get prediction by thres -----//
//-> for binary classification, we mark '1' when the prob > thres
void Get_Prediction(vector <double> &value, double thres, vector <int> &pred_reso)
{
	int i;
	int size=(int)value.size();
	pred_reso.resize(size,0);
	for(i=0;i<size;i++)
	{
		int pred;
		if(value[i]<thres)pred=0;
		else pred=1;
		pred_reso[i]=pred;
	}
}

//--------------- select segment ----------//
void Select_Segment(vector <int> &input, int label,vector < pair<int,int> > & output)
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
				output.push_back(pair<int,int>(start,end));
			}
		}
	}
	if(first==0)
	{
		first=1;
		end=(int)input.size()-1;
		output.push_back(pair<int,int>(start,end));
	}
}

//------ calculate overlap betwee two segments ---------//
int Check_Overlap(int s1_start,int s1_end,int s2_start,int s2_end)
{
	if(s1_start<=s2_end && s1_end>=s2_start)
	{
		//-> prepare two ranges
		vector <int> first_range;
		for(int i=s1_start;i<=s1_end;i++)first_range.push_back(i);
		vector <int> second_range;
		for(int i=s2_start;i<=s2_end;i++)second_range.push_back(i);
		//-> calculate overlap number
		vector <int> common_data;
		set_intersection(first_range.begin(), first_range.end(), 
			second_range.begin(), second_range.end(), 
			std::back_inserter(common_data) );
		//-> return overlap number
		return (int)common_data.size();
	}
	else return 0;
}


//=============== Part 1: calculate binary SOV_score =================//

//--------- overlap for SOV_score -------//
double Calculate_Overlap(int s1_start,int s1_end,int s2_start,int s2_end)
{
	int start1=s1_start<s2_start?s1_start:s2_start;
	int end1=s1_end>s2_end?s1_end:s2_end;
	int start2=s1_start>s2_start?s1_start:s2_start;
	int end2=s1_end<s2_end?s1_end:s2_end;
	int max_s1_s2=end1-start1+1;
	int min_s1_s2=end2-start2+1;
	int l1=(int)((s1_end-s1_start+1)/2);
	int l2=(int)((s2_end-s2_start+1)/2);
	int l1_rel=s1_end-s1_start+1;
	int min1=(max_s1_s2-min_s1_s2)<min_s1_s2?(max_s1_s2-min_s1_s2):min_s1_s2;
	int min2=l1<l2?l1:l2;
	int delta_s1_s2=min1<min2?min1:min2;
	return 1.0*l1_rel*(min_s1_s2+delta_s1_s2)/max_s1_s2;
}

//---------- SOV_score for a given label -------------//
double SOV_score_single(vector < pair<int,int> > & seg1, vector < pair<int,int> > & seg2, int &n_len)
{
	double sov_score=0;
	n_len=0;
	for(int i=0;i<(int)seg1.size();i++)
	{
		int s1_start=seg1[i].first;
		int s1_end=seg1[i].second;
		int hit=0;
		for(int j=0;j<(int)seg2.size();j++)
		{
			int s2_start=seg2[j].first;
			int s2_end=seg2[j].second;
			int check=Check_Overlap(s1_start,s1_end,s2_start,s2_end);
			if(check==0)continue;
			double score=Calculate_Overlap(s1_start,s1_end,s2_start,s2_end);
			sov_score+=score;
			n_len+=s1_end-s1_start+1;
			hit=1;
		}
		//not hit
		if(hit==0)n_len+=s1_end-s1_start+1;
	}
	return sov_score;
}

//-------- calculate Q2 SOV score -----//
//-> by default, we use FIX mode according to the following paper:
/*
A Modified Definition of Sov, a Segment-Based Measure for Protein Secondary Structure Prediction Assessment
Adam Zemla, Ceslovas Venclovas, Krzysztof Fidelis, and Burkhard Rost
PROTEINS: Structure, Function, and Genetics 34:220?23 (1999)
*/
double Calculate_SOV_Score(vector <int> &seq1,vector <int> &seq2,int ORIorFIX=1)
{
	//calculate overall SOV score
	double overall_score=0;
	int overall_nlen=0;
	for(int k=0;k<=1;k++)
	{
		//-> extract fragment according to the given label
		vector < pair<int,int> > seg1;
		vector < pair<int,int> > seg2;
		Select_Segment(seq1,k,seg1);
		Select_Segment(seq2,k,seg2);
		//-> calculate single SOV score for the given label
		int n_len;
		overall_score+=SOV_score_single(seg1, seg2,n_len);
		overall_nlen+=n_len;
	}
	//output score
	int len;
	if(ORIorFIX==0)len=(int)seq1.size();
	else len=overall_nlen;
	return 1.0*overall_score/len;
}



//========== Part 2: calculate TP/FP related measurements ==========//

//-------- calculate TP/FP -----------//
void Calculate_TP_FP_Value(vector <int> &label, vector <int> &pred_reso,
	int &TP, int &FP, int &TN, int &FN)
{
	int i;
	int size=(int)label.size();
	TP=0;
	FP=0;
	TN=0;
	FN=0;
	for(i=0;i<size;i++)
	{
		if(pred_reso[i]==label[i])
		{
			if(pred_reso[i]==1)TP++;
			else TN++;
		}
		else
		{
			if(pred_reso[i]==1)FP++;
			else FN++;
		}
	}
}

//---- predictive -----//
//-> precision (or, positive predictive value, PPV)
double Precision(int TP, int FP, int TN, int FN)
{
	if(TP==0)return 0;
	return 1.0*TP/(TP+FP);
}
//-> negative_predictive_value ( NPV)
double negative_predictive_value(int TP, int FP, int TN, int FN)
{
	if(TN==0)return 0;
	return 1.0*TN/(TN+FN);
}

//----- data -----//
//-> recall (or, sensitivity, true positive rate, TPR)
double Recall(int TP, int FP, int TN, int FN)
{
	if(TP==0)return 0;
	return 1.0*TP/(TP+FN);
}
//-> specificity (or, true negative rate, TNR)
double specificity(int TP, int FP, int TN, int FN)
{
	if(TN==0)return 0;
	return 1.0*TN/(TN+FP);
}

//----- accuracy ----//
//-> accuracy
double Accuracy(int TP, int FP, int TN, int FN)
{
	if(TP+TN==0)return 0;
	return 1.0*(TP+TN)/(TP+TN+FP+FN);
}
//-> balanced accurcy
double balanced_accurcy(int TP, int FP, int TN, int FN)
{
	if(TP==0 && TN==0)return 0;
	if(TP==0)return 0.5*TN/(TN+FP);
	if(TN==0)return 0.5*TP/(TP+FN);
	return 0.5*( 1.0*TP/(TP+FN) + 1.0*TN/(TN+FP) );
}
//-> F1 score
double F1_score(int TP, int FP, int TN, int FN)
{
	if(TP==0)return 0;
	return 2.0*TP/(2.0*TP+FP+FN);
}

//------ Matthews correlation coefficient (MCC) ---//
double MCC_Value(int TP, int FP, int TN, int FN)
{
	if(TP*TN==0 && FP*FN==0)return 0;
	return 1.0*( 1.0*TP*TN-1.0*FP*FN )/sqrt( 1.0*(TP+FP)*(TN+FP)*(TP+FN)*(TN+FN) );
}


//========== Part 3: calculate protein-level and segment-level measurement ==========//

//---- protein-level accuracy -----//
//-> we only consider segment marked with '1'
int Protein_Accuracy_Strict(vector <int> &seq1,vector <int> &seq2, int mark=1,int len=5)
{
	vector < pair<int,int> > seg1;
	vector < pair<int,int> > seg2;
	Select_Segment(seq1,mark,seg1);
	Select_Segment(seq2,mark,seg2);
	//-- calculate --//
	if(seg1.size() != seg2.size())return 0;
	for(int i=0;i<(int)seg1.size();i++)
	{
		//-> overlap_1
		int s1_start=seg1[i].first;
		int s1_end=seg1[i].second;
		//-> overlap_2
		int s2_start=seg2[i].first;
		int s2_end=seg2[i].second;
		//-> check overlap
		int check=Check_Overlap(s1_start,s1_end,s2_start,s2_end);
		if(check<len)return 0;
	}
	//return
	return 1;
}
//-> here we only consider TM and non-TM
//[note]: we ALWAYS put 'ground-truth' to the first, and 'prediction' to the second
int Protein_Accuracy_Loose(vector <int> &seq1,vector <int> &seq2, int mark=1,int len=5)
{
	vector < pair<int,int> > seg1;
	vector < pair<int,int> > seg2;
	Select_Segment(seq1,mark,seg1);
	Select_Segment(seq2,mark,seg2);
	//-- calculate --//
	if(seg1.size()==0) //-> non-TM
	{
		if(seg2.size()==0)return 1;
		else return 0;
	}
	else               //-> TM
	{
		if(seg2.size()!=0)return 1;
		else return 0;
	}
}

//---- segment-level acccuracy -----//
//-> we only consider segment marked with '1'
double Segment_Accuracy(vector <int> &seq1,vector <int> &seq2, int mark=1,int len=5)
{
	vector < pair<int,int> > seg1;
	vector < pair<int,int> > seg2;
	Select_Segment(seq1,mark,seg1);
	Select_Segment(seq2,mark,seg2);
	//-- calculate --//
	int correct=0;
	int total=0;
	for(int i=0;i<(int)seg1.size();i++)
	{
		//-> overlap_1
		int s1_start=seg1[i].first;
		int s1_end=seg1[i].second;
		for(int j=0;j<(int)seg2.size();j++)
		{
			//-> overlap_2
			int s2_start=seg2[j].first;
			int s2_end=seg2[j].second;
			//-> check overlap
			int check=Check_Overlap(s1_start,s1_end,s2_start,s2_end);
			if(check>=len)
			{
				correct++;
				break;
			}
		}
		//-> total
		total++;
	}
	//return
	if(total==0)return 0;
	else return 1.0*correct/total;
}


//================= main function =======================//
void TM2_Evaluation(string &label_value_file, double thres)
{
	//---- load data -----//
	vector <int> label;
	vector <double> value;
	vector <int> pred_reso;
	Load_Data(label_value_file,label,value);
	Get_Prediction(value,thres,pred_reso);

	//---- Part 1: Q2_SOV_Score ----//
	double sov_score=Calculate_SOV_Score(label,pred_reso);

	//---- Part 2: TP/FP relevant score ----//
	//-> calculate TP/FP	
	int TP,FP,TN,FN;
	Calculate_TP_FP_Value(label,pred_reso,TP,FP,TN,FN);
	//-> calculate TP/FP related values
	double Prec=Precision(TP, FP, TN, FN);
	double NPV=negative_predictive_value(TP, FP, TN, FN);
	double Reca=Recall(TP, FP, TN, FN);
	double Spec=specificity(TP, FP, TN, FN);
	double Acc=Accuracy(TP, FP, TN, FN);
	double Bacc=balanced_accurcy(TP, FP, TN, FN);
	double F1=F1_score(TP, FP, TN, FN);
	double MCC=MCC_Value(TP, FP, TN, FN);

	//---- Part 3: protein-level and segment-level score ----//
	//-> protein-level score
	int strict_acc=Protein_Accuracy_Strict(label,pred_reso);
	int loose_acc=Protein_Accuracy_Loose(label,pred_reso);
	//-> segment-level score
	double segment_recall=Segment_Accuracy(label,pred_reso);
	double segment_precision=Segment_Accuracy(pred_reso,label);

	//-> printf
	string name;
	getBaseName(label_value_file,name,'/','.');
	printf("%s : pAccS %d pAccL %d | sReca %f sPrec %f | SOV %f | TP %d FP %d TN %d FN %d -> Prec %f NPV %f Reca %f Spec %f Acc %f Bacc %f F1 %f MCC %f \n",
		name.c_str(),strict_acc,loose_acc,segment_recall,segment_precision,sov_score,TP,FP,TN,FN,Prec,NPV,Reca,Spec,Acc,Bacc,F1,MCC);
}



//----------- main -------------//
int main(int argc,char **argv)
{

	//---- TM2_Evaluation ----//__190404__//
	{
		if(argc<3)
		{
			fprintf(stderr,"Version 1.01 \n");
			fprintf(stderr,"TM2_Evaluation <label_value_file> <threshold> \n");
			fprintf(stderr,"[note]: <label_value_file> should be <label, value> \n");
			fprintf(stderr,"        <label> should be 0 or 1, with 1 be true label \n");
			fprintf(stderr,"        <value> between 0 to 1, more higher more true label \n");
			fprintf(stderr,"        for example:    1   0.923 \n");
			exit(-1);
		}
		//---- read argument ----//
		string label_value_file=argv[1];
		double threshold=atof(argv[2]);
		//---- process -----//
		TM2_Evaluation(label_value_file,threshold);
		//exit
		exit(0);		
	}
}

