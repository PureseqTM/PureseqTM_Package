#!/bin/bash


# ----- usage ------ #
usage()
{
	echo "PureseqTM v0.19 [Apr-30-2019] "
	echo "    Predict TransMembrane (TM) topology from a given sequence in FASTA format. "
	echo ""
	echo "USAGE:  ./PureseqTM.sh <-i input_fasta | input_tgt> [-o out_root] [-p signal_pep] [-P plot] [-K remove_tmp] [-H home] "
	echo "                       [-l truth_label] [-s high_thres] [-S low_thres] [-c high_len] [-C low_len] [-L seg_len] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_fasta  : Query protein sequence in FASTA format. "
	echo "(or)"
	echo "-i input_tgt    : Query protein profile file in TGT format"
	echo ""
	echo "***** optional arguments *****"
	echo "-o out_root     : Default output would the current directory. [default = './\${input_name}_PureTM'] "
	echo ""
	echo "-p signal_pep   : Output signal peptide label or not. [default = 0 NOT to output] "
	echo ""
	echo "-P plot         : Plot posterior probabilities or not. [default = 0 NOT to plot] "
	echo ""
	echo "-K remove_tmp   : Remove temporary folder or not. [default = 1 to remove] "
	echo ""
	echo "-H home         : home directory of PureseqTM_Package. "
	echo "                  [default = `dirname $0`] "
	echo "***** additional arguments *****"
	echo "-l truth_label  : The ground-truth 2-state TM label (with 1 for TM and 0 for non-TM) in FASTA format "
	echo "                  [default = null] "
	echo "-s high_thres   : The high threshold to discriminate 'TM' and 'non-TM' regions. [default = 0.5] "
	echo ""
	echo "-S low_thres    : The low threshold to detect 'TM' regions. [default = 0.5] "
	echo ""
	echo "-c high_len     : The high cutoff length to discriminate 'TM' and 'non-TM' regions. [default = 10] "
	echo ""
	echo "-C low_len      : The low cutoff length to detect 'TM' regions. [default = 10] "
	echo ""
	echo "-L seg_len      : The segment len cutoff for 'cut' operation. [default = 30] "
	echo ""
	exit 1
}


#------------------------------------------------------------#
##### ===== get pwd and check BlastSearchHome ====== #########
#------------------------------------------------------------#

#------ current directory ------#
curdir="$(pwd)"

#-------- check usage -------#
if [ $# -lt 1 ];
then
        usage
fi


#---------------------------------------------------------#
##### ===== All arguments are defined here ====== #########
#---------------------------------------------------------#


# ----- get arguments ----- #
#-> required arguments
input=""
input_fasta=""
input_tgt=""
amino_only=0
#-> others
out_root=""
plot=0              #-> default :NOT plot posterior probabilities
signal_pep=0        #-> default: NOT output signal peptide
kill_tmp=1          #-> default: kill temporary root
home=`dirname $0`   #-> default: dirname of the 0-th argument
#-> addi
truth_label=""
high_thres="0.5"
low_thres="0.5"
high_len="10"
low_len="10"
seg_len="30"


#-> parse arguments
while getopts ":i:o:p:P:K:H:l:s:S:c:C:L:" opt;
do
	case $opt in
	#-> required arguments
	i)
		input=$OPTARG
		;;
	#-> optional arguments
	o)
		out_root=$OPTARG
		;;
	p)
		signal_pep=$OPTARG
		;;
	P)
		plot=$OPTARG
		;;
	K)
		kill_tmp=$OPTARG
		;;
	H)
		home=$OPTARG
		;;
	#-> additional arguments
	l)
		truth_label=$OPTARG
		;;
	s)
		high_thres=$OPTARG
		;;
	S)
		low_thres=$OPTARG
		;;
	c)
		high_len=$OPTARG
		;;
	C)
		low_len=$OPTARG
		;;
	L)
		seg_len=$OPTARG
		;;
	#-> others
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done



#---------------------------------------------------------#
##### ===== Part 0: initial argument check ====== #########
#---------------------------------------------------------#

# ------ check home directory ---------- #
if [ ! -d "$home" ]
then
	echo "home directory $home not exist " >&2
	exit 1
fi
home=`readlink -f $home`

# ------ check input fasta ------#
if [ ! -s "$input" ]
then
	echo "input $input not found !!" >&2
	exit 1
fi
input=`readlink -f $input`
fulnam=`basename $input`
relnam=${fulnam%.*}

# ------ determine FASTA or TGT -------#
filename=`basename $input`
extension=${filename##*.}
filename=${filename%.*}
if [ "$extension" == "tgt" ]
then
	input_tgt=$input
	amino_only=0
else
	input_fasta=$input
	amino_only=1
fi

# ------ check output directory ------#
if [ "$out_root" == "" ]
then
	if [ $amino_only -eq 1 ]
	then
		out_root=${relnam}_PureTM
	else
		out_root=${relnam}_ProfTM
	fi
fi
mkdir -p $out_root
out_root=`readlink -f $out_root`



#---------------- initialization ----------------#

# --- default directories --#
bin=$home/bin
util=$home/util
param=$home/param

# --- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp_root="${out_root}/TMP_PureTM_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp_root

# --- truth_label --#
if [ "$truth_label" != "" ]
then
	tail -n1 $truth_label | \
	    awk '{ split($0,chars,""); for (i=1;i<=length($0);i++){printf("%s\n",chars[i])} }' \
	    > $tmp_root/$relnam.label
	truth_label_file=$tmp_root/$relnam.label
	truth_output=1
else
	rm -f $tmp_root/$relnam.label
	touch $tmp_root/$relnam.label
	truth_label_file=$tmp_root/$relnam.label
	truth_output=0
fi

# --- copy initial file to tmp ----#
if [ $amino_only -eq 1 ]
then
	cp $input_fasta $tmp_root/$relnam.seq
else
	cp $input_tgt $tmp_root/$relnam.tgt
	echo ">$relnam" > $tmp_root/$relnam.seq
	head -n4 $tmp_root/$relnam.tgt | tail -n1 | awk '{print $3}' >> $tmp_root/$relnam.seq
fi



#------------------------------------------------------#
##### ===== Part 1: Phobius and Philius ====== #########
#------------------------------------------------------#

#----------- Part 1.1 generate Phobius features ---------------#
params=$param/phobius
runcom="$bin/decodeanhmm -f $params/phobius.options $params/phobius.model"
$runcom $tmp_root/$relnam.seq -plp 1> $tmp_root/$relnam.plp 2> $tmp_root/$relnam.err
grep -v "#" $tmp_root/$relnam.plp | \
            awk '{a=$2;b=$3+$4;c=$5;d=$6+$7+$8+$9;print a" "b" "c" "d}' \
            > $tmp_root/$relnam.phobius
#-> get the number of TM segment from Phobius
$runcom -GFF $tmp_root/$relnam.seq 2> $tmp_root/$relnam.err | grep -v "SEQUENCE" \
             | grep -v "#" > $tmp_root/$relnam.pbseg
rm -f $tmp_root/$relnam.err
tmseg=`awk '{if($3=="M"){print $0}}' $tmp_root/$relnam.pbseg | wc | awk '{print $1}'`
#-> check the existance of N-termi signal peptide from Phobius
ntsp=`awk '{if($3=="h"){print $0}}' $tmp_root/$relnam.pbseg | wc | awk '{print $1}'`

#----------- Part 1.2 generate Philius features ---------------#
params=$param/philius
$util/AA_to_CODE $tmp_root/$relnam.seq $tmp_root/$relnam.obs
echo "$tmp_root/$relnam.obs" > $tmp_root/$relnam.list
$bin/gmtkJT -strFile $params/model_test1.str \
            -triFile $params/model_test1.mytrifile \
            -of1 $tmp_root/$relnam.list -nf1 0 -ni1 2 -fmt1 ascii \
            -inputMasterFile $params/model_test1.master.Viterbi \
            -inputTrainableParameters $params/learnedParamsALL.out \
            -pCliquePrintRange 1:2 \
            -cCliquePrintRange 1:1 \
            -eCliquePrintRange 1:2 \
            -doDist \
            -cptNormThreshold 10. \
            -verbosity 1 > $tmp_root/$relnam.runJT
rm -f jt_info.txt
$util/Proc_runJT $tmp_root/$relnam.runJT $tmp_root/$relnam.philius


#--------------------------------------------------------------#
##### ===== Part 2: AA onehot and AAindex features ====== ######
#--------------------------------------------------------------#

$util/pureseq_dump $tmp_root/$relnam.seq 1 > $tmp_root/$relnam.pureseq
paste $tmp_root/$relnam.phobius $tmp_root/$relnam.philius \
            $tmp_root/$relnam.pureseq > $tmp_root/$relnam.feat_prev
$bin/DeepCNF_FeatMake $tmp_root/$relnam.feat_prev $truth_label_file \
            > $tmp_root/$relnam.feat

#-> has profile or not
if [ $amino_only -eq 0 ]
then
	#-> profile
	$util/profile_dump $tmp_root/$relnam.tgt \
	    > $tmp_root/$relnam.profile
	#-> one-hot
	$util/pureseq_dump $tmp_root/$relnam.seq 0 \
	    > $tmp_root/$relnam.onehot
	#-> AAindex
	awk '{print $1" "$2" "$3" "$4" "$5}' $tmp_root/$relnam.pureseq \
	    > $tmp_root/$relnam.aaindex
	#---- paste all features ----#
	paste $tmp_root/$relnam.profile $tmp_root/$relnam.onehot $tmp_root/$relnam.aaindex \
	    > $tmp_root/$relnam.feat_prof
	#---- feature make ----------#
	$bin/DeepCNF_FeatMake $tmp_root/$relnam.feat_prof $truth_label_file \
            > $tmp_root/$relnam.prof_feat
fi


#--------------------------------------------------------------#
##### ===== Part 3: DeepCNF model to predict TM ====== #########
#--------------------------------------------------------------#

labeldim="2"
featdim="36"
profdim="65"
wind="5,5,5,5,5"
node="100,100,100,100,100"
detect_model=$param/detect.model
predict_model=$param/puretm.model
profile_model=$param/proftm.model
addi=""
#---- discriminate TM and nonTM -----#
$bin/DeepCNF_Pred -i $tmp_root/$relnam.feat -w $wind -d $node -s $labeldim -l $featdim \
            -m $detect_model $addi > $tmp_root/$relnam.detect_out 2> $tmp_root/$relnam.err2
grep -v "#" $tmp_root/$relnam.detect_out | awk '{if($4>thres){print $0}}' thres=$high_thres \
            > $tmp_root/$relnam.detect_reso
rm -f $tmp_root/$relnam.err2
#---- if TM region is detected, then run predict -----#
len=`wc $tmp_root/$relnam.detect_reso | awk '{print $1}'`
if [ $len -gt $high_len ]
then
	if [ $amino_only -eq 1 ]   #-> pureseq mode
	then
		$bin/DeepCNF_Pred -i $tmp_root/$relnam.feat -w $wind -d $node -s $labeldim -l $featdim \
		    -m $predict_model $addi > $tmp_root/$relnam.out 2> $tmp_root/$relnam.err3
	else                       #-> profile mode
		$bin/DeepCNF_Pred -i $tmp_root/$relnam.prof_feat -w $wind -d $node -s $labeldim -l $profdim \
		    -m $profile_model $addi > $tmp_root/$relnam.out 2> $tmp_root/$relnam.err3
	fi
else
	cp $tmp_root/$relnam.detect_out $tmp_root/$relnam.out
fi
rm -f $tmp_root/$relnam.err3

#---- make output readable -----#
$util/Verify_FASTA $tmp_root/$relnam.seq $out_root/$relnam.fasta
$util/TM2_Trans $out_root/$relnam.fasta $tmp_root/$relnam.out $low_thres $truth_output \
            $low_len $tmseg $seg_len > ${out_root}/$relnam.prob
cp $out_root/$relnam.fasta ${out_root}/$relnam.top
grep -v "#" ${out_root}/$relnam.prob | awk '{printf $3}END{printf "\n"}' \
            >> ${out_root}/$relnam.top



#-------------------------------------------------------------#
##### ===== Part 4: Check Signal Peptide or TM ====== #########
#-------------------------------------------------------------#

#-> process signal peptide
if [ $ntsp -gt 0 ] && [ $signal_pep -eq 1 ]
then
	#--- merge SP and TM ---#
	$util/Merge_SP_TM $tmp_root/$relnam.plp ${out_root}/$relnam.top 2 \
	    > $tmp_root/$relnam.top_sptm
	#--- final output ------#
	cp $out_root/$relnam.fasta ${out_root}/$relnam.topsp
	cat $tmp_root/$relnam.top_sptm >> ${out_root}/$relnam.topsp
fi


#-------------------------------------------------------------#
##### ===== Part 5: Plot Posterior Probabilities ====== #######
#-------------------------------------------------------------#

#-> plot posterior probabilities
if [ $plot -eq 1 ]
then
	$util/GnuPlot_Script ${out_root}/$relnam.prob $low_thres $truth_output \
	    ${out_root}/$relnam.png $tmp_root/$relnam.gff > $tmp_root/$relnam.gnuplot_script
	gnuplot $tmp_root/$relnam.gnuplot_script
	sort -n -k 3 $tmp_root/$relnam.gff > $out_root/$relnam.gff
	rm -f $tmp_root/$relnam.gff
fi


#----------- remove tmp -------------#
cp $tmp_root/$relnam.seq $out_root/$relnam.fasta_raw
echo "$amino_only" > $out_root/$relnam.pred_mode
if [ $kill_tmp -eq 1 ]
then
	rm -rf $tmp_root
else
	rm -rf $out_root/"TMP_PureTM_"${relnam}
	mv $tmp_root $out_root/"TMP_PureTM_"${relnam}
fi

#=============== exit ===============#
exit


