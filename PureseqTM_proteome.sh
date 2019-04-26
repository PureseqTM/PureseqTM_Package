#!/bin/bash

# ----- usage ------ #
usage()
{
	echo "PureseqTM_proteome v0.10 [Apr-04-2019] "
	echo "    Predict TransMembrane (TM) topology from a given proteome in FASTA format. "
	echo ""
	echo "USAGE:  ./PureseqTM_proteome.sh <-i in_proteome> [-o out_file] [-K remove_tmp] "
	echo "                                [-c CPU_num] [-H home] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i in_proteome  : Query proteome in FASTA format. "
	echo ""
	echo "***** optional arguments *****"
	echo "-o out_file     : Default output file would be './\${input_name}_topo' "
	echo ""
	echo "-K remove_tmp   : Remove temporary folder or not. [default = 1 to remove] "
	echo ""
	echo "-c CPU_num      : Number of processors. [default = 4] "
	echo ""
	echo "-H home         : home directory of PureseqTM_proteome."
	echo "                  [default = `dirname $0`] "
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
in_proteome=""
#-> others
out_file=""         #-> default: '${input_name}_topo'
kill_tmp=1          #-> default: kill temporary root
cpu_num=4           #-> default : use 4 CPUs
home=`dirname $0`   #-> default: dirname of the 0-th argument


#-> parse arguments
while getopts ":i:o:K:c:H:" opt;
do
	case $opt in
	#-> required arguments
	i)
		in_proteome=$OPTARG
		;;
	#-> optional arguments
	o)
		out_file=$OPTARG
		;;
	K)
		kill_tmp=$OPTARG
		;;
	c)
		cpu_num=$OPTARG
		;;
	H)
		home=$OPTARG
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
if [ ! -s "$in_proteome" ]
then
	echo "in_proteome $in_proteome not found !!" >&2
	exit 1
fi
in_proteome=`readlink -f $in_proteome`
fulnam=`basename $in_proteome`
relnam=${fulnam%.*}

# ------ check output file ------#
if [ "$out_file" == "" ]
then
	out_file=${relnam}_topo
fi
out_file=`readlink -f $out_file`



#---------------- initialization ----------------#

# --- default directories --#
util=$home/util
bin=$home/bin

# --- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp_root="/tmp/PureTM_Proteome_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp_root


#------------------------------------------------#
##### ===== Part 1: Run PureseqTM ====== #########
#------------------------------------------------#

#------ extract each sequences ----#
mkdir -p $tmp_root/${relnam}_fasta
$util/MSA_To_FASTA_num $in_proteome $tmp_root/${relnam}_fasta
ls $tmp_root/${relnam}_fasta | cut -d '.' -f 1 > $tmp_root/${relnam}_list

#------ batch process -----#
mkdir -p $tmp_root/${relnam}_puretm
rm -f $tmp_root/${relnam}_puretm_proc
for i in `cat $tmp_root/${relnam}_list`
do
	echo "$home/PureseqTM.sh -i $tmp_root/${relnam}_fasta/$i.fasta -o $tmp_root/${relnam}_puretm/${i} -K $kill_tmp -H $home" >> $tmp_root/${relnam}_puretm_proc
done
$bin/distribute_ii_openmp_cpu $tmp_root/${relnam}_puretm_proc $cpu_num 

#------ collect result ----#
rm -f $out_file
for i in `cat $tmp_root/${relnam}_list`
do
	cat $tmp_root/${relnam}_puretm/${i}/$i.fasta_raw >> $out_file
	tail -n1 $tmp_root/${relnam}_puretm/${i}/$i.top >> $out_file
done


#----------- remove tmp -------------#
if [ $kill_tmp -eq 1 ]
then
	rm -rf $tmp_root
else
	mv $tmp_root ${relnam}.PureTM_Proteome
fi


#-------------- exit ---------------#
exit 0



