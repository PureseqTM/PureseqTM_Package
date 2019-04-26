

# make the subdir first
SHELL :=/bin/bash
CXX = g++
FLAG = -O3

all:
	#----- for Philius -----#
	$(CXX) AA_to_CODE.cpp $(FLAG) -o AA_to_CODE
	$(CXX) Proc_runJT.cpp $(FLAG) -o Proc_runJT
	#----- for DeepCNF -----#
	cd profile_dump_src&&make
	$(CXX) pureseq_dump.cpp $(FLAG) -o pureseq_dump
	#----- for evaluation --#
	$(CXX) Verify_FASTA.cpp $(FLAG) -o Verify_FASTA
	$(CXX) TM2_Trans.cpp $(FLAG) -o TM2_Trans
	$(CXX) TM2_Evaluation.cpp $(FLAG) -o TM2_Evaluation
	#----- for misc --------#
	$(CXX) MSA_To_FASTA_num.cpp $(FLAG) -o MSA_To_FASTA_num
	$(CXX) Merge_SP_TM.cpp $(FLAG) -o Merge_SP_TM
	$(CXX) GnuPlot_Script.cpp $(FLAG) -o GnuPlot_Script
	#-> for Philius
	mv AA_to_CODE Proc_runJT ../util
	#-> for DeepCNF
	mv profile_dump pureseq_dump ../util
	#-> for evaluation
	mv Verify_FASTA TM2_Trans TM2_Evaluation ../util
	#-> for misc
	mv MSA_To_FASTA_num Merge_SP_TM GnuPlot_Script ../util

# clean all the temp files
clean:

