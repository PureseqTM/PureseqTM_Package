# PureseqTM
TransMembrane (TM) topology prediction from amino acid sequence only

<p align="center">
<img src="https://github.com/PureseqTM/PureseqTM_Package/blob/master/example/pipeline_v0.19.png"/>
</p>

## Server
Users may submit their sequences to the PureseqTM Server at: http://pureseqtm.predmp.com

### Train and test data
Users may find the train, test data of PureseqTM at: https://github.com/PureseqTM/PureseqTM_Dataset 

### Paper
https://www.biorxiv.org/content/10.1101/627307v2.abstract

### Contact email
pureseqtm@gmail.com


# Install
### Download the PureseqTM package
```
git clone https://github.com/PureseqTM/PureseqTM_Package.git
cd ./PureseqTM_Package
```

### Compilation if necessary
```
cd source_code
    make
cd ../
```

# Examples

## single protein prediction

#### prediction with fast mode (i.e., NOT use DBN features from Philius)
```
./PureseqTM.sh -i example/4j7cK.fasta -m 0
```

#### prediction without ground-truth label
```
./PureseqTM.sh -i example/4j7cK.fasta
```

#### given ground-truth 2-state TM label
```
./PureseqTM.sh -i example/4j7cK.fasta -l example/4j7cK.top
```

## evaluate the prediction result
If the ground-truth label is provided, run below commands to evaluate the prediction accuracy:
```
grep -v "#" 4j7cK_PureTM/4j7cK.prob | awk '{print $NF" "$3}' > 4j7cK.pred_reso
util/TM2_Evaluation 4j7cK.pred_reso 0.5
rm -f 4j7cK.pred_reso
```

## plot the posterior probabilities
```
./PureseqTM.sh -i example/4j7cK.fasta -P 1
```
Note that gnuplot is required to be installed in the local system.

## prediction with given profile
```
./PureseqTM.sh -i example/1bhaA.tgt
```
Note that to generate the evolutionary profile, users are suggested to use the [TGT_Package](https://github.com/realbigws/TGT_Package) to generate the TGT file:
```
./A3M_TGT_Gen.sh -i <input_fasta> -d uniprot20_2016_02 -h hhsuite2 -n 3 -E 10
```
For more details on evolutionary profile, please goto https://github.com/realbigws/TGT_Package

## proteome prediction
As it takes about 2 hours for the whole Human proteome, below please find a toy example:
```
./PureseqTM_proteome.sh -i example/test_proteome.fasta
```

## transfer PDBTM label
Transfer PDBTM 9-state label to 0/1 TransMembrane (TM) label:
```
python util/pdbtm2binary.py example/1bhaA.pdbtm
```


# Output files
There shall be 5 to 7 output files in XXX_PureTM by default, where XXX is the input name:

| File name     | Description   | Option |
| ------------- | ------------- | ------ |
| XXX.fasta_raw | original input sequence file in FASTA format. | |
| XXX.fasta     | canonical sequence file without non-canonical characters. | |
| XXX.top       | simple 2-state TransMembrane (TM) prediction in FASTA format. | |
| XXX.prob      | detailed 2-state TransMembrane (TM) probability prediction. | |
| XXX.pred_mode | prediction mode (1 for sequence and 0 for profile). | |
| XXX.png       | posterior probabilities plotted by GNUPLOT | if option -P 1 is set |
| XXX.gff       | segment-level output | if option -P 1 is set |


# References

### [Phobius](http://phobius.sbc.su.se/)
```
Title:
     A Combined Transmembrane Topology and Signal Peptide Prediction Method
Authors:
     Lukas Kall, Anders Krogh and Erik L. L. Sonnhammer
Journal:
     J. Mol. Biol. (2004) 338, 1027–1036
```

### [Philius](http://www.yeastrc.org/philius/pages/philius/runPhilius.jsp)
```
Title:
     Transmembrane Topology and Signal Peptide Prediction Using Dynamic Bayesian Networks
Authors:
     Sheila M. Reynolds, Lukas Kall, Michael E. Riffle, Jeff A. Bilmes, William Stafford Noble
Journal:
     PLoS computational biology (2008) 4(11):e1000213
```

### [DeepCNF](https://github.com/realbigws/DeepCNF_AUC)
```
Title:
     Protein secondary structure prediction using deep convolutional neural fields
Authors:
     Sheng Wang, Jian Peng, Jianzhu Ma, Jinbo Xu
Journal:
     Scientific reports (2016) 6:18962
```

