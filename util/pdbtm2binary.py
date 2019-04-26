#!/usr/bin/env python
import sys


#------- usage -------#
def Usage():
    print 'python pdbtm2binary.py <input_pdbtm>'


#------- find segment -----#
def find_seg(vec, key, length=10):
    flag = 0
    segment = []
    for i in xrange(len(vec)):
        if (vec[i] == key) and (flag == 0):
            s = i
            flag = 1
        if (vec[i] != key) and (flag == 1):
            e = i - 1
            flag = 0
            if (e-s+1 >= length):
                segment.append((s,e))
        if (i == len(vec) - 1) and (flag == 1):
            if (i-s+1 > length):
                segment.append((s,i))
    return segment


#------- main --------#
def main(argv):
    if len(argv) != 1:
        Usage()
        sys.exit(-1)

    #-> load PDBTM label file
    inputFile = argv[0]
    with open(inputFile) as f:
        pName = f.readline().strip()
        pSeq = f.readline().strip()
        pLabel = f.readline().strip()

    #-> find segment
    segH = find_seg(pLabel, 'H')
    segC = find_seg(pLabel, 'C')
    binaryLabel = ['0'] * len(pSeq)
    for s, e in segH:
        binaryLabel[s: e+1] = ['1'] * (e-s+1)
    for s, e in segC:
        binaryLabel[s: e+1] = ['1'] * (e-s+1)
    binaryLabel = ''.join(binaryLabel)

    #-> output
    print pName
    print pSeq
    print binaryLabel


#-------------- python main ---------------#
if __name__ == "__main__":
    main(sys.argv[1:])

