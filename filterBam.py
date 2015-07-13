#!/usr/bin/env python3
#script to discard reads without sufficient quality information
#Julin Maloof
#December 1, 2013

import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument("file_in", help="bam file")
parser.add_argument("file_out", help="output file")

args = parser.parse_args()

file_in = pysam.Samfile(args.file_in,"rb")
file_out = pysam.Samfile(args.file_out,"wb",template=file_in)

i = 0
discard = 0

for read in file_in.fetch(until_eof=True):
    i += 1
    if i % 1000000 == 0:
        print("processed ",i," reads")
    try:
        qual=read.qual
    except ValueError: #occurs if quality length not the same as readlength
        discard +=1
        print("discarding read because of ValueError")
        print(read)
        continue
    try:
        if(len(qual) != read.rlen):
            discard += 1
            print("discarding read because lengths don't match",len(qual),read.rlen)
            print(read)
        else:
            file_out.write(read)
    except TypeError:
        discard += 1
        print("discarding read because can't obtain lengths")

file_in.close()
file_out.close()

print("Done. Processed ",i," reads, discarded ",discard," reads.")


