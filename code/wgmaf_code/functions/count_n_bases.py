import sys
import pysam

bamFile = sys.argv[1]
fragmentFile = sys.argv[2]
outFile = sys.argv[3]

# Getting the list of fragments to get the number of 'N' bases in
# and converting to a set for faster lookup
ff = open(fragmentFile, "r")
fragment_list = ff.read().splitlines()
ff.close()
fragment_set = set(fragment_list)

# Opening the bam file
alnFile = pysam.AlignmentFile(bamFile, "rb")

# Opening the output file to write to
f = open(outFile, "w")
f.write("fragment" + "\t" + "read" + "\t" + "n_bases" + "\n")

iter = alnFile.fetch()

# Iterating through each read
for i in iter:

  # Get fragment name
  qname = i.qname
 
  # If qname is in fragment_set, get number of 'N' bases
  if qname in fragment_set:
    # Get read 1/2 status
    if i.is_read1 == True:
      read = "1"
    elif i.is_read2 == True:
      read = "2"
    else:
      print("Error: read must be 1 or 2.")
      exit()

    # Get number of 'N' bases in read
    n_bases = i.query.count("N")  

    f.write(qname + "\t" + read + "\t" + str(n_bases) + "\n") 

f.close()
alnFile.close()
