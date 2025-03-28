import sys
import pysam

bamFile = sys.argv[1]
regionFile = sys.argv[2]
outFile = sys.argv[3]

def getReadPileup(read, region):

  if read == 1:
    flag_filter = 1152
  elif read == 2:
    flag_filter = 1088
  else:
    print("Error: read must be 1 or 2.")
    exit()

  iter = alnFile.pileup(region = region,
                        max_depth = 1e9,
                        ignore_overlaps = False,
                        min_base_quality = 0,
                        min_mapping_quality = 0,
                        flag_require = 2,
                        flag_filter = flag_filter,
                        truncate = True)

  for pu in iter:

    qnames = pu.get_query_names()
    chrom = pu.reference_name
    pos = pu.reference_pos + 1 # Convert to 1-based coordinate system
    obs = pu.get_query_sequences(add_indels = False) # Do not include insertions in the pileup, deletions are an empty string
    mapq = pu.get_mapping_qualities()
    phred = pu.get_query_qualities()

    for i in range(pu.get_num_aligned()):
      obs_i = obs[i]
      mapq_i = mapq[i]
      phred_i = phred[i]
      if obs_i != "": # Exclude deletions
        if obs_i.islower():
          strand = "-"
          obs_i = obs_i.upper()
        elif obs_i.isupper():
          strand = "+"
        else:
          print("Error: the observed base is neither uppercase nor lowercase..")
          exit()

        if obs_i != "N":
          f.write(qnames[i] + "\t" + str(read) + "\t" + strand + "\t" + chrom + "\t" + str(pos) + "\t" + obs_i + "\t" + str(mapq_i) + "\t" + str(phred_i) + "\n")


# Opening the bam file
alnFile = pysam.AlignmentFile(bamFile, "rb")

# Opening the output file to write to
f = open(outFile, "w")
f.write("fragment" + "\t" + "read" + "\t" + "strand" + "\t" + "chrom" + "\t" + "pos" + "\t" + "obs" + "\t" + "mapq" + "\t" + "phred" + "\n")

# Getting the list of regions to perform the pileup on
rf = open(regionFile, "r")
region_list = rf.read().splitlines()
rf.close()

# For each region, performing a pileup for read 1 and read 2 separately
for region in region_list:
  getReadPileup(read = 1, region = region)
  getReadPileup(read = 2, region = region)

f.close()
alnFile.close()








