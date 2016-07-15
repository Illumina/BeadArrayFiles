from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype
import sys
import os
from datetime import datetime

delim = "\t"

if len(sys.argv) < 4:
	sys.stderr.write("Generate a final report from a directory of GTC files\n")
	sys.stderr.write("usage: python gtc_final_report.py <BPM manifest file> <GTC directory> <output file>\n")
	sys.exit(-1)

try:
	names = BeadPoolManifest(sys.argv[1]).names
except:
	sys.stderr.write("Failed to read loci names from manifest\n")
	sys.exit(-1)

output_file = sys.argv[3]

if os.path.isfile( output_file ):
	sys.stderr.write("Output file already exists, please delete and re-run\n")
	sys.exit(-1)

with open(output_file, "w") as output_handle:
	output_handle.write("[Header]\n")
	output_handle.write(delim.join(["Processing Date", datetime.now().strftime("%m/%d/%Y %I:%M %p") ] )+ "\n")
	output_handle.write(delim.join(["Content", os.path.basename( sys.argv[1]) ]) + "\n")
	output_handle.write(delim.join(["Num SNPs", str( len(names) )]) + "\n")
	output_handle.write(delim.join(["Total SNPs", str( len(names))]) + "\n")

	samples = []
	for file in os.listdir(sys.argv[2]):
		if file.lower().endswith(".gtc"):
			samples.append( file )

	output_handle.write(delim.join(["Num Samples", str(len(samples))]) + "\n")
	output_handle.write(delim.join(["Total Samples", str(len(samples))]) + "\n")

	output_handle.write("[Data]\n")
	output_handle.write(delim.join([ "SNP Name", "Sample ID", "Alleles - AB"]) + "\n")
	for file in samples:
		sys.stderr.write("Processing " + file + "\n")
		gtc_file = os.path.join( sys.argv[2], file )
		genotypes = GenotypeCalls(gtc_file).get_genotypes()
		assert len(genotypes) == len(names)
		for (name, genotype) in zip( names, genotypes):
			output_handle.write(delim.join( [name, file[:-4],code2genotype[genotype], ]  ) + "\n")
