import sys
import argparse
import os
from datetime import datetime
from IlluminaBeadArrayFiles import BeadPoolManifest, GenotypeCalls, code2genotype
from pandas import DataFrame
from numpy import around

"""
This example will allow the user to create matrix report files identical to the ones exported from GenomeStudio
by chosing the genotype with or without the genotype scores.
"""

def build_dict(matrix = None, samplename = None, names = None, genotypes = None, genotype_scores = None):
    if genotype_scores is not None:
        for (name, genotype, genotype_score) in zip(names, genotypes, genotype_scores):
            if samplename in matrix:
                matrix[samplename][name] = genotype + "|" + str(genotype_score)
            else:
                matrix[samplename] = {}
                matrix[samplename][name] = genotype + "|" + str(genotype_score)
    else:
        for (name, genotype) in zip(names, genotypes):
            if samplename in matrix:
                matrix[samplename][name] = genotype
            else:
                matrix[samplename] = {}
                matrix[samplename][name] = genotype
    return matrix

delim = "\t"
NUM_ARGUMENTS = 6

parser = argparse.ArgumentParser("Generate a final report from a directory of GTC files")
parser.add_argument("manifest", help="BPM manifest file")
parser.add_argument("gtc_directory", help="Directory containing GTC files")
parser.add_argument("output_file", help="Location to write report")
parser.add_argument("--forward", help="python gtc_final_report_matrix.py <path_to_manifest> <path_to_gtc_directory> <path_to_output_file> --forward 1, print matrix with forward alleles")
parser.add_argument("--forward_GC", help="python gtc_final_report_matrix.py <path_to_manifest> <path_to_gtc_directory> <path_to_output_file> --forward_GC 1, print matrix with forward alleles including genotype scores")
parser.add_argument("--top", help="python gtc_final_report_matrix.py <path_to_manifest> <path_to_gtc_directory> <path_to_output_file> --top 1, print matrix with top alleles")
parser.add_argument("--top_GC", help="python gtc_final_report_matrix.py <path_to_manifest> <path_to_gtc_directory> <path_to_output_file> --top_GC 1, print matrix with top alleles including genotype scores")
parser.add_argument("--AB", help="python gtc_final_report_matrix.py <path_to_manifest> <path_to_gtc_directory> <path_to_output_file> --forward 1, print matrix with forward alleles")
parser.add_argument("--AB_GC", help="python gtc_final_report_matrix.py <path_to_manifest> <path_to_gtc_directory> <path_to_output_file> --forward_GC 1, print matrix with forward alleles including genotype scores")
parser.add_argument("--plus", help="python gtc_final_report_matrix.py <path_to_manifest> <path_to_gtc_directory> <path_to_output_file> --top 1, print matrix with top alleles")
parser.add_argument("--plus_GC", help="python gtc_final_report_matrix.py <path_to_manifest> <path_to_gtc_directory> <path_to_output_file> --top_GC 1, print matrix with top alleles including genotype scores")

args = parser.parse_args()

if len(sys.argv) != NUM_ARGUMENTS:
    sys.stderr.write("For matrix report user needs to provide either forward or top strand parameter with or without genotype score, can only build one report!\n")
    sys.exit(-1)

if os.path.isfile(args.output_file):
    sys.stderr.write("Output file already exists, please delete and re-run\n")
    sys.exit(-1)

try:
    manifest = BeadPoolManifest.BeadPoolManifest(args.manifest)
except:
    sys.stderr.write("Failed to read data from manifest\n")
    sys.exit(-1)

with open(args.output_file, "w") as output_handle:
    output_handle.write("[Header]\n")
    output_handle.write(delim.join(["Processing Date", datetime.now().strftime("%m/%d/%Y %I:%M %p")])+ "\n")
    output_handle.write(delim.join(["Content", os.path.basename(args.manifest)]) + "\n")
    output_handle.write(delim.join(["Num SNPs", str(len(manifest.names))]) + "\n")
    output_handle.write(delim.join(["Total SNPs", str(len(manifest.names))]) + "\n")
    
    samples = []
    for gtc_file in os.listdir(args.gtc_directory):
        if gtc_file.lower().endswith(".gtc"):
            samples.append(gtc_file)
    
    output_handle.write(delim.join(["Num Samples", str(len(samples))]) + "\n")
    output_handle.write(delim.join(["Total Samples", str(len(samples))]) + "\n")
    output_handle.write("[Data]\n")
    
    matrix_forward = {}
    matrix_forward_GC = {}
    matrix_top = {}
    matrix_top_GC = {}
    matrix_plus = {}
    matrix_plus_GC = {}
    matrix_AB = {}
    matrix_AB_GC = {}
    
    for gtc_file in samples:
        samplename = os.path.basename(gtc_file)[:-4]
        sys.stderr.write("Processing " + gtc_file + "\n")
        gtc_file = os.path.join(args.gtc_directory, gtc_file)
        gtc = GenotypeCalls.GenotypeCalls(gtc_file)
        genotypes = [code2genotype[genotype] for genotype in gtc.get_genotypes()]
        assert len(genotypes) == len(manifest.names)
        if args.plus or args.plus_GC:
            plus_strand_genotypes = gtc.get_base_calls_plus_strand(manifest.snps, manifest.ref_strands)
        if args.forward or args.forward_GC:
            forward_strand_genotypes = gtc.get_base_calls_forward_strand(manifest.snps, manifest.source_strands)
        if args.top or args.top_GC:
            top_strand_genotypes = gtc.get_base_calls_TOP_strand(manifest.snps, manifest.ilmn_strands)
        if args.forward_GC or args.top_GC or args.plus_GC or args.AB_GC:
            genotype_scores = around(gtc.get_genotype_scores(), decimals = 4)
        #build dictionary for pandas
        if args.forward_GC:
            matrix_forward_GC = build_dict(matrix = matrix_forward_GC, samplename = samplename, \
                names = manifest.names, genotypes = forward_strand_genotypes, genotype_scores = genotype_scores)
        elif args.forward:
            matrix_forward = build_dict(matrix = matrix_forward, samplename = samplename, \
                names = manifest.names, genotypes = forward_strand_genotypes)
        elif args.top:
            matrix_top = build_dict(matrix = matrix_top, samplename = samplename, \
                names = manifest.names, genotypes = top_strand_genotypes)
        elif args.top_GC:
            matrix_top_GC = build_dict(matrix = matrix_top_GC, samplename = samplename, \
                names = manifest.names, genotypes = top_strand_genotypes, genotype_scores = genotype_scores)
        elif args.plus_GC:
            matrix_plus_GC = build_dict(matrix = matrix_plus_GC, samplename = samplename, \
                names = manifest.names, genotypes = plus_strand_genotypes, genotype_scores = genotype_scores)
        elif args.plus:
            matrix_plus = build_dict(matrix = matrix_plus, samplename = samplename, \
                names = manifest.names, genotypes = plus_strand_genotypes)
        elif args.AB:
            matrix_AB = build_dict(matrix = matrix_AB, samplename = samplename, \
                names = manifest.names, genotypes = genotypes)
        elif args.AB_GC:
            matrix_AB_GC = build_dict(matrix = matrix_AB_GC, samplename = samplename, \
                names = manifest.names, genotypes = genotypes, genotype_scores = genotype_scores)
    #create pandas dataframe from dictionaryand append to file
    if args.forward_GC:
        df = DataFrame.from_dict(matrix_forward_GC)
    elif args.forward:
        df = DataFrame.from_dict(matrix_forward)
    elif args.top:
        df = DataFrame.from_dict(matrix_top)
    elif args.top_GC:
        df = DataFrame.from_dict(matrix_top_GC)
    elif args.plus_GC:
        df = DataFrame.from_dict(matrix_plus_GC)
    elif args.plus:
        df = DataFrame.from_dict(matrix_plus)
    elif args.AB:
        df = DataFrame.from_dict(matrix_AB)
    elif args.AB_GC:
        df = DataFrame.from_dict(matrix_AB_GC)
    df = df.reindex(manifest.names)
    df.to_csv(output_handle, sep = delim)
