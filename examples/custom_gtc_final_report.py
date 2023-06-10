import sys
import argparse
import os
from datetime import datetime
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, ClusterFile, __version__

delim = "\t"

parser = argparse.ArgumentParser("Generate a final report from a directory of GTC files")
parser.add_argument("manifest", help="BPM manifest file")
parser.add_argument("cluster", help="EGT cluster file")
parser.add_argument("gtc_directory", help="Directory containing GTC files")
parser.add_argument("samplesheet", help="Genome Studio formatted samplesheet to pull Sample_ID and Sample_Name")
parser.add_argument("output_file", help="Location to write report")
parser.add_argument("--exclude-samples", dest="exclude_samples", default=[], nargs="*", help="List of sample names to exclude (optional)")
parser.add_argument("--gc-cutoff", dest="gc_cutoff", type=float, default=0.15, help="GenCall score cutoff to set GT/GC scores to 0")
parser.add_argument("--exclude-zeroed-snps", dest="exclude_zeroed_snps", action="store_true", default=False, help="Exclude SNPs with a cluster total score of 0.")
args = parser.parse_args()

if os.path.isfile(args.output_file):
    sys.stderr.write("Output file already exists, please delete and re-run\n")
    sys.exit(-1)

with open(args.cluster, "rb") as cluster_handle:
    egt = ClusterFile.read_cluster_file(cluster_handle)

try:
    manifest = BeadPoolManifest(args.manifest)
except:
    sys.stderr.write("Failed to read data from manifest\n")
    sys.exit(-1)

try:
    sample_map = {}
    with open(args.samplesheet, "r") as f:
        found_data = False
        for line in f.readlines():
            if line.startswith("Sample_ID"):
                found_data = True
                continue
            if not found_data:
                continue
            (sample_id, sample_name, sentrix_barcode, sentrix_position, _) = line.split(',')
            sample_map[sentrix_barcode+'_'+sentrix_position] = (sample_id, sample_name)
except:
    sys.stderr.write("Failed to read data from the samplesheet, is it formatted correctly with a Sample_ID column in the header?\n")
    sys.exit(-1)

with open(args.output_file, "w") as output_handle:
    output_handle.write("[Header]\n")
    output_handle.write(delim.join(["BeadArrayFiles Version", __version__])+ "\n")
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
    output_handle.write(delim.join(["Sample ID", "Sample Name", "SNP Name", "Allele1 - Top", "Allele2 - Top", "GC Score", "GT Score"]) + "\n")
    for gtc_file in samples:
        sys.stderr.write("Processing " + gtc_file + "\n")
        sample = os.path.basename(gtc_file)[:-4]
        try:
            (sample_id, sample_name) = sample_map[sample]
        except:
            sys.stderr.write("Failed to find: "+sample+" in the samplesheet.\n")
            sys.exit(-1)
        gtc_file = os.path.join(args.gtc_directory, gtc_file)
        gtc = GenotypeCalls(gtc_file)
        if sample_name in args.exclude_samples:
            continue
        top_strand_genotypes = gtc.get_base_calls()
        gc_scores = gtc.get_genotype_scores()

        assert len(top_strand_genotypes) == len(manifest.names)
        for (name, genotype, score) in zip(manifest.names, top_strand_genotypes, gc_scores):
            if args.exclude_zeroed_snps and egt.get_record(name).cluster_score.total_score == 0.0:
                continue
            top_allele1, top_allele2 = genotype.decode('ascii')

            gc_score = '{:.4f}'.format(0)
            gt_score = '{:.4f}'.format(0)
            if score >= args.gc_cutoff:
                gc_score = '{:.4f}'.format(round(score, 4))
                gt_score = '{:.4f}'.format(round(egt.get_record(name).cluster_score.total_score, 4))
            output_handle.write(delim.join([sample_id, sample_name, name, top_allele1, top_allele2, gc_score, gt_score]) + "\n")
