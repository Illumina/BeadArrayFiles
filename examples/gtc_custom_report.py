import sys
import argparse
import os
from datetime import datetime
#import GenotypeCalls
from IlluminaBeadArrayFiles import GenotypeCalls

delim = "\t"

parser = argparse.ArgumentParser("Generate a custom report from a directory of GTC files")
parser.add_argument("gtc_directory", help="Directory containing GTC files")
parser.add_argument("output_file", help="Location to write report")

args = parser.parse_args()

if os.path.isfile(args.output_file):
    sys.stderr.write("Output file already exists, please delete and re-run\n")
    sys.exit(-1)

with open(args.output_file, "w") as output_handle:
    output_handle.write("[Header]\n")
    output_handle.write(delim.join(["Processing Date", datetime.now().strftime("%m/%d/%Y %I:%M %p")])+ "\n")
    samples = []
    for gtc_file in os.listdir(args.gtc_directory):
        if gtc_file.lower().endswith(".gtc"):
            samples.append(gtc_file)
    index = 1
    output_handle.write(delim.join(["Num Samples", str(len(samples))]) + "\n")
    output_handle.write(delim.join(["Total Samples", str(len(samples))]) + "\n")
    output_handle.write("[Data]\n")
    output_handle.write(delim.join(["Index","Sample ID", "Call Rate", "p05 Grn", "p50 Grn", "p90 Grn", "p05 Red", "p50 Red", \
            "p90 Red","p10 GC", "p50 GC"]) + "\n")
    for gtc_file in samples:
        sys.stderr.write("Processing " + gtc_file + "\n")
        gtc_file = os.path.join(args.gtc_directory, gtc_file)
        gtc = GenotypeCalls.GenotypeCalls(gtc_file)
        call_rate = gtc.get_call_rate()
        percentiles_grn = gtc.get_percentiles_x()
        p05_grn = percentiles_grn[0]
        p50_grn = percentiles_grn[1]
        p90_grn = percentiles_grn[2]
        percentiles_red = gtc.get_percentiles_y()
        p05_red = percentiles_red[0]
        p50_red = percentiles_red[1]
        p90_red = percentiles_red[2]
        p10_GC = gtc.get_gc10()
        p50_GC = gtc.get_gc50()
        output_handle.write(delim.join([str(index), os.path.basename(gtc_file)[:-4], str(call_rate), str(p05_grn), \
            str(p50_grn), str(p90_grn), str(p05_red), str(p50_red), str(p90_red), str(p10_GC), str(p50_GC)]) + "\n")
        index +=1
