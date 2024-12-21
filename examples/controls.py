import sys
import argparse
import os
from datetime import datetime
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype

delim = "\t"

parser = argparse.ArgumentParser("Generate a control dashboard report")
parser.add_argument("manifest", help="BPM manifest file")
parser.add_argument("gtc_directory", help="Directory containing GTC files")
parser.add_argument("output_file", help="Location to write report")

args = parser.parse_args()

try:
    manifest = BeadPoolManifest(args.manifest)
except:
    sys.stderr.write("Failed to read data from manifest\n")
    sys.exit(-1)

with open(args.output_file, "w") as output_handle:
    controls = manifest.control_config.split('\n')
    samples = []
    for gtc_file in os.listdir(args.gtc_directory):
        if gtc_file.lower().endswith(".gtc"):
            samples.append(gtc_file)

    for gtc_file in samples:
        sys.stderr.write("Processing " + gtc_file + "\n")
        gtc_file = os.path.join(args.gtc_directory, gtc_file)
        gtc = GenotypeCalls(gtc_file)
        print(gtc_file)
        controls_x = gtc.get_control_x_intensities()
        print(controls_x)
        controls_y = gtc.get_control_y_intensities()
        print(controls_y)


