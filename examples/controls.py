import os
import sys
import argparse
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest

parser = argparse.ArgumentParser('Generate a control dashboard report')
parser.add_argument('manifest', help='BPM manifest file')
parser.add_argument('gtc_directory', help='Directory containing GTC files')
parser.add_argument('output_file', help='Location to write report')

args = parser.parse_args()

try:
    manifest = BeadPoolManifest(args.manifest)
except:
    sys.stderr.write('Failed to read data from manifest\n')
    sys.exit(-1)

with open(args.output_file, 'w') as output_handle:
    # Write header
    output_handle.write(','.join(['Category','Control','BeadType','Sample_ID','Sentrix_Label']))
    controls = manifest.control_config.split('\n')
    num_sections = len(controls[0].split(',')[0].split(':'))
    for i in range(1, num_sections+1):
        output_handle.write(',')
        output_handle.write(','.join([f'Section {i} X', f'Section {i} Y']))
    output_handle.write('\n')

    samples = []
    for filename in os.listdir(args.gtc_directory):
        if filename.lower().endswith('.gtc'):
            samples.append(filename)

    for gtc_file in samples:
        sys.stderr.write('Processing ' + gtc_file + '\n')
        gtc_file = os.path.join(args.gtc_directory, gtc_file)
        gtc = GenotypeCalls(gtc_file)
        print(gtc_file)
        controls_x = gtc.get_control_x_intensities()
        print(controls_x)
        controls_y = gtc.get_control_y_intensities()
        print(controls_y)
