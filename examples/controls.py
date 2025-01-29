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
    # Write header and grab control_config from manifest
    output_handle.write(','.join(['Category','Control','BeadType','Sample_ID','Sentrix_Label']))
    controls = manifest.control_config.split('\n')
    # Trim empty control if exists
    if controls[-1] == '':
        controls.pop(len(controls)-1)
    num_sections = len(controls[0].split(',')[0].split(':'))
    for i in range(1, num_sections+1):
        output_handle.write(',')
        output_handle.write(','.join([f'Section {i} X', f'Section {i} Y']))
    output_handle.write('\n')

    # Build controls dictionary from gtcs
    samples = {}
    for filename in os.listdir(args.gtc_directory):
        if filename.lower().endswith('.gtc'):
            base, ext = os.path.splitext(filename)
            samples[base] = {'gtc': filename}
    for sample in samples:
        sys.stdout.write('Processing ' + sample + '\n')
        gtc = GenotypeCalls(os.path.join(args.gtc_directory, samples[sample]['gtc']))
        samples[sample]['controls_x'] = gtc.get_control_x_intensities()
        samples[sample]['controls_y'] = gtc.get_control_y_intensities()

    # Write out controls in genome studio "ControlDashboard.csv" format
    sys.stdout.write(f'Writing to {args.output_file}\n')
    control_offset = 0
    for control_data in controls:
        beadtypes, category, color, control = control_data.split(',')
        # assumes all same beadtype, converts to int then back to str to trim off preceeding zeroes
        beadtype = str(int(beadtypes.split(':')[0]))
        for sample in sorted(samples):
            output_handle.write(','.join([category, control, beadtype, sample, sample]))
            controls_x = samples[sample]['controls_x'][control_offset:control_offset + num_sections]
            controls_y = samples[sample]['controls_y'][control_offset:control_offset + num_sections]
            for x, y in zip(controls_x, controls_y):
                output_handle.write(',')
                output_handle.write(f'{x},{y}')
            output_handle.write('\n')
        control_offset = control_offset + num_sections

