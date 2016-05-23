# IlluminaBeadArrayFiles
Library to parse file formats related to Illumina bead arrays

## Build instructions
The IlluminaBeadArrayFiles repository supports building a package with the python distutils module (https://docs.python.org/2/distutils/setupscript.html). To build a source distribution, run the included setup.py script supplying the "sdist" command

>python setup.py sdist

## Installation instructions
After unpacking the installation package, run the setup.py script supplying the "install" command

>python setup.py install

If the user prefers not to use the python distutils framework, it is also possible to copy the IlluminaBeadArrayFiles.py source file into a location referenced by the PYTHONPATH environment variable.

## GTC File Format
The specification for the GTC file format is provided in docs/GTC_File_Format_v5.pdf

## Description of classes and objects in exposed in package
For further details on specific class methods, please consult the built-in docstring documentation

### code2genotype
Dictionary mapping from genotype byte code (see GTC file format specification) to a string representing the genotype (e.g., "AA")

### NC, AA, AB, and BB
Constants representing byte values for associated genotypes

### GenotypeCalls
Class to parse GTC files as produced by Illumina AutoConvert and AutoCall software.

### NormalizationTransform
Class to encapsulate a normalization transform for conversion of raw intensity data to normalized intensity data

### ScannerData
Class to encapsulate the meta-data collected in the GTC file for a scanner instrument

### BeadPoolManifest
Class for parsing a binary (BPM) manifest file. This class can be used to recover the in-order list of loci names in a given manifest, which is necessary to associate the genotype information present in the GTC file to a specific locus name. 





