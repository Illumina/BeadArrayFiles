# IlluminaBeadArrayFiles
Library to parse binary file formats related to Illumina bead arrays. The IlluminaBeadArrayFiles library provides a parser to extract information from these binary files.

## Generating GTC files
If you have intensity data files (IDATs) for which GTC files are not currently available, it is possible to generate these files using Illumina software. We recommend using [DRAGEN Array](https://support.illumina.com/array/array_software/dragen-array-secondary-analysis.html).

## Build instructions
The IlluminaBeadArrayFiles repository supports building a package with the python distutils module (https://docs.python.org/2/distutils/setupscript.html). To build a source distribution, run the included setup.py script supplying the "sdist" command

>python setup.py sdist

## Installation instructions
After unpacking the installation package, run the setup.py script supplying the "install" command

>python setup.py install

## Dependencies
The library depends on the availability of the numpy package in the python installation (http://www.numpy.org/)

## Example usage

```python
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype

gtc_file = "path_to_genotypes.gtc"
manifest_file = "path_to_manifest.bpm"
names = BeadPoolManifest( manifest_file ).names
genotypes = GenotypeCalls( gtc_file ).get_genotypes()

c = 0  # a counter to only show the first 10 genotypes
for (locus, genotype) in zip( names, genotypes ):
    print( locus + "," + code2genotype[genotype] )
    if c >= 10:
        break
    c += 1
```

Also, see examples/* for additional examples of usage.
These scripts are based on common Genome Studio (https://support.illumina.com/array/array_software/genomestudio.html) reports.

**NOTE:**
For the DNA summary report, it does not exclude zeroed loci from the cluster file as it does in Genome Studio, so overall call rates may look lower.
There is an open issue to address this: https://github.com/Illumina/BeadArrayFiles/issues/22

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

### RefStrand, SourceStrand
Represents different strand conventions in manifest

### BeadPoolManifest
Read information from a binary BPM file

### LocusAggregate
Aggregate information across many samples for a given loci

### ClusterFile
Read information from a binary EGT file

## Compatibility
This library is compatible with Python 3 and was tested with Python 3.8.5

## License

>Copyright (c) 2020, Illumina
> All rights reserved.
>
> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are met:
>
>1. Redistributions of source code must retain the above copyright notice, this
>list of conditions and the following disclaimer.
>2. Redistributions in binary form must reproduce the above copyright notice,
>this list of conditions and the following disclaimer in the documentation
>and/or other materials provided with the distribution.
>
>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
>ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
>WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
>DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
>ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
>(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
>LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
>ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
>(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
>SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
>
>The views and conclusions contained in the software and documentation are those
>of the authors and should not be interpreted as representing official policies,
>either expressed or implied, of the FreeBSD Project.

