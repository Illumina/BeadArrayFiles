##    Copyright (c) 2016, Illumina
##    All rights reserved.
##
##    Redistribution and use in source and binary forms, with or without
##    modification, are permitted provided that the following conditions are met:
##
##    1. Redistributions of source code must retain the above copyright notice, this
##       list of conditions and the following disclaimer.
##    2. Redistributions in binary form must reproduce the above copyright notice,
##       this list of conditions and the following disclaimer in the documentation
##       and/or other materials provided with the distribution.
##
##    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
##    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
##    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
##    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
##    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
##    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
##    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
##    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
##    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
##    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##    The views and conclusions contained in the software and documentation are those
##    of the authors and should not be interpreted as representing official policies,
##    either expressed or implied, of the FreeBSD Project.

__version__ = "1.0.0"
import struct
from math import cos,sin,pi,atan2


code2genotype = [
"NC",
"AA",
"AB",
"BB",
"NULL",
"A",
"B",
"AAA",
"AAB",
"ABB",
"BBB",
"AAAA",
"AAAB",
"AABB",
"ABBB",
"BBBB",
"AAAAA",
"AAAAB",
"AAABB",
"AABBB",
"ABBBB",
"BBBBB",
"AAAAAA",
"AAAAAB",
"AAAABB",
"AAABBB",
"AABBBB",
"ABBBBB",
"BBBBBB",
"AAAAAAA",
"AAAAAAB",
"AAAAABB",
"AAAABBB",
"AAABBBB",
"AABBBBB",
"ABBBBBB",
"BBBBBBB",
"AAAAAAAA",
"AAAAAAAB",
"AAAAAABB",
"AAAAABBB",
"AAAABBBB",
"AAABBBBB",
"AABBBBBB",
"ABBBBBBB",
"BBBBBBBB"
]

NC = chr(0)
AA = chr(1)
AB = chr(2)
BB = chr(3)


class GenotypeCalls:
    """Class to parse gtc files as produced by Illumina AutoConvert
    and AutoCall software.

    Attributes:
        supported_versions: Supported file versions as a list of integers
    """
    __ID_NUM_SNPS = 1
    __ID_PLOIDY = 2
    __ID_PLOIDY_TYPE = 3
    __ID_SAMPLE_NAME = 10
    __ID_SAMPLE_PLATE = 11
    __ID_SAMPLE_WELL = 12
    __ID_CLUSTER_FILE = 100
    __ID_SNP_MANIFEST = 101
    __ID_IMAGING_DATE = 200
    __ID_AUTOCALL_DATE = 201
    __ID_AUTOCALL_VERSION = 300
    __ID_NORMALIZATION_TRANSFORMS = 400    
    __ID_CONTROLS_X = 500
    __ID_CONTROLS_Y = 501
    __ID_RAW_X = 1000
    __ID_RAW_Y = 1001
    __ID_GENOTYPES = 1002
    __ID_BASE_CALLS = 1003
    __ID_GENOTYPE_SCORES = 1004
    __ID_SCANNER_DATA = 1005
    __ID_CALL_RATE = 1006
    __ID_GENDER = 1007
    __ID_LOGR_DEV = 1008
    __ID_GC10 = 1009
    __ID_GC50 = 1011
    __ID_B_ALLELE_FREQS = 1012
    __ID_LOGR_RATIOS = 1013
    __ID_PERCENTILES_X = 1014
    __ID_PERCENTILES_Y = 1015
    __ID_SLIDE_IDENTIFIER = 1016
    
    supported_version = [3,4,5];
    def __init__(self, filename, ignore_version = False):
        """Constructor

        Args:
            filename (string): GTC filename 
            ignore_version (bool): boolean to ignore automated checks on
                            file version, not recommended (default: False)
        Returns:
            GenotypeCalls object
        """
        self.filename = filename
        with open(self.filename, "rb") as gtc_handle:
            identifier = gtc_handle.read(3)
            if identifier != "gtc":
                raise Exception("GTC format error: bad format identifier")
            self.version = read_byte(gtc_handle)
            if self.version not in GenotypeCalls.supported_version and not ignore_version:
                raise Exception("Unsupported GTC File version ("+str(self.version)+")")
            number_toc_entries = read_int(gtc_handle)

            #
            # Parse the table of contents and map the toc entry
            # to the lookup
            #
            self.toc_table = {}
            for toc_idx in xrange(number_toc_entries):
                (id, offset) = struct.unpack("<hI",gtc_handle.read(6))
                self.toc_table[id] = offset
    
    def __get_generic(self, toc_entry, parse_function):
        """Internal helper function to access a data element in a
        generic fashion.

        Args:
            toc_entry (int): Identifier for entry in table of contents
            parse_function (function): Used to parse the value
                                         from a file handle
        Returns:
            A single value parsed from the file (type dependent on
            parse_function)
        """
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[toc_entry])
            return parse_function(gtc_handle)
    
    def __get_generic_array(self, toc_entry, parse_function):
        """Internal helper function to access a data array in a generic
        fashion.

        Args:
            toc_entry (int): Identifier for entry in table of contents
            parse_function (function): A function used to parse the value
                                         from a file handle
        Returns:
            An array parsed from the file (type dependent on parse_function)
        """
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[toc_entry])
            num_entries = read_int(gtc_handle)
            result = []
            for idx in xrange(num_entries):
                result.append( parse_function( gtc_handle ) )
            return result
    

    def get_num_snps(self):
        """Returns:
            The number of SNPs in the file as an integer 
        """
        return self.toc_table[GenotypeCalls.__ID_NUM_SNPS]

    def get_ploidy(self):
        """Returns:
            The ploidy of the sample 
        """
        return self.toc_table[GenotypeCalls.__ID_PLOIDY]

    def get_ploidy_type(self):
        """Returns:
            The ploidy of the sample 
        """
        return self.toc_table[GenotypeCalls.__ID_PLOIDY_TYPE]
        
    
    def get_sample_name(self):
        """Returns:
            The name of the sample as a string
        """
        return self.__get_generic(GenotypeCalls.__ID_SAMPLE_NAME, read_string)

    def get_slide_identifier(self):
        """Returns:
            The name of the sample as a string
        """
        return self.__get_generic(GenotypeCalls.__ID_SLIDE_IDENTIFIER, read_string)


    def get_sample_plate(self):
        """Returns:
            The name of the sample plate as a string
        """
        return self.__get_generic(GenotypeCalls.__ID_SAMPLE_PLATE, read_string)

    def get_sample_well(self):
        """Returns:
            The name of the sample well as a string
        """        
        return self.__get_generic(GenotypeCalls.__ID_SAMPLE_WELL, read_string)
    
    def get_cluster_file(self):
        """Returns:
            The name of the cluster file used for genotyping as a string
        """    
        return self.__get_generic(GenotypeCalls.__ID_CLUSTER_FILE, read_string)

    def get_snp_manifest(self):
        """Returns:
            The name of the manifest used for genotyping as a string
        """   
        return self.__get_generic(GenotypeCalls.__ID_SNP_MANIFEST, read_string)

    def get_imaging_date(self):
        """Returns:
            The imaging date of scanning as a string
            For example
                Monday, December 01, 2014 4:51:47 PM
        """         
        return self.__get_generic(GenotypeCalls.__ID_IMAGING_DATE, read_string)

    def get_autocall_date(self):
        """Returns:
            The imaging date of scanning as a string
            For example
                2/17/2015 1:47 PM
        """         
        return self.__get_generic(GenotypeCalls.__ID_AUTOCALL_DATE, read_string)

    def get_autocall_version(self):
        """Returns:
            The version of AutoCall used for genotyping as a string
            For example
                1.6.2.2
        """
        return self.__get_generic(GenotypeCalls.__ID_AUTOCALL_VERSION, read_string)

            
    def get_genotypes(self):
        """ Returns:
            A byte list (string) of genotypes. See code2genotype for mapping
        """
        return self.__get_generic_array(GenotypeCalls.__ID_GENOTYPES, read_byte)

    def get_base_calls(self):
        """Returns:
            The genotype basecalls as a list of strings. 
            The characters are A, C, G, T, or - for a no-call/null.
            The calls are relative to the top strand.
        """
        try:
            ploidy_type = self.get_ploidy_type()
        except:
            ploidy_type = 1
            
	if ploidy_type != 1:
		genotypes = self.get_genotypes()
		
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_BASE_CALLS])
            num_entries = read_int(gtc_handle)
            result = []
            for idx in xrange(num_entries):
            	if ploidy_type == 1:
                	result.append( gtc_handle.read(2) )
                else:
                	byte_string = gtc_handle.read(2)
                	ab_genotype = code2genotype[genotypes[idx]]
                	if ab_genotype == "NC" or ab_genotype == "NULL":
                		result.append("-")
                	else:
                		top_genotype = "".join([ byte_string[0] if allele == "A" else byte_string[1] for allele in ab_genotype])
                		result.append( top_genotype )
            return result

    def get_genotype_scores(self):
        """Returns:
            The genotype scores as a list of floats
        """
        return self.__get_generic_array(GenotypeCalls.__ID_GENOTYPE_SCORES, read_float)

    def get_scanner_data(self):
        """Returns:
            Information about scanner as ScannerData object
        """
        return self.__get_generic(GenotypeCalls.__ID_SCANNER_DATA, read_scanner_data)

    def get_control_x_intensities(self):
        """Returns:
            The x intensities of control bead types as a list of integers
        """
        return self.__get_generic_array(GenotypeCalls.__ID_CONTROLS_X, read_ushort)
    
    def get_control_y_intensities(self):
        """Returns:
            The y intensities of control bead types as a list of integers
        """
        return self.__get_generic_array(GenotypeCalls.__ID_CONTROLS_Y, read_ushort)

    def get_raw_x_intensities(self):
        """Returns:
            The raw x intensities of assay bead types as a list of integers
        """
        return self.__get_generic_array(GenotypeCalls.__ID_RAW_X, read_ushort)
    
    def get_raw_y_intensities(self):
        """Returns:
            The raw y intensities of assay bead types as a list of integers
        """
        return self.__get_generic_array(GenotypeCalls.__ID_RAW_Y, read_ushort)

    def get_call_rate(self):
        """Returns:
            The call rate as a float
        """
        return self.__get_generic(GenotypeCalls.__ID_CALL_RATE, read_float)

    def get_gender(self):
        """Returns:
            The gender as a char
            M - Male, F - Female, U-Unknown
        """
        return self.__get_generic(GenotypeCalls.__ID_GENDER, read_char)
        
    def get_logr_dev(self):
        """Returns:
            The logR deviation as a float
        """        
        return self.__get_generic(GenotypeCalls.__ID_LOGR_DEV, read_float)

    def get_gc10(self):
        """Returns:
            The GC10 (GenCall score - 10th percentile) as a float
        """
        return self.__get_generic(GenotypeCalls.__ID_GC10, read_float)

    def get_gc50(self):
        """Returns:
            The GC50 (GenCall score - 50th percentile) as a float
        """
        return self.__get_generic(GenotypeCalls.__ID_GC50, read_float)

    def get_num_calls(self):
        """Returns:
            The number of calls as an integer
        """
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_GC50] + 4)
            return read_int(gtc_handle)

    def get_num_no_calls(self):
        """Returns:
            The number of no calls as an integer
        """
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_GC50] + 8)
            return read_int(gtc_handle)

    def get_num_intensity_only(self):
        """Returns:
            The number of intensity only SNPs
        """
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_GC50] + 12)
            return read_int(gtc_handle)
    
    def get_ballele_freqs(self):
        """Returns:
            The B allele frequencies as a list of floats
        """
        if self.version < 4:
            raise Exception("B allele frequencies unavailable in GTC File version ("+str(self.version)+")")
        return self.__get_generic_array(GenotypeCalls.__ID_B_ALLELE_FREQS, read_float)

    def get_logr_ratios(self):
        """Returns:
            The logR ratios as a list of floats
        """        
        if self.version < 4:
            raise Exception("LogR ratios unavailable in GTC File version ("+str(self.version)+")")
        return self.__get_generic_array(GenotypeCalls.__ID_LOGR_RATIOS, read_float)

    def get_percentiles_x(self):
        """Returns:
            An array of length three representing 5th, 50th and 95th percentiles for
            x intensity
        """
        if self.version < 5:
            raise Exception("Percentile intensities unavailable in GTC File version ("+str(self.version)+")")
        result = []
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_PERCENTILES_X])
            result = []
            for idx in xrange(3):
                result.append(read_ushort(gtc_handle))
            return result

    def get_percentiles_y(self):
        """Returns:
            An array of length three representing 5th, 50th and 95th percentiles for
            y intensity
        """
        if self.version < 5:
            raise Exception("Percentile intensities unavailable in GTC File version ("+str(self.version)+")")
        result = []
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_PERCENTILES_Y])
            result = []
            for idx in xrange(3):
                result.append(read_ushort(gtc_handle))
            return result
       
    def get_normalized_intensities(self, normalization_lookups):
        """Calculate and return the normalized intensities
        Args:
            normalization_lookups (list of ints): Map from each SNP to a normalization transform.
                                   This list can be obtained from the BeadPoolManifest object.
                                   
        Return:
            The normalized intensities for the sample as a list of (x,y) float tuples
        """
        normalization_transforms = self.get_normalization_transforms()
        return [normalization_transforms[lookup].normalize_intensities(x_raw, y_raw) for (x_raw, y_raw, lookup) in zip(self.get_raw_x_intensities(), self.get_raw_y_intensities(), normalization_lookups)]  
    
    def get_normalization_transforms(self):
        """Returns:
            The normalization transforms used during genotyping (as a lit of NormalizationTransforms)
        """
        return self.__get_generic_array(GenotypeCalls.__ID_NORMALIZATION_TRANSFORMS, NormalizationTransform.read_normalization_transform)
    	



class NormalizationTransform:
    def __init__(self, buffer):
        """Constructor
        Args:
            buffer (string): byte string containing binary data for transform 
        Returns:
            NormalizationTransform object
        """
        (self.version, self.offset_x, self.offset_y, self.scale_x, self.scale_y, self.shear, self.theta,) = struct.unpack("<iffffff", buffer[:28])
    
    @staticmethod
    def read_normalization_transform(handle):
        """Static helper function to read normalization transform from file handle
        Args:
            handle (file handle): File handle with position at start of normalization transform entry
        Returns:
            NormalizationTransform object
        """        
        return NormalizationTransform(handle.read(52))

    @staticmethod
    def rect_to_polar((x,y)):
        """Converts x,y intensities to (pseudo) polar co-ordinates (R, theta)
        Args:
            x,y (int,int): Raw x,y intensities for probe
        Returns:
            (R,theta) polar values as tuple of floats
        """
        if x == 0 and y == 0:
            return (float("nan"), float("nan"))
        return (x + y, atan2(y,x) * 2.0 / pi)

    def normalize_intensities(self, x, y, threshold = True):
        """Apply this normalization transform to raw intensities
        Args:
            x (int): Raw x intensities
            y (int): Raw y intensities
        Returns:
            (xn, yn) normalized intensities as tuple of floats
        """
        tempx = x - self.offset_x
        tempy = y - self.offset_y

        tempx2 =  cos(self.theta) * tempx + sin(self.theta) * tempy
        tempy2 = -sin(self.theta) * tempx + cos(self.theta) * tempy

        tempx3 = tempx2 - self.shear * tempy2
        tempy3 = tempy2

        xn = tempx3 / self.scale_x
        yn = tempy3 / self.scale_y

        if threshold:
            xn = max(0, xn)
            yn = max(0, yn)
            
        return (xn, yn)


class ScannerData:
    def __init__(self, name, pmt_green, pmt_red, version, user):
        """Constructor
        Args:
            name (string): scanner identifier
            pmt_green (int): gain setting (green channel)
            pmt_red (int): gain setting (red channel)
            version (string): version of scanner software
            user (string): user of the scanner software
        Returns:
            ScannerData object
        """
        self.name = name
        self.pmt_green = pmt_green
        self.pmt_red = pmt_red
        self.version = version
        self.user = user


class BeadPoolManifest:
    """Class for parsing binary (BPM) manifest file.
    Attributes:
        names (list of strings): Names of loci from manifest
        addresses (list of ints): AddressA IDs of loci from manifest
        normalization_lookups (list of ints): Normalization lookups from manifest. This indexes into
                                         list of normalization transforms read from GTC file
        num_loci (int): Number of loci in manifest
        manifest_name (string): Name of manifest
        control_config (string): Control description from manifest
    """
    def __init__(self, filename):
        """Constructor
        Args:
            filename (string): Locations of BPM (bead pool manifest) file 
        Returns:
            BeadPoolManifest
        """        
        self.names = []
        self.addresses = []
        self.normalization_ids = []
        self.assay_types = []
        self.normalization_lookups = []
        self.num_loci = 0
        self.manifest_name = ""
        self.control_config = ""
        self.__parse_file(filename)

    def __parse_file(self, manifest_file):
        """Helper function to initialize this object from a file.
        Args:
            manifest_file (string): Location of BPM (bead pool manifest) file
        Returns:
            None
        """
        with open(manifest_file, "rb") as manifest_handle:
            header = manifest_handle.read(3)
            if len(header) != 3 or header != "BPM":
                raise Exception("Invalid BPM format")
            version = read_byte(manifest_handle)
            if version != 1:
                raise Exception("Unknown BPM version (" + str(ord(version)) + ")")

            version = read_int(manifest_handle)
            version_flag = 0x1000
            if version & version_flag == version_flag:
                version = version ^ version_flag
            if version > 5 or version < 3:
                raise Exception("Unsupported BPM version (" + str(version) + ")")
            self.manifest_name = read_string(manifest_handle)
            
            if version > 1:
                self.control_config = read_string(manifest_handle)

            self.num_loci = read_int(manifest_handle)
            manifest_handle.seek(4 * self.num_loci, 1)
            name_lookup = {}
            for idx in xrange(self.num_loci):
                self.names.append(read_string(manifest_handle))
                name_lookup[self.names[-1]] = idx

            for idx in xrange(self.num_loci):
                normalization_id = read_byte(manifest_handle)
                if normalization_id >= 100:
                    raise Exception("Manifest format error: read invalid normalization ID")
                self.normalization_ids.append(normalization_id)
            
            self.assay_types = [0] * self.num_loci
            self.addresses = [0] * self.num_loci
            for idx in xrange(self.num_loci):
                locus_entry = LocusEntry(manifest_handle)
                self.assay_types[ name_lookup[locus_entry.name] ] = locus_entry.assay_type
                self.addresses[ name_lookup[locus_entry.name] ] = locus_entry.address_a
            
            if len(self.normalization_ids) != len(self.assay_types):
                raise Exception("Manifest format error: read invalid number of assay entries")

            all_norm_ids = set()
            for locus_idx in xrange(self.num_loci):
                self.normalization_ids[locus_idx] += 100 * self.assay_types[locus_idx]
                all_norm_ids.add(self.normalization_ids[locus_idx])
            sorted_norm_ids = sorted(all_norm_ids)
            lookup_dictionary = {}
            for idx in xrange(len(sorted_norm_ids)):
                lookup_dictionary[sorted_norm_ids[idx]] = idx
            self.normalization_lookups =  [lookup_dictionary[normalization_id] for normalization_id in self.normalization_ids]

class LocusEntry:
    """Helper class representing a locus entry within a bead pool manifest. Current only support version
    6,7, and 8.
    Attributes:
        ilmn_id (string) : IlmnID (probe identifier) of locus
        name (string): Name (variant identifier) of locus
        assay_type (int) : Identifies type of assay (0 - Infinium II , 1 - Infinium I (A/T), 2 - Infinium I (G/C)
        address_a (int) : AddressA ID of locus
        address_b (int) : AddressB ID of locus (0 if none)
    """
    def __init__(self, handle):
        """Constructor
        Args:
            handle (file handle):  File handle at start of locus entry record
        Returns:
            LocusEntry
        """ 
        self.ilmn_id = ""
        self.name = ""  
        self.assay_type = -1
        self.address_a = -1
        self.address_b = -1
        self.__parse_file(handle)

    def __parse_file(self, handle):
        """Helper function to initialize this object from a file handle
        Args:
            handle (file handle): File handle at start of locus entry record
        Returns:
            None
        """     
        version = read_int(handle)
        if version == 6:
            self.__parse_locus_version_6(handle)
        elif version == 7:
            self.__parse_locus_version_7(handle)
        elif version == 8:
            self.__parse_locus_version_8(handle)
        else:
            raise Exception("Manifest format error: unknown version for locus entry (" + str(version) + ")")

    def __parse_locus_version_6(self,handle):
        """Helper function to parse version 6 locus entry
        Args:
            handle (file handle): File handle at start of locus entry record
        Returns:
            None
        """  
        self.ilmn_id = read_string(handle)
        self.name = read_string(handle)
        for idx in xrange(3):
            read_string(handle)
        handle.read(4)
        for idx in xrange(9):
            read_string(handle)
        self.address_a = read_int(handle)
        self.address_b = read_int(handle)
        for idx in xrange(7):
            read_string(handle)
        handle.read(3)
        self.assay_type = read_byte(handle)
        if self.assay_type not in [0,1,2]:
            raise Exception("Format error in reading assay type from locus entry")
        if self.address_b == 0:
            if self.assay_type != 0:
                raise Exception("Manifest format error: Assay type is inconsistent with address B")
        else:
            if self.assay_type == 0:
                raise Exception("Manifest format error: Assay type is inconsistent with address B")

    def __parse_locus_version_7(self, handle):
        """Helper function to parse version 7 locus entry
        Args:
            handle (file handle): File handle at start of locus entry record
        Returns:
            None
        """  
        self.__parse_locus_version_6(handle)
        handle.read(4 * 4)

    def __parse_locus_version_8(self, handle):
        """Helper function to parse version 8 locus entry
        Args:
            handle (file handle): File handle at start of locus entry record
        Returns:
            None
        """         
        self.__parse_locus_version_7(handle)
        read_string(handle)
    
def read_char(handle):
    """Helper function to parse character from file handle
    Args:
        handle (file handle): File handle
    Returns:
        char value read from handle
    """ 
    return handle.read(1)

def read_ushort(handle):
    """Helper function to parse ushort from file handle
    Args:
        handle (file handle): File handle
    Returns:
        ushort value read from handle
    """     
    return struct.unpack("<H", handle.read(2))[0]

def read_int(handle):
    """Helper function to parse int from file handle
    Args:
        handle (file handle): File handle
    Returns:
        int value read from handle
    """     
    return struct.unpack("<i", handle.read(4))[0]

def read_float(handle):
    """Helper function to parse float from file handle
    Args:
        handle (file handle): File handle
    Returns:
        float value read from handle
    """ 
    return struct.unpack("<f",handle.read(4))[0]

def read_byte(handle):
    """Helper function to parse byte from file handle
    Args:
        handle (file handle): File handle
    Returns:
        byte value read from handle
    """ 
    return struct.unpack("<B",handle.read(1))[0]

def read_string(handle):
    """Helper function to parse string from file handle. See https://msdn.microsoft.com/en-us/library/yzxa6408(v=vs.100).aspx
    for additional details on string format.
    Args:
        handle (file handle): File handle
    Returns:
        string value read from handle
    """    
    total_length = 0
    partial_length = read_byte(handle)
    num_bytes = 0
    while partial_length & 0x80 > 0:
        total_length += (partial_length & 0x7F) << ( 7 * num_bytes)
        partial_length = ord(struct.unpack("c", handle.read(1))[0])
        num_bytes += 1
    total_length += partial_length << ( 7 * num_bytes)
    return handle.read(total_length)

def read_scanner_data(handle):
    """Helper function to parse ScannerData object from file handle. 
    Args:
        handle (file handle): File handle
    Returns:
        ScannerData value read from handle
    """    
    name = read_string(handle)
    pmt_green = read_int(handle)
    pmt_red = read_int(handle)
    scanner_version = read_string(handle)
    imaging_user = read_string(handle)
    return ScannerData(name, pmt_green, pmt_red, scanner_version, imaging_user)

