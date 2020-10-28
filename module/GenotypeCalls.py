import struct
from numpy import cos, sin, pi, arctan2, float32, uint16, int32, seterr, frombuffer, dtype
from .BeadArrayUtility import read_int, read_string, read_byte, read_float, read_char, read_ushort, complement
from .BeadPoolManifest import RefStrand, SourceStrand

seterr(divide='ignore', invalid='ignore')
nan = float32('nan')

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

class GenotypeCalls(object):
    """
    Class to parse gtc files as produced by Illumina AutoConvert
    and AutoCall software.

    Attributes:
        supported_versions (list(int)): Supported file versions as a list of integers
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

    supported_version = [3, 4, 5]

    def __init__(self, filename, ignore_version=False, check_write_complete=True):
        """
        Constructor

        Args:
            filename (string): GTC filename
            ignore_version (bool): boolean to ignore automated checks on
                            file version, not recommended (default: False)

        Returns:
            GenotypeCalls
        """
        self.filename = filename
        with open(self.filename, "rb") as gtc_handle:
            identifier = gtc_handle.read(3)
            identifier = identifier.decode("utf-8")
            if identifier != "gtc":
                raise Exception("GTC format error: bad format identifier")
            self.version = read_byte(gtc_handle)
            if self.version not in GenotypeCalls.supported_version and not ignore_version:
                raise Exception(
                    "Unsupported GTC File version (" + str(self.version) + ")")
            number_toc_entries = read_int(gtc_handle)

            #
            # Parse the table of contents and map the toc entry
            # to the lookup
            #
            self.toc_table = {}
            for toc_idx in range(number_toc_entries):
                (id, offset) = struct.unpack("<hI", gtc_handle.read(6))
                self.toc_table[id] = offset
        if check_write_complete and not self.is_write_complete():
            raise Exception("GTC file is incomplete")

    def __get_generic(self, toc_entry, parse_function):
        """
        Internal helper function to access a data element in a
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

    def __get_generic_array(self, toc_entry, parse_function, item_size, offset, count):
        """
        Internal helper function to access a data array in a generic
        fashion.

        Args:
            toc_entry (int): Identifier for entry in table of contents
            parse_function (function): A function used to parse the value
                                         from a file handle
            item_size (int): Size (in bytes) of individual entry
            offset (int): Offset (in elements counts) to start reading
            count (int): Number of entries to read (None is read all remaining entries)

        Returns:
            list(type): An array parsed from the file (type dependent on parse_function)
        """
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[toc_entry])
            num_entries = read_int(gtc_handle) - offset
            if count is not None:
                num_entries = min(num_entries, count)
            if offset > 0:
                gtc_handle.seek(
                    self.toc_table[toc_entry] + 4 + offset * item_size)
            result = []
            for idx in range(num_entries):
                result.append(parse_function(gtc_handle))
            return result

    def __get_generic_array_numpy(self, toc_entry, numpy_type, offset=0, count=None):
        """
        Internal helper function to access a data array in a generic
        fashion.

        Args:
            toc_entry (int): Identifier for entry in table of contents
            numpy_type (numpy.dtype): Data type to read into array
            offset (int): Offset (in element counts) to start reading
            count (int): Number of entries to read (None will read remaining entries)

        Returns:
            list(type): An array parsed from the file (type dependent on parse_function)
        """
        numpy_type = dtype(numpy_type)
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[toc_entry])
            num_entries = read_int(gtc_handle) - offset
            if count is not None:
                num_entries = min(num_entries, count)
            if offset > 0:
                gtc_handle.seek(
                    self.toc_table[toc_entry] + 4 + offset * numpy_type.itemsize)
            return frombuffer(gtc_handle.read(num_entries * numpy_type.itemsize), dtype=numpy_type)

    def get_num_snps(self):
        """
        Returns:
            int: The number of SNPs in the file
        """
        return self.toc_table[GenotypeCalls.__ID_NUM_SNPS]

    def get_ploidy(self):
        """
        Returns:
            int: The ploidy of the sample
        """
        return self.toc_table[GenotypeCalls.__ID_PLOIDY]

    def get_ploidy_type(self):
        """
        Returns:
            int: The ploidy type of the sample
        """
        return self.toc_table[GenotypeCalls.__ID_PLOIDY_TYPE]

    def get_sample_name(self):
        """
        Returns:
            string: The name of the sample
        """
        return self.__get_generic(GenotypeCalls.__ID_SAMPLE_NAME, read_string)

    def get_slide_identifier(self):
        """
        Returns:
            string: The sentrix identifier for the slide
        """
        return self.__get_generic(GenotypeCalls.__ID_SLIDE_IDENTIFIER, read_string)

    def get_sample_plate(self):
        """
        Returns:
            string: The name of the sample plate
        """
        return self.__get_generic(GenotypeCalls.__ID_SAMPLE_PLATE, read_string)

    def get_sample_well(self):
        """
        Returns:
            string: The name of the sample well
        """
        return self.__get_generic(GenotypeCalls.__ID_SAMPLE_WELL, read_string)

    def get_cluster_file(self):
        """
        Returns:
            string: The name of the cluster file used for genotyping
        """
        return self.__get_generic(GenotypeCalls.__ID_CLUSTER_FILE, read_string)

    def get_snp_manifest(self):
        """
        Returns:
            string: The name of the manifest used for genotyping
        """
        return self.__get_generic(GenotypeCalls.__ID_SNP_MANIFEST, read_string)

    def get_imaging_date(self):
        """
        Returns:
            string: The imaging date of scanning
            For example
                Monday, December 01, 2014 4:51:47 PM
        """
        return self.__get_generic(GenotypeCalls.__ID_IMAGING_DATE, read_string)

    def get_autocall_date(self):
        """
        Returns:
            string: The imaging date of autocall execution
            For example
                2/17/2015 1:47 PM
        """
        return self.__get_generic(GenotypeCalls.__ID_AUTOCALL_DATE, read_string)

    def get_autocall_version(self):
        """
        Returns:
            string: the version of AutoCall used for genotyping
            For example
                1.6.2.2
        """
        return self.__get_generic(GenotypeCalls.__ID_AUTOCALL_VERSION, read_string)

    def get_genotypes(self, offset=0, count=None):
        """
        Returns:
            string: A byte list (string) of genotypes. See code2genotype for mapping
        """
        return self.__get_generic_array(GenotypeCalls.__ID_GENOTYPES, read_byte, 1, offset, count)

    def get_base_calls_generic(self, snps, strand_annotations, report_strand, unknown_annotation):
        """
        Get base calls on arbitrary strand

        Args:
            snps (list(string)) : A list of string representing the snp on the design strand for the loci (e.g. [A/C])
            strand_annotations (list(int)) : A list of strand annotations for the loci
            report_strand (int) : The strand to use for reporting (must match encoding of strand_annotations)
            unknown_annotation (int) : The encoding used in strand annotations for an unknown strand

        Returns:
            list(string): The genotype basecalls on the report strand
            The characters are A, C, G, T, or - for a no-call/null.

        Raises:
            ValueError: Number of SNPs or ref strand annotations not matched to entries in GTC file
        """
        genotypes = self.get_genotypes()
        if len(genotypes) != len(snps):
            raise ValueError(
                "The number of SNPs must match the number of loci in the GTC file")

        if len(genotypes) != len(strand_annotations):
            raise ValueError(
                "The number of reference strand annotations must match the number of loci in the GTC file")

        result = []
        for (genotype, snp, strand_annotation) in zip(genotypes, snps, strand_annotations):
            ab_genotype = code2genotype[genotype]
            a_nucleotide = snp[1]
            b_nucleotide = snp[-2]
            if a_nucleotide == "N" or b_nucleotide == "N" or strand_annotation == unknown_annotation or ab_genotype == "NC" or ab_genotype == "NULL":
                result.append("-")
            else:
                report_strand_nucleotides = []
                for ab_allele in ab_genotype:
                    nucleotide_allele = a_nucleotide if ab_allele == "A" else b_nucleotide
                    report_strand_nucleotides.append(
                        nucleotide_allele if strand_annotation == report_strand else complement(nucleotide_allele))
                result.append("".join(report_strand_nucleotides))
        return result

    def get_base_calls_plus_strand(self, snps, ref_strand_annotations):
        """
        Get base calls on plus strand of genomic reference. If you only see no-calls returned from this method,
        please verify that the reference strand annotations passed as argument are not unknown (RefStrand.Unknown)

        Args:
            snps (list(string)) : A list of string representing the snp on the design strand for the loci (e.g. [A/C])
            ref_strand_annotations (list(int)) : A list of strand annotations for the loci (e.g., RefStrand.Plus)

        Returns:
            list(string): The genotype basecalls on the report strand
            The characters are A, C, G, T, or - for a no-call/null.
        """
        return self.get_base_calls_generic(snps, ref_strand_annotations, RefStrand.Plus, RefStrand.Unknown)

    def get_base_calls_forward_strand(self, snps, forward_strand_annotations):
        """
        Get base calls on the forward strand.

        Args:
            snps (list(string)) : A list of string representing the snp on the design strand for the loci (e.g. [A/C])
            forward_strand_annotations (list(int)) : A list of strand annotations for the loci (e.g., SourceStrand.Forward)

        Returns:
            The genotype basecalls on the report strand as a list of strings.
            The characters are A, C, G, T, or - for a no-call/null.
        """
        return self.get_base_calls_generic(snps, forward_strand_annotations, SourceStrand.Forward, RefStrand.Unknown)

    def get_base_calls(self):
        """
        Returns:
            list(string): The genotype basecalls
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
            for idx in range(num_entries):
                if ploidy_type == 1:
                    result.append(gtc_handle.read(2))
                else:
                    byte_string = gtc_handle.read(2)
                    ab_genotype = code2genotype[genotypes[idx]]
                    if ab_genotype == "NC" or ab_genotype == "NULL":
                        result.append("-")
                    else:
                        top_genotype = "".join(
                            [byte_string[0] if allele == "A" else byte_string[1] for allele in ab_genotype])
                        result.append(top_genotype)
            return result

    def get_genotype_scores(self, offset=0, count=None):
        """
        Returns:
            list(numpy.float32): The genotype scores
        """
        return self.__get_generic_array_numpy(GenotypeCalls.__ID_GENOTYPE_SCORES, float32, offset, count)

    def get_scanner_data(self):
        """
        Returns:
            ScannerData: Information about scanner
        """
        return self.__get_generic(GenotypeCalls.__ID_SCANNER_DATA, ScannerData.read_scanner_data)

    def get_control_x_intensities(self):
        """
        Returns:
            list(int): The x intensities of control bead types
        """
        return self.__get_generic_array_numpy(GenotypeCalls.__ID_CONTROLS_X, uint16)

    def get_control_y_intensities(self):
        """
        Returns:
            list(int): The y intensities of control bead types
        """
        return self.__get_generic_array_numpy(GenotypeCalls.__ID_CONTROLS_Y, uint16)

    def get_raw_x_intensities(self, offset=0, count=None):
        """
        Returns:
            list(int): The raw x intensities of assay bead types
        """
        return self.__get_generic_array_numpy(GenotypeCalls.__ID_RAW_X, uint16, offset, count)

    def get_raw_y_intensities(self, offset=0, count=None):
        """
        Returns:
            list(int): The raw y intensities of assay bead types
        """
        return self.__get_generic_array_numpy(GenotypeCalls.__ID_RAW_Y, uint16, offset, count)

    def get_call_rate(self):
        """
        Returns:
            float: the call rate
        """
        return self.__get_generic(GenotypeCalls.__ID_CALL_RATE, read_float)

    def get_gender(self):
        """
        Returns:
            char: the gender
            M - Male, F - Female, U-Unknown
        """
        return self.__get_generic(GenotypeCalls.__ID_GENDER, read_char)

    def get_logr_dev(self):
        """
        Returns:
            float: the logR deviation
        """
        return self.__get_generic(GenotypeCalls.__ID_LOGR_DEV, read_float)

    def get_gc10(self):
        """
        Returns:
            float: the GC10 (GenCall score - 10th percentile)
        """
        return self.__get_generic(GenotypeCalls.__ID_GC10, read_float)

    def get_gc50(self):
        """
        Returns:
            float: the GC50 (GenCall score - 50th percentile)
        """
        return self.__get_generic(GenotypeCalls.__ID_GC50, read_float)

    def get_num_calls(self):
        """
        Returns:
            int: The number of calls
        """
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_GC50] + 4)
            return read_int(gtc_handle)

    def get_num_no_calls(self):
        """
        Returns:
            int: The number of no calls
        """
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_GC50] + 8)
            return read_int(gtc_handle)

    def get_num_intensity_only(self):
        """
        Returns:
            int: The number of intensity only SNPs
        """
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_GC50] + 12)
            return read_int(gtc_handle)

    def get_ballele_freqs(self, offset=0, count=None):
        """
        Returns:
            list(numpy.float32): The B allele frequencies
        """
        if self.version < 4:
            raise Exception(
                "B allele frequencies unavailable in GTC File version (" + str(self.version) + ")")
        return self.__get_generic_array_numpy(GenotypeCalls.__ID_B_ALLELE_FREQS, float32, offset, count)

    def get_logr_ratios(self, offset=0, count=None):
        """
        Returns:
            list(numpy.float32): The logR ratios
        """
        if self.version < 4:
            raise Exception(
                "LogR ratios unavailable in GTC File version (" + str(self.version) + ")")
        return self.__get_generic_array_numpy(GenotypeCalls.__ID_LOGR_RATIOS, float32, offset, count)

    def get_percentiles_x(self):
        """
        Returns:
            An array of length three representing 5th, 50th and 95th percentiles for
            x intensity
        """
        if self.version < 5:
            raise Exception(
                "Percentile intensities unavailable in GTC File version (" + str(self.version) + ")")
        result = []
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_PERCENTILES_X])
            result = []
            for idx in range(3):
                result.append(read_ushort(gtc_handle))
            return result

    def get_percentiles_y(self):
        """
        Returns:
            An array of length three representing 5th, 50th and 95th percentiles for
            y intensity
        """
        if self.version < 5:
            raise Exception(
                "Percentile intensities unavailable in GTC File version (" + str(self.version) + ")")
        result = []
        with open(self.filename, "rb") as gtc_handle:
            gtc_handle.seek(self.toc_table[GenotypeCalls.__ID_PERCENTILES_Y])
            result = []
            for idx in range(3):
                result.append(read_ushort(gtc_handle))
            return result

    def get_normalized_intensities(self, normalization_lookups):
        """
        Calculate and return the normalized intensities

        Args:
            normalization_lookups (list(int)): Map from each SNP to a normalization transform.
                                   This list can be obtained from the BeadPoolManifest object.

        Return:
            (float,float): The normalized intensities for the sample as a list of (x,y) float tuples
        """
        normalization_transforms = self.get_normalization_transforms()
        return [normalization_transforms[lookup].normalize_intensities(x_raw, y_raw) for (x_raw, y_raw, lookup) in zip(self.get_raw_x_intensities(), self.get_raw_y_intensities(), normalization_lookups)]

    def get_normalization_transforms(self):
        """
        Returns:
            List(NormalizationTransform): The normalization transforms used during genotyping 
        """
        return self.__get_generic_array(GenotypeCalls.__ID_NORMALIZATION_TRANSFORMS, NormalizationTransform.read_normalization_transform, 52, 0, None)

    def is_write_complete(self):
        """
        Check for last item written to GTC file to verify that write
        has successfully completed

        Args:
            None

        Returns:
            bool: Whether or not write is complete
        """
        if self.version == 3:
            try:
                self.get_num_intensity_only()
                return True
            except:
                return False
        elif self.version == 4:
            try:
                self.get_logr_ratios()
                return True
            except:
                return False
        elif self.version == 5:
            try:
                self.get_slide_identifier()
                return True
            except:
                return False
        else:
            raise Exception(
                "Unable to test for write completion on version " + str(self.version) + " GTC file")


class NormalizationTransform:
    def __init__(self, buffer):
        """
        Constructor

        Args:
            buffer (string): byte string containing binary data for transform

        Returns:
            NormalizationTransform
        """
        (self.version, ) = struct.unpack("<i", buffer[:4])
        (self.offset_x, self.offset_y, self.scale_x, self.scale_y,
         self.shear, self.theta,) = frombuffer(buffer[4:28], dtype=float32)

    @staticmethod
    def read_normalization_transform(handle):
        """
        Static helper function to read normalization transform from file handle

        Args:
            handle (file): File handle with position at start of normalization transform entry

        Returns:
            NormalizationTransform object
        """
        return NormalizationTransform(handle.read(52))

    @staticmethod
    def rect_to_polar(x, y):
        """
        Converts normalized x,y intensities to (pseudo) polar co-ordinates (R, theta)

        Args:
            x,y (float, float): Normalized x,y intensities for probe

        Returns:
            (float, float): (R,theta) polar values as tuple of floats
        """
        if x == 0 and y == 0:
            return (nan, nan)
        return (x + y, arctan2(y, x) * 2.0 / pi)

    def normalize_intensities(self, x, y, threshold=True):
        """
        Apply this normalization transform to raw intensities

        Args:
            x (int): Raw x intensities
            y (int): Raw y intensities

        Returns:
            (float, float): (xn, yn) normalized intensities as tuple of floats
        """
        if x <= 0 and y <= 0:
            return (nan, nan)

        tempx = x - self.offset_x
        tempy = y - self.offset_y

        tempx2 = cos(self.theta) * tempx + sin(self.theta) * tempy
        tempy2 = -sin(self.theta) * tempx + cos(self.theta) * tempy

        tempx3 = tempx2 - self.shear * tempy2
        tempy3 = tempy2

        xn = tempx3 / self.scale_x
        yn = tempy3 / self.scale_y

        if threshold:
            xn = 0 if 0 > xn else xn
            yn = 0 if 0 > yn else yn

        return (xn, yn)


class ScannerData(object):
    """
    Represents information about a scanner

    Attributes:
        name (string): scanner identifier
        pmt_green (int): gain setting (green channel)
        pmt_red (int): gain setting (red channel)
        version (string): version of scanner software
        user (string): user of the scanner software
    """

    def __init__(self, name, pmt_green, pmt_red, version, user):
        """
        Constructor

        Args:
            name (string): scanner identifier
            pmt_green (int): gain setting (green channel)
            pmt_red (int): gain setting (red channel)
            version (string): version of scanner software
            user (string): user of the scanner software

        Returns:
            ScannerData
        """
        self.name = name
        self.pmt_green = pmt_green
        self.pmt_red = pmt_red
        self.version = version
        self.user = user

    @staticmethod
    def read_scanner_data(handle):
        """
        Helper function to parse ScannerData object from file handle.

        Args:
            handle (file): File handle

        Returns:
            ScannerData
        """
        name = read_string(handle)
        pmt_green = read_int(handle)
        pmt_red = read_int(handle)
        scanner_version = read_string(handle)
        imaging_user = read_string(handle)
        return ScannerData(name, pmt_green, pmt_red, scanner_version, imaging_user)
