from .BeadArrayUtility import read_int, read_string, read_byte

class BeadPoolManifest(object):
    """
    Class for parsing binary (BPM) manifest file.

    Attributes:
        names (list(strings): Names of loci from manifest
        snps (list(strings): SNP values of loci from manifest
        chroms (list(strings): Chromosome values for loci
        map_infos (list(int)): Map infor values for loci
        addresses (list(int)): AddressA IDs of loci from manifest
        normalization_lookups (list(int)): Normalization lookups from manifest. This indexes into
                                         list of normalization transforms read from GTC file
        ref_strands (list(int)): Reference strand annotation for loci (see RefStrand class)
        source_strands (list(int)): Source strand annotations for loci (see SourceStrand class)
        num_loci (int): Number of loci in manifest
        manifest_name (string): Name of manifest
        control_config (string): Control description from manifest
    """

    def __init__(self, filename):
        """
        Constructor

        Args:
            filename (string): Locations of BPM (bead pool manifest) file

        Returns:
            BeadPoolManifest
        """
        self.names = []
        self.snps = []
        self.chroms = []
        self.map_infos = []
        self.addresses = []
        self.normalization_ids = []
        self.assay_types = []
        self.normalization_lookups = []
        self.ref_strands = []
        self.source_strands = []
        self.num_loci = 0
        self.manifest_name = ""
        self.control_config = ""
        self.__parse_file(filename)

    def __parse_file(self, manifest_file):
        """
        Helper function to initialize this object from a file.

        Args:
            manifest_file (string): Location of BPM (bead pool manifest) file

        Returns:
            None

        Raises:
            Exception: Unsupported or unknown BPM version
            Exception: Manifest format error
        """
        with open(manifest_file, "rb") as manifest_handle:
            header = manifest_handle.read(3)
            header = header.decode("utf-8")
            if len(header) != 3 or header != "BPM":
                raise Exception("Invalid BPM format")
            version = read_byte(manifest_handle)
            if version != 1:
                raise Exception(
                    "Unknown BPM version (" + str(ord(version)) + ")")

            version = read_int(manifest_handle)
            version_flag = 0x1000
            if version & version_flag == version_flag:
                version = version ^ version_flag
            if version > 5 or version < 3:
                raise Exception(
                    "Unsupported BPM version (" + str(version) + ")")
            self.manifest_name = read_string(manifest_handle)

            if version > 1:
                self.control_config = read_string(manifest_handle)

            self.num_loci = read_int(manifest_handle)
            manifest_handle.seek(4 * self.num_loci, 1)
            name_lookup = {}
            for idx in range(self.num_loci):
                self.names.append(read_string(manifest_handle))
                name_lookup[self.names[-1]] = idx

            for idx in range(self.num_loci):
                normalization_id = read_byte(manifest_handle)
                if normalization_id >= 100:
                    raise Exception(
                        "Manifest format error: read invalid normalization ID")
                self.normalization_ids.append(normalization_id)

            self.assay_types = [0] * self.num_loci
            self.addresses = [0] * self.num_loci
            self.snps = [""] * self.num_loci
            self.chroms = [""] * self.num_loci
            self.map_infos = [0] * self.num_loci
            self.ref_strands = [RefStrand.Unknown] * self.num_loci
            self.source_strands = [SourceStrand.Unknown] * self.num_loci
            for idx in range(self.num_loci):
                locus_entry = LocusEntry(manifest_handle)
                self.assay_types[name_lookup[locus_entry.name]] = locus_entry.assay_type
                self.addresses[name_lookup[locus_entry.name]] = locus_entry.address_a
                self.snps[name_lookup[locus_entry.name]] = locus_entry.snp
                self.chroms[name_lookup[locus_entry.name]] = locus_entry.chrom
                self.map_infos[name_lookup[locus_entry.name]] = locus_entry.map_info
                self.ref_strands[name_lookup[locus_entry.name]] = locus_entry.ref_strand
                self.source_strands[name_lookup[locus_entry.name]] = locus_entry.source_strand

            if len(self.normalization_ids) != len(self.assay_types):
                raise Exception(
                    "Manifest format error: read invalid number of assay entries")

            all_norm_ids = set()
            for locus_idx in range(self.num_loci):
                self.normalization_ids[locus_idx] += 100 * \
                    self.assay_types[locus_idx]
                # To mimic the byte-wrapping behavior from GenomeStudio, AutoCall, IAAP take the mod of 256
                self.normalization_ids[locus_idx] %= 256
                all_norm_ids.add(self.normalization_ids[locus_idx])
            sorted_norm_ids = sorted(all_norm_ids)
            lookup_dictionary = {}
            for idx in range(len(sorted_norm_ids)):
                lookup_dictionary[sorted_norm_ids[idx]] = idx
            self.normalization_lookups = [
                lookup_dictionary[normalization_id] for normalization_id in self.normalization_ids]


class SourceStrand(object):
    Unknown = 0
    Forward = 1
    Reverse = 2

    @staticmethod
    def to_string(source_strand):
        """Get an integer representation of source strand annotation

        Args:
            source_strand (str) : string representation of source strand annotation (e.g., "F")

        Returns:
            int : int representation of source strand annotation (e.g. SourceStrand.Forward)

        Raises:
            ValueError: Unexpected value for source strand
        """
        if source_strand == SourceStrand.Unknown:
            return "U"
        elif source_strand == SourceStrand.Forward:
            return "F"
        elif source_strand == SourceStrand.Reverse:
            return "R"
        else:
            raise ValueError(
                "Unexpected value for source strand " + source_strand)

    @staticmethod
    def from_string(source_strand):
        """
        Get a string representation of source strand annotation

        Args:
            source_strand (int) : int representation of source strand (e.g., SourceStrand.Forward)

        Returns:
            str : string representation of source strand annotation

        Raises:
            ValueError: Unexpected value for source strand
        """
        if source_strand == "U" or source_strand == "":
            return SourceStrand.Unknown
        if source_strand == "F":
            return SourceStrand.Forward
        elif source_strand == "R":
            return SourceStrand.Reverse
        else:
            raise ValueError(
                "Unexpected value for source strand " + source_strand)

class RefStrand(object):
    Unknown = 0
    Plus = 1
    Minus = 2

    @staticmethod
    def to_string(ref_strand):
        """
        Get a string reprensetation of ref strand annotation

        Args:
            ref_strand (int) : int representation of ref strand (e.g., RefStrand.Plus)

        Returns:
            str : string representation of reference strand annotation

        Raises:
            ValueError: Unexpected value for reference strand
        """
        if ref_strand == RefStrand.Unknown:
            return "U"
        elif ref_strand == RefStrand.Plus:
            return "+"
        elif ref_strand == RefStrand.Minus:
            return "-"
        else:
            raise ValueError(
                "Unexpected value for reference strand " + ref_strand)

    @staticmethod
    def from_string(ref_strand):
        """Get an integer representation of ref strand annotation

        Args:
            ref_strand (str) : string representation of reference strand annotation (e.g., "+")

        Returns:
            int : int representation of reference strand annotation (e.g. RefStrand.Plus)

        Raises:
            ValueError: Unexpected value for reference strand
        """
        if ref_strand == "U" or ref_strand == "":
            return RefStrand.Unknown
        if ref_strand == "+":
            return RefStrand.Plus
        elif ref_strand == "-":
            return RefStrand.Minus
        else:
            raise ValueError(
                "Unexpected value for reference strand " + ref_strand)

class LocusEntry(object):
    """
    Helper class representing a locus entry within a bead pool manifest. Current only support version
    6,7, and 8.

    Attributes:
        ilmn_id (string) : IlmnID (probe identifier) of locus
        name (string): Name (variant identifier) of locus
        snp (string) : SNP value for locus (e.g., [A/C])
        chrom (string) : Chromosome for the locus (e.g., XY)
        map_info (int) : Mapping location of locus
        assay_type (int) : Identifies type of assay (0 - Infinium II , 1 - Infinium I (A/T), 2 - Infinium I (G/C)
        address_a (int) : AddressA ID of locus
        address_b (int) : AddressB ID of locus (0 if none)
        ref_strand (int) : See RefStrand class
        source_strand (int) : See SourceStrand class
    """

    def __init__(self, handle):
        """
        Constructor

        Args:
            handle (file):  File handle at start of locus entry record

        Returns:
            LocusEntry
        """
        self.ilmn_id = ""
        self.name = ""
        self.snp = ""
        self.chrom = ""
        self.map_info = -1
        self.assay_type = -1
        self.address_a = -1
        self.address_b = -1
        self.ref_strand = RefStrand.Unknown
        self.source_strand = SourceStrand.Unknown
        self.__parse_file(handle)

    def __parse_file(self, handle):
        """
        Helper function to initialize this object from a file handle

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
            raise Exception(
                "Manifest format error: unknown version for locus entry (" + str(version) + ")")

    def __parse_locus_version_6(self, handle):
        """
        Helper function to parse version 6 locus entry

        Args:
            handle (file): File handle at start of locus entry record

        Returns:
            None

        Raises:
            Exception: Manifest format error
        """
        self.ilmn_id = read_string(handle)
        self.source_strand = SourceStrand.from_string(
            self.ilmn_id.split("_")[-2])
        self.name = read_string(handle)
        for idx in range(3):
            read_string(handle)
        handle.read(4)
        for idx in range(2):
            read_string(handle)
        self.snp = read_string(handle)
        self.chrom = read_string(handle)
        for idx in range(2):
            read_string(handle)
        self.map_info = int(read_string(handle))
        for idx in range(2):
            read_string(handle)
        self.address_a = read_int(handle)
        self.address_b = read_int(handle)
        for idx in range(7):
            read_string(handle)
        handle.read(3)
        self.assay_type = read_byte(handle)
        if self.assay_type not in [0, 1, 2]:
            raise Exception(
                "Format error in reading assay type from locus entry")
        if self.address_b == 0:
            if self.assay_type != 0:
                raise Exception(
                    "Manifest format error: Assay type is inconsistent with address B")
        else:
            if self.assay_type == 0:
                raise Exception(
                    "Manifest format error: Assay type is inconsistent with address B")

    def __parse_locus_version_7(self, handle):
        """
        Helper function to parse version 7 locus entry

        Args:
            handle (file): File handle at start of locus entry record

        Returns:
            None
        """
        self.__parse_locus_version_6(handle)
        handle.read(4 * 4)

    def __parse_locus_version_8(self, handle):
        """
        Helper function to parse version 8 locus entry

        Args:
            handle (file): File handle at start of locus entry record

        Returns:
            None
        """
        self.__parse_locus_version_7(handle)
        self.ref_strand = RefStrand.from_string(read_string(handle))
