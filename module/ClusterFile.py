from .BeadArrayUtility import read_int, read_string, read_byte, read_float

class ClusterFile(object):
    """
    Represents an EGT cluster file

    Attributes:
        gencall_version (string): The GenCall version
        cluster_version (string): The clustering algorithm version
        call_version (string): The genotyping algorithm version
        normalization_version (string): The normalization algorithm version
        date_created (string): The date the cluster file was created (e.g., 3/9/2017 2:18:30 PM)
        manifest_name (string): The manifest name used to build this cluster file
    """

    def __init__(self, gencall_version, cluster_version, call_version, normalization_version, date_created, manifest_name):
        """
        Constructor

        Args:
            gencall_version (string): The GenCall version
            cluster_version (string): The clustering algorithm version
            call_version (string): The genotyping algorithm version
            normalization_version (string): The normalization algorithm version
            date_created (string): The date the cluster file was created (e.g., 3/9/2017 2:18:30 PM)
            manifest_name (string): The manifest name used to build this cluster file

        Returns:
            ClusterFile
        """
        self.gencall_version = gencall_version
        self.cluster_version = cluster_version
        self.call_version = call_version
        self.normalization_version = normalization_version
        self.date_created = date_created
        self.manifest_name = manifest_name
        self.name2cluster_record = {}

    def add_record(self, name, record):
        """
        Add a new record to the cluster file

        Args:
            name (string): Locus name
            record (ClusterRecord): Record for the locus
        """
        self.name2cluster_record[name] = record

    def get_record(self, name):
        """
        Get the record for a locus

        Args:
            name (string): Locus name

        Returns:
            LocusRecord: The record associated with the locus

        Raises:
            KeyError: Locus name not present in cluster file
        """
        return self.name2cluster_record[name]

    @staticmethod
    def read_array(handle, num_entries, read_record):
        """
        Helper method for reading an array of data from a file handle

        Args:
            handle (file): File handle
            num_entries (int): Number of entries to read
            read_record (func(file)): Function (taking file handle as argument) to read a single entry

        Returns:
            list(value): A list of type value returned by read_record function
        """
        result = []
        for idx in range(num_entries):
            result.append(read_record(handle))
        return result

    @staticmethod
    def read_cluster_file(handle):
        """
        Read a cluster file

        Args:
            file: EGT cluster file handle

        Returns:
            ClusterFile

        Raises:
            Exception: Incompatible cluster file format
        """
        version = read_int(handle)
        if version != 3:
            raise Exception("Cluster file version " +
                            str(version) + " not supported")

        gencall_version = read_string(handle)
        cluster_version = read_string(handle)
        call_version = read_string(handle)
        normalization_version = read_string(handle)
        date_created = read_string(handle)

        is_wgt = read_byte(handle) == 1
        if not is_wgt:
            raise Exception("Only WGT cluster file version supported")

        manifest_name = read_string(handle)

        result = ClusterFile(gencall_version, cluster_version, call_version,
                             normalization_version, date_created, manifest_name)
        data_block_version = read_int(handle)
        if data_block_version not in [8, 9]:
            raise Exception("Data block version in cluster file " +
                            str(data_block_version) + " not  supported")
        # opa
        _ = read_string(handle)

        num_records = read_int(handle)
        cluster_records = ClusterFile.read_array(
            handle, num_records, lambda handle: ClusterRecord.read_record(handle, data_block_version))
        cluster_scores = ClusterFile.read_array(
            handle, num_records, ClusterScore.read_record)

        # genotypes
        _ = ClusterFile.read_array(handle, num_records, read_string)

        loci_names = ClusterFile.read_array(handle, num_records, read_string)
        addresses = ClusterFile.read_array(handle, num_records, read_int)

        # cluster counts
        cluster_counts = []
        for idx in range(num_records):
            # 3 corresponds to number genotypes (AA, AB, BB)
            cluster_counts.append(ClusterFile.read_array(handle, 3, read_int))

        for (cluster_record, count_record) in zip(cluster_records, cluster_counts):
            assert cluster_record.aa_cluster_stats.N == count_record[0]
            assert cluster_record.ab_cluster_stats.N == count_record[1]
            assert cluster_record.bb_cluster_stats.N == count_record[2]

        for (locus_name, address, cluster_record, cluster_score) in zip(loci_names, addresses, cluster_records, cluster_scores):
            cluster_record.address = address
            cluster_record.cluster_score = cluster_score
            result.add_record(locus_name, cluster_record)

        return result

class ClusterRecord(object):
    """
    Store clustering information for a single locus

    Attributes:
        aa_cluster_stats (ClusterRecord): Describes AA genotype cluster
        ab_cluster_stats (ClusterRecord): Describes AB genotype cluster
        bb_cluster_stats (ClusterRecord): Describes BB genotype cluster
        intensity_threshold (float): Intensity threshold for no-calll
        cluster_score (ClusterScore): Various scores for cluster
        address (int): Bead type identifier for probe A
    """
    def __init__(self, aa_cluster_stats, ab_cluster_stats, bb_cluster_stats, intensity_threshold, cluster_score, address):
        """
        Constructor

        Args:
            aa_cluster_stats (ClusterRecord): Describes AA genotype cluster
            ab_cluster_stats (ClusterRecord): Describes AB genotype cluster
            bb_cluster_stats (ClusterRecord): Describes BB genotype cluster
            intensity_threshold (float): Intensity threshold for no-calll
            cluster_score (ClusterScore): Various scores for cluster
            address (int): Bead type identifier for probe A

        Returns:
            ClusterRecord
        """

        self.aa_cluster_stats = aa_cluster_stats
        self.ab_cluster_stats = ab_cluster_stats
        self.bb_cluster_stats = bb_cluster_stats
        self.intensity_threshold = intensity_threshold
        self.cluster_score = cluster_score
        self.address = address

    @staticmethod
    def read_record(handle, version):
        """
        Read a cluster record from the file handle

        Args:
            handle (file): The file handle
            version (int): The cluster record version (from header)

        Returns:
            ClusterRecord: Result will not be populated with either address or scores (read separately)

        Raises:
            Exception : Unsupported cluster record version
        """
        (aa_n, ab_n, bb_n) = ClusterFile.read_array(handle, 3, read_int)
        (aa_r_dev, ab_r_dev, bb_r_dev) = ClusterFile.read_array(
            handle, 3, read_float)
        (aa_r_mean, ab_r_mean, bb_r_mean) = ClusterFile.read_array(
            handle, 3, read_float)
        (aa_theta_dev, ab_theta_dev, bb_theta_dev) = ClusterFile.read_array(
            handle, 3, read_float)
        (aa_theta_mean, ab_theta_mean, bb_theta_mean) = ClusterFile.read_array(
            handle, 3, read_float)

        if version == 9:
            intensity_threshold = read_float(handle)
        elif version == 8:
            _ = read_float(handle)
            intensity_threshold = 0
        else:
            raise Exception(
                "Unsupported cluster record version " + str(version))

        # read through unused fields
        for idx in range(14):
            _ = read_float(handle)

        aa_cluster_stats = ClusterStats(
            aa_theta_mean, aa_theta_dev, aa_r_mean, aa_r_dev, aa_n)
        ab_cluster_stats = ClusterStats(
            ab_theta_mean, ab_theta_dev, ab_r_mean, ab_r_dev, ab_n)
        bb_cluster_stats = ClusterStats(
            bb_theta_mean, bb_theta_dev, bb_r_mean, bb_r_dev, bb_n)

        return ClusterRecord(aa_cluster_stats, ab_cluster_stats, bb_cluster_stats, intensity_threshold, None, None)


class ClusterScore(object):
    """
    All scores for a given locus clustering.

    Attributes:
        cluster_separation (float): A score measure the separation between genotype clusters
        total_score (float): The GenTrain score
        original_score (float): The original score before editing this cluster
        edited (bool): Whether this cluster has been manually manipulated
    """

    def __init__(self, cluster_separation, total_score, original_score, edited):
        """
        Constructor

        Args:
            cluster_separation (float): A score measure the separation between genotype clusters
            total_score (float): The GenTrain score
            original_score (float): The original score before editing this cluster
            edited (bool): Whether this cluster has been manually manipulated

        Returns:
            ClusterScore
        """
        self.cluster_separation = cluster_separation
        self.total_score = total_score
        self.original_score = original_score
        self.edited = edited

    @staticmethod
    def read_record(handle):
        """
        Read a ClusterScore from a file handle

        Args:
            handle (file): The file handle

        Returns:
            ClusterScore
        """
        cluster_separation = read_float(handle)
        total_score = read_float(handle)
        original_score = read_float(handle)
        edited = read_byte(handle) != 0
        return ClusterScore(cluster_separation, total_score, original_score, edited)


class ClusterStats(object):
    """
    Represents statistics for a single genotype cluster.

    Attributes:
        theta_mean (float): Theta mean value
        theta_dev (float): Theta std devation value
        r_mean (float): R (intensity) mean value
        r_dev (float): R (intensity) std deviation value
        N (int): Number of samples assigned to cluster during training
    """

    def __init__(self, theta_mean, theta_dev, r_mean, r_dev, N):
        """
        Constructor

        Args:
            theta_mean (float): Theta mean value
            theta_dev (float): Theta std devation value
            r_mean (float): R (intensity) mean value
            r_dev (float): R (intensity) std deviation value
            N (int): Number of samples assigned to cluster during training

        Returns:
            ClusterStats
        """
        self.theta_mean = theta_mean
        self.theta_dev = theta_dev
        self.r_mean = r_mean
        self.r_dev = r_dev
        self.N = N
