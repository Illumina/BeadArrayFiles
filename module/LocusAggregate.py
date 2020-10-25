class Loader(object):
    """
    A function object which handles loading a contiguous slice of loci
    data for a given sample
    """
    def __init__(self, locus_offset, loci_buffer_size, normalization_lookups):
        """
        Constructor

        Args:
            locus_offset (int): Start loading data at this locus index
            loci_buffer_size (int): Number of loci to load (starting at offset)
            normalization_lookups (list(int)): List of normalization lookups (for all loci in sample data, not just slice).

        Returns:
            Loader: The new Loader object
        """
        self.locus_offset = locus_offset
        self.loci_buffer_size = loci_buffer_size
        self.normalization_lookups = normalization_lookups


    def __call__(self, sample_data):
        """
        Based on the configuration of the Loader function, load a slice of loci data
        from the given sample. If the sample data is not version 4 or above, the log_r_ratios
        and b_allele_freqs will be populated with an array of "None"

        Args:
            sample_data (GenotypeCalls): The sample data to provide loci information

        Returns:
            LocusAggregate: slice of locus data from this sample
        """
        locus_offset = self.locus_offset
        loci_buffer_size = self.loci_buffer_size
        locus_aggregate = LocusAggregate()
        locus_aggregate.genotypes = sample_data.get_genotypes(
            locus_offset, loci_buffer_size)
        locus_aggregate.scores = sample_data.get_genotype_scores(
            locus_offset, loci_buffer_size)
        locus_aggregate.b_allele_freqs = sample_data.get_ballele_freqs(
            locus_offset, loci_buffer_size) if sample_data.version >= 4 else [None] * loci_buffer_size
        locus_aggregate.log_r_ratios = sample_data.get_logr_ratios(
            locus_offset, loci_buffer_size) if sample_data.version >= 4 else [None] * loci_buffer_size
        locus_aggregate.x_intensities = sample_data.get_raw_x_intensities(
            locus_offset, loci_buffer_size)
        locus_aggregate.y_intensities = sample_data.get_raw_y_intensities(
            locus_offset, loci_buffer_size)
        transforms = sample_data.get_normalization_transforms()
        locus_aggregate.transforms = [
            transforms[idx] for idx in self.normalization_lookups[locus_offset:locus_offset + loci_buffer_size]]
        return locus_aggregate


class GenerateLocusAggregate(object):
    """
    Function object to generate aggregate data for a given
    locus
    """

    def __init__(self, buffer, relative_offset):
        """
        Constructor

        Args:
            buffer (list(LocusAggregate)): List where each element is slice of locus data from single sample
            relative_offset (int): Offset that was used when loading buffer

        Returns:
            GenerateLocusAggregate
        """
        self.buffer = buffer
        self.relative_offset = relative_offset

    def __call__(self, locus_idx):
        """
        Create a new LocusAggregate representing data from a single locus across all samples

        Args:
            locus_idx (int): Global locus index (will be automatically adjusted by relative offset in constructor)

        Returns:
            LocusAggregate : Data for a single locus aggregated across all samples
        """
        locus_aggregate = LocusAggregate()
        relative_locus_idx = locus_idx - self.relative_offset
        for sample_buffer in self.buffer:
            locus_aggregate.genotypes.append(
                sample_buffer.genotypes[relative_locus_idx])
            locus_aggregate.scores.append(
                sample_buffer.scores[relative_locus_idx])
            locus_aggregate.b_allele_freqs.append(
                sample_buffer.b_allele_freqs[relative_locus_idx])
            locus_aggregate.log_r_ratios.append(
                sample_buffer.log_r_ratios[relative_locus_idx])
            locus_aggregate.x_intensities.append(
                sample_buffer.x_intensities[relative_locus_idx])
            locus_aggregate.y_intensities.append(
                sample_buffer.y_intensities[relative_locus_idx])
            locus_aggregate.transforms.append(
                sample_buffer.transforms[relative_locus_idx])
        return locus_aggregate

class LocusAggregate(object):
    """
    Class to contain aggregated data for a single locus
    across many samples. For each attribute list, the individual lists
    are those values as a slice across the specified samples

    Attributes:
        genotypes (list(int)): List of integer genotypes
        scores (list(float)): List of Gencall scores
        b_allele_freqs (list(float)): List of b allele frequencies
        log_r_ratios (list(float)): List of log R ratios
        x_intensities (list(int)): List of raw x intensities
        y_intensities (list(int)): List of raw y intensities
        transforms (list(NormalizationTransform)): List of normalization transforms
    """

    def __init__(self):
        """
        Constructor
            All attributes initialized to be empty
        """
        self.genotypes = []
        self.scores = []
        self.b_allele_freqs = []
        self.log_r_ratios = []
        self.x_intensities = []
        self.y_intensities = []
        self.transforms = []

    @staticmethod
    def load_buffer(samples, locus_offset, loci_buffer_size, normalization_lookups):
        """
        Load a subset of continguous loci across a number of samples

        Args:
            samples (list(GenotypeCalls)): The samples to aggregate
            locus_offset (int): Index of the first locus to load
            loci_buffer_size (int): Total number of loci to load for each sample:
            process_pool (multiprocessing.Pool): Process pool for running parallel operations

        Returns:
            list(LocusAggregate): A locus aggregate object for each sample
        """
        return map(Loader(locus_offset, loci_buffer_size, normalization_lookups), samples)

    @staticmethod
    def group_loci(loci, loci_batch_size):
        """
        Group loci indices such that the first and last loci
        are not separated by more than the batch size

        Args:
            loci (list(int)): List of loci indices (must be sorted)
            loci_batch_size (int): Maximum difference to allow between first and last grouped index

        Yields:
            list(int) : Grouped loci indexes
        """
        current_indices = []
        for locus in loci:
            if len(current_indices) != 0:
                if current_indices[-1] > locus:
                    raise ValueError("Loci indices are not sorted")
                if locus - loci_batch_size >= current_indices[0]:
                    yield current_indices
                    current_indices = []
            current_indices.append(locus)
        if len(current_indices) > 0:
            yield current_indices

    @staticmethod
    def aggregate_samples(samples, loci, callback, normalization_lookups, bin_size=100000000):
        """
        Generate LocusAggregate information from a collection of samples. Will call the callback
        function for a LocusAggregate object for each specified locus index and yield the result.

        Args:
            samples (list(GenotypeCalls)): The samples to aggregate for each locus
            loci (iter(int)): Enumerates the loci indices of interest (must be sorted in ascending order)
            callback (func): A function that takes a LocusAggregate and return a new result
            bin_size (int): Used to determine how much data will be loaded into memory at one time. Larger bin size will use more memory and (generally) run faster. This bin_size already accounts for how many samples are being handled.

        Yields:
            Result of callback function
        """
        # figure out how many loci to load at once
        loci_batch_size = int(bin_size / float(len(list(samples)))) + 1

        for loci_group in LocusAggregate.group_loci(loci, loci_batch_size):

            # read in the buffer for this group of loci
            buffer = LocusAggregate.load_buffer(
                samples, loci_group[0], loci_group[-1] - loci_group[0] + 1, normalization_lookups)

            # generate corresponding locus aggregates
            aggregates = map(GenerateLocusAggregate(
                buffer, loci_group[0]), loci_group)

            for result in map(callback, aggregates):
                yield result
