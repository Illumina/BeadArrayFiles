import sys
import os
import argparse
import logging
import pandas as pd
from math import isnan
from numpy import percentile
from scipy.stats import norm


from IlluminaBeadArrayFiles import LocusAggregate, BeadPoolManifest, GenotypeCalls, ClusterFile, RefStrand

nan = float("nan")


class LocusSummary(object):
	def __init__(self, genotype_counts, score_stats):
		self.genotype_counts = genotype_counts
		self.score_stats = score_stats

class GenotypeCounts(object):
	"""
	Summarize information about genotype counts for diploid genotyping counting
	"""

	def __init__(self, genotypes):
		self.no_calls = 0
		self.aa_count = 0
		self.ab_count = 0
		self.bb_count = 0


		for genotype in genotypes:
			if genotype == 0:
				self.no_calls += 1
			elif genotype == 1:
				self.aa_count += 1
			elif genotype == 2:
				self.ab_count += 1
			elif genotype == 3:
				self.bb_count += 1

	def get_num_calls(self):
		"""
		Get the number of calls (i.e., not no-calls)

		Returns:
			int: The number of calls
		"""
		return self.aa_count + self.ab_count + self.bb_count

	def get_call_frequency(self):
		"""
		Get the call rate

		Returns:
			float: The frequency of calls
		"""
		num_calls = self.get_num_calls()
		return num_calls / float(num_calls + self.no_calls) if num_calls + self.no_calls > 0 else nan

	def get_aa_frequency(self):
		"""
		Frequency of AA genotype (as fraction of all calls)

		Returns:
			float: AA genotype frequency
		"""
		return self.aa_count / float(self.get_num_calls()) if self.get_num_calls() > 0 else nan

	def get_ab_frequency(self):
		"""
		Frequency of AB genotype (as fraction of all calls)

		Returns:
			float: AB genotype frequency
		"""
		return self.ab_count / float(self.get_num_calls()) if self.get_num_calls() > 0 else nan

	def get_bb_frequency(self):
		"""
		Frequency of BB genotype (as fraction of all calls)

		Returns:
			float: BB genotype frequency
		"""
		return self.bb_count / float(self.get_num_calls()) if self.get_num_calls() > 0 else nan

	def get_minor_frequency(self):
		"""
		Comoputes and return the minor allele frequency. If no calls, will be NaN

		Returns:
			float
		"""
		a_allele_count = self.aa_count * 2 + self.ab_count
		a_frequency = a_allele_count / \
			float(2 * self.get_num_calls()) if self.get_num_calls() > 0 else nan
		return min(a_frequency, 1.0 - a_frequency) if not isnan(a_frequency) else nan

	def compute_hardy_weinberg(self):
		"""
		Computes and returns statistics related to HW equilibrium

		Returns:
			(float, float): Het excess and ChiSq 100 statistics, respectively
		"""
		num_calls = self.get_num_calls()
		if num_calls == 0:
			return (0.0, 0.0)

		if self.aa_count + self.ab_count == 0 or self.ab_count + self.bb_count == 0:
			return (1.0, 0.0)

		num_calls = float(num_calls)

		q = self.get_minor_frequency()
		p = 1 - q

		temp = 0.013 / q
		k = temp * temp * temp * temp
		dh = ((self.ab_count / num_calls + k) / (2 * p * q + k)) - 1
		if dh < 0:
			hw = (2 * norm.cdf(dh, 0, 1 / 10.0))
		else:
			hw = (2 * (1 - norm.cdf(dh, 0, 1 / 10.0)))

		return (hw, dh)


class ScoreStatistics(object):
	"""
	Capture statistics related to the gencall score distribution

	Attributes:
		gc_10 : 10th percentile of Gencall score distribution
		gc_50 : 50th percentile of Gencall score distribution
	"""

	def __init__(self, scores, genotypes):
		"""
		Create new ScoreStatistics object

		Args:
			score (list(float)): A list of gencall scores
			genotypes (list(int)): A list of genotypes

		Returns:
			ScoreStatistics
		"""
		called_scores = sorted([score for (score, genotype) in zip(scores, genotypes) if genotype != 0])
		self.gc_10 = ScoreStatistics.percentile(called_scores, 10)
		self.gc_50 = ScoreStatistics.percentile(called_scores, 50)

	@staticmethod
	def percentile(scores, percentile):
		"""
		Percentile as calculated in GenomeStudio

		Args:
			scores (list(float)): list of scores (typically for called genotypes)
			percentile (int): percentile to calculate
		
		Returns:
			float
		"""
		num_scores = len(scores)
		if num_scores == 0:
			return nan

		idx = int(num_scores*percentile/100)
		fractional_index = num_scores*percentile/100.0 - idx
		if fractional_index < 0.5 :
			idx -= 1

		if idx < 0:
			return scores[0]

		if idx >= num_scores - 1:
			return scores[-1]

		x1 = 100 * (idx + 0.5)/float(num_scores)
		x2 = 100 * (idx + 1 + 0.5)/float(num_scores)
		y1 = float(scores[idx])
		y2 = float(scores[idx+1])

		return y1 + (y2 - y1) / (x2 - x1) * (percentile - x1)



def summarize_locus(snp_wise_genotypes,snp_wise_scores):
	"""
	Generate a locus summary based on aggregated locus information

	Args:
		LocusAggregate : Aggregated information for a locus
	
	Returns
		LocusSummary
	"""
	genotype_counts = GenotypeCounts(snp_wise_genotypes)
	score_stats = ScoreStatistics(snp_wise_scores, snp_wise_genotypes)
	return LocusSummary(genotype_counts, score_stats)


def get_logger():
	# set up log file
	# create logger
	logger = logging.getLogger('Locus Summary Report')
	logger.setLevel(logging.DEBUG)

	# create console handler and set level to debug
	handler = logging.StreamHandler()
	handler.setLevel(logging.INFO)

	# create formatter
	formatter = logging.Formatter(
		'%(asctime)s - %(name)s - %(levelname)s - %(message)s')

	# add formatter to ch
	handler.setFormatter(formatter)
	logger.addHandler(handler)
	return logger


def driver(gtc_dir, manifest_filename, cluster_filename, output_filename, project_name, delim, logger):
	logger.info("Reading cluster file")
	with open(cluster_filename, "rb") as cluster_handle:
		egt = ClusterFile.read_cluster_file(cluster_handle)

	logger.info("Reading manifest file")
	bpm = BeadPoolManifest(manifest_filename)
	samples = []

	logger.info("Initializing genotype data")
	gtc_files = []
	for gtc_file in os.listdir(gtc_dir):
		if gtc_file.endswith(".gtc"):
			gtc_files.append(os.path.join(gtc_dir, gtc_file))

	samples = map(GenotypeCalls, gtc_files)


	ls_genotypes = []
	ls_genotype_scores = []
	ls_sample_names = []
	ls_snps = bpm.names
	for sample in samples:
		genotypes = sample.get_genotypes()
		assert len(genotypes) == len(bpm.names)
		ls_genotypes.append(genotypes)
		ls_genotype_scores.append(sample.get_genotype_scores())
		ls_sample_names.append(sample.get_sample_name())
			
	logger.info("Generating report")
	loci = range(len(bpm.normalization_lookups))
	row = 0
	with open(output_filename, "w") as output_handle:
		output_handle.write("Locus Summary on " + os.path.abspath(output_filename) + "\n")
		header = [""]
		header.append("# LOCI = {}".format(len(loci)))
		header.append("# DNAs = {}".format(len(gtc_files)))
		header.append("ProjectName = {}".format(project_name))
		header.append("GenCall Version = {}".format(egt.gencall_version))
		header.append("Low GenCall Score Cutoff = NaN")

		output_handle.write(delim.join(header) + "\n")

		output_handle.write(delim.join("Row,Locus_Name,Illumicode_Name,#No_Calls,#Calls,Call_Freq,A/A_Freq,A/B_Freq,B/B_Freq,Minor_Freq,Gentrain_Score,50%_GC_Score,10%_GC_Score,Het_Excess_Freq,ChiTest_P100,Cluster_Sep,AA_T_Mean,AA_T_Std,AB_T_Mean,AB_T_Std,BB_T_Mean,BB_T_Std,AA_R_Mean,AA_R_Std,AB_R_Mean,AB_R_Std,BB_R_Mean,BB_R_Std,Plus/Minus Strand".split(",")) + "\n")
		for i in range(0,len(ls_snps)):
			row += 1
			snp_wise_genotypes = [item[i] for item in ls_genotypes]
			snp_wise_scores = [item[i] for item in ls_genotype_scores]
			locus_summary = summarize_locus(snp_wise_genotypes,snp_wise_scores)
			cluster_record = egt.get_record(ls_snps[i])
			row_data = []
			row_data.append(row)
			row_data.append(ls_snps[i])
			row_data.append(cluster_record.address)
			row_data.append(locus_summary.genotype_counts.no_calls)
			row_data.append(locus_summary.genotype_counts.get_num_calls())
			row_data.append(locus_summary.genotype_counts.get_call_frequency())
			row_data.append(locus_summary.genotype_counts.get_aa_frequency())
			row_data.append(locus_summary.genotype_counts.get_ab_frequency())
			row_data.append(locus_summary.genotype_counts.get_bb_frequency())
			row_data.append(locus_summary.genotype_counts.get_minor_frequency())
			row_data.append(cluster_record.cluster_score.total_score)
			row_data.append(locus_summary.score_stats.gc_50)
			row_data.append(locus_summary.score_stats.gc_10)

			(hw_equilibrium, het_excess) = locus_summary.genotype_counts.compute_hardy_weinberg()
			row_data.append(het_excess)
			row_data.append(hw_equilibrium)

			row_data.append(cluster_record.cluster_score.cluster_separation)

			for cluster_stats in (cluster_record.aa_cluster_stats, cluster_record.ab_cluster_stats, cluster_record.bb_cluster_stats):
				row_data.append(cluster_stats.theta_mean)
				row_data.append(cluster_stats.theta_dev)

			for cluster_stats in (cluster_record.aa_cluster_stats, cluster_record.ab_cluster_stats, cluster_record.bb_cluster_stats):
				row_data.append(cluster_stats.r_mean)
				row_data.append(cluster_stats.r_dev)

			if len(bpm.ref_strands) > 0:
				row_data.append(RefStrand.to_string(bpm.ref_strands[i]))
			else:
				row_data.append("U")
			output_handle.write(delim.join(map(str, row_data)) + "\n")
		logger.info("Report generation complete")


def main():
	parser = argparse.ArgumentParser(
		"Generate a locus summary report from a directory of GTC files")

	parser.add_argument("gtc_directory", help="Directory containing GTC files")
	parser.add_argument("manifest_file", help="BPM manifest file")
	parser.add_argument("cluster_file", help="EGT cluster file")
	parser.add_argument("output_file", help="Location to write report")
	parser.add_argument("-p", "--project-name", dest="project_name",
						default="Project", help="A project name to report in the output header")

	args = parser.parse_args()
	logger = get_logger()
	driver(args.gtc_directory, args.manifest_file, args.cluster_file,
		   args.output_file, args.project_name, ",", logger)

if __name__ == "__main__":
	main()
