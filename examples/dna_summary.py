import sys
import os
import argparse
import logging
from math import isnan
from numpy import percentile
from scipy.stats import norm

from IlluminaBeadArrayFiles import BeadPoolManifest, GenotypeCalls, RefStrand

nan = float("nan")


def compute_genotypes(genotypes):
    """
    Summarize information about genotype counts for diploid genotyping counting

    Returns:
        list: genotype counts and frequencies
    """
    # Count genotypes
    no_calls = genotypes.count(0)
    aa_count = genotypes.count(1)
    ab_count = genotypes.count(2)
    bb_count = genotypes.count(3)

    # Get the number of calls (i.e., not no-calls)
    num_calls = aa_count + ab_count + bb_count

    # Get the call rate
    if (num_calls + no_calls) > 0: 
        call_rate = num_calls / float(num_calls + no_calls)
    else:
        call_rate = nan
    
    # Genotype frequencies (as fraction of all calls)
    if num_calls > 0:
        # Frequency of AA genotype
        aa_freq = round(aa_count / float(num_calls),4)
        # Frequency of AB genotype
        ab_freq = round(ab_count / float(num_calls),4)
        # Frequency of BB genotype
        bb_freq = round(bb_count / float(num_calls),4)
    else:
        aa_freq = nan
        ab_freq = nan
        bb_freq = nan
    
    # Computes and return the minor allele frequency. If no calls, will be NaN
    a_allele_count = aa_count * 2 + ab_count
    if num_calls > 0:
        a_frequency = a_allele_count / float(2 * num_calls)
        minor_freq = round(min(a_frequency, 1.0 - a_frequency),4)
    else:
        minor_freq = nan
    return [no_calls, num_calls, call_rate, aa_freq, ab_freq, bb_freq, minor_freq]



def ScoreStatistics(scores, percentile):
    """
    Capture statistics related to the gencall score distribution

    Args:
        scores (list(float)): A list of gencall scores
        percentile (int): percentile to calculate
            gc_10 : 10th percentile of Gencall score distribution
            gc_50 : 50th percentile of Gencall score distribution
        
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
    score_stat = y1 + (y2 - y1) / (x2 - x1) * (percentile - x1)
    score_stat = round(score_stat,4)
    return score_stat




def get_logger():
    # set up log file
    # create logger
    logger = logging.getLogger('DNA Report')
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


def driver(gtc_dir, manifest_filename, output_filename, project_name, delim, logger):
    logger.info("Reading manifest file")
    bpm = BeadPoolManifest(manifest_filename)

    samples = []
    logger.info("Initializing genotype data")
    gtc_files = []
    for gtc_file in os.listdir(gtc_dir):
        if gtc_file.endswith(".gtc"):
            gtc_files.append(os.path.join(gtc_dir, gtc_file))




    logger.info("Generating report")
    loci = range(len(bpm.normalization_lookups))
    with open(output_filename, "w") as output_handle:
        output_handle.write("DNA Report on " + os.path.abspath(output_filename) + "\n")
        header = [""]
        header.append("# LOCI = {}".format(len(loci)))
        header.append("# DNAs = {}".format(len(samples)))
        header.append("ProjectName = {}".format(project_name))
        header.append("GenCall Version = NaN")
        header.append("Low GenCall Score Cutoff = NaN")

        output_handle.write(delim.join(header) + "\n")
        output_handle.write(delim.join("Row,DNA_ID,#No_Calls,#Calls,Call_Rate,A/A_Freq,A/B_Freq,B/B_Freq,Minor_Freq,50%_GC_Score,10%_GC_Score".split(",")) + "\n")
        row = 0
        for gtc_file in gtc_files:
            row += 1
            gtc = GenotypeCalls(gtc_file)
            genotypes = gtc.get_genotypes()
            scores = gtc.get_genotype_scores()
            assert len(genotypes) == len(bpm.names)
            row_data = []
            row_data.append(row)
            row_data.append(gtc.get_sample_name())
            row_data += compute_genotypes(genotypes)
            row_data.append(ScoreStatistics(scores, 50))
            row_data.append(ScoreStatistics(scores, 10))
            output_handle.write(delim.join(map(str, row_data)) + "\n")
        logger.info("Report generation complete")


def main():
    parser = argparse.ArgumentParser(
        "Generate a DNA report from a directory of GTC files")

    parser.add_argument("gtc_directory", help="Directory containing GTC files")
    parser.add_argument("manifest_file", help="BPM manifest file")
    parser.add_argument("output_file", help="Location to write report")
    parser.add_argument("-p", "--project-name", dest="project_name",
                        default="Project", help="A project name to report in the output header")

    args = parser.parse_args()
    logger = get_logger()
    driver(args.gtc_directory, args.manifest_file, args.output_file, args.project_name, ",", logger)

if __name__ == "__main__":
    main()
