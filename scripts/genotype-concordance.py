#!/usr/bin/python

import sys
from collections import defaultdict
from collections import namedtuple
import argparse
import math
from decimal import Decimal 
import vcf
from enum import Enum

Position = namedtuple("Position", "chrom position")

class VariantType(Enum):
	snp = 0
	small_insertion = 1
	small_deletion = 2
	small_complex = 3
	midsize_insertion = 4
	midsize_deletion = 5
	midsize_complex = 6
	large_insertion = 7
	large_deletion = 8
	large_complex = 9

class Variant:
	def __init__(self, genotype, binary_genotype, variant_type : VariantType, quality=None, allele_frequency=None, unique_kmers=None):
		"""
		Represents a genotyped variant.
		"""
		self._genotype=genotype
		self._binary_genotype=binary_genotype
		self._variant_type=variant_type
		self._quality=quality
		self._allele_frequency=allele_frequency
		self._unique_kmers=unique_kmers

	def get_genotype(self):
		return self._genotype

	def get_binary_genotype(self):
		return self._binary_genotype

	def is_type(self, vartype : VariantType):
		if vartype == self._variant_type:
			return True
		else:
			return False

	def get_type(self):
		return self._variant_type

	def get_quality(self):
		return self._quality

	def get_allele_frequency(self):
		return self._allele_frequency

	def get_unique_kmers(self):
		return self._unique_kmers


def determine_variant_type(record):
	"""
	Determines the variant type.
	"""
	alleles = [record.REF] + record.ALT
	varlen = max([len(a) for a in alleles])

	if record.is_snp:
		return VariantType.snp

	is_deletion = record.var_subtype == 'del'
	is_insertion = record.var_subtype == 'ins'

	if varlen < 20:
		if is_insertion:
			return VariantType.small_insertion
		if is_deletion:
			return VariantType.small_deletion
		return VariantType.small_complex

	if varlen >= 20 and varlen <= 50:
		if is_insertion:
			return VariantType.midsize_insertion
		if is_deletion:
			return VariantType.midsize_deletion
		return VariantType.midsize_complex

	if varlen > 50:
		if is_insertion:
			return VariantType.large_insertion
		if is_deletion:
			return VariantType.large_deletion
		return VariantType.large_complex


def extract_call(record, read_gl=False, read_qual=False):
	"""
	Extract genotype information from a VCF Record.
	"""

	# genotype alleles
	genotype_sequences = None
	# binary genotype
	binary_genotype = 3
	# genotype quality
	likelihood = 0
	# allele frequency of least frequent allele
	allele_frequency = 0.0
	unique_kmers = 0
	# determine the type of variant
	variant_type = determine_variant_type(record)

	if (record.samples[0]['GT'][0] != '.') and (record.samples[0]['GT'][-1] != '.'):

		# extract genotype information
		genotype_string = record.samples[0]['GT'].replace('/', '|')
		genotype_list = genotype_string.split('|')
		genotype = (int(genotype_list[0]), int(genotype_list[1]))
		variant_alleles = [record.REF] + record.ALT
		genotype_sequences = set([str(variant_alleles[genotype[0]]), str(variant_alleles[genotype[1]])])

		# compute likelihood of genotype
		if read_gl:
			# check if GQ field is present
			if 'GQ' in record.FORMAT.split(':'):
				likelihood = int(record.samples[0]['GQ'])
		if read_qual:
			if record.QUAL is not None:
				likelihood = int(record.QUAL)

		# in case of bi-allelic SNP store binary gt
		if len(variant_alleles) < 3:
			binary_genotype = genotype[0] + genotype[1]
		else:
			binary_genotype = -1

		# determine the allele frequency (of the least frequent genotype allele)
		if 'AF' in record.INFO:
			info_fields = record.INFO['AF']
			# add frequency of reference allele
			frequencies = [1.0 - sum(info_fields)] + info_fields
			allele_frequency = min(frequencies)

		# determine minimum number of kmers that cover each allele
		if 'AK' in record.INFO:
			info_fields = record.INFO['AK']
			# determine minimum of unique kmers (exclude alleles not covered by paths)
			unique_kmers = min([i for i in info_fields if i >= 0])

	return Position(record.CHROM, record.POS), Variant(genotype_sequences, binary_genotype, variant_type, likelihood, allele_frequency, unique_kmers)


class GenotypingStatistics:

	def __init__(self, quality, allele_frequency, unique_kmers):

		# thresholds used to filter variant set
		self.quality = quality
		self.allele_frequency = allele_frequency
		self.unique_kmers = unique_kmers

		# numbers of variants in baseline and callset
		self.total_baseline = 0
		self.total_baseline_biallelic = 0
		self.total_intersection = 0

		# fractions of correct, wrong and untyped variants
		self.correct_all = 0
		self.wrong_all = 0
		self.not_typed_all = 0

		self.not_in_callset_all = 0
		self.not_in_callset_biallelic = 0

		# confusion matrix
		self.confusion_matrix = [[0 for x in range(4)] for y in range(4)]

	def print_stats_to_file(self, tsv_file):
		typed_all = max(self.correct_all + self.wrong_all, 1)
		assert self.total_baseline == self.correct_all + self.wrong_all + self.not_in_callset_all + self.not_typed_all
		correct_biallelic = self.confusion_matrix[0][0] + self.confusion_matrix[1][1] + self.confusion_matrix[2][2]
		wrong_biallelic = self.confusion_matrix[0][1] + self.confusion_matrix[0][2] + self.confusion_matrix[1][0] + self.confusion_matrix[1][2] + self.confusion_matrix[2][1] + self.confusion_matrix[2][0]
		not_typed_biallelic = self.total_baseline_biallelic - correct_biallelic - wrong_biallelic - self.not_in_callset_biallelic
		typed_biallelic = max(correct_biallelic + wrong_biallelic, 1)

		tsv_file.write('\t'.join([str(self.quality), # quality
																	str(self.allele_frequency), # allele freq
																	str(self.unique_kmers), # unique kmer count
																	str(self.total_baseline), # total baseline variants
																	str(self.total_baseline_biallelic), # total biallelic variants
																	str(self.total_intersection), # intersection of callsets
																	str(self.correct_all/float(typed_all)), # correct
																	str(self.wrong_all/float(typed_all)), # wrong
																	str((self.not_typed_all+self.not_in_callset_all)/float(self.total_baseline if self.total_baseline is not 0 else 1)), # not typed
																	str(correct_biallelic/ float(typed_biallelic)), # correct biallelic
																	str(wrong_biallelic/ float(typed_biallelic)), # wrong biallelic
																	str((not_typed_biallelic+self.not_in_callset_biallelic)/float(self.total_baseline_biallelic if self.total_baseline_biallelic is not 0 else 1)),  # not typed biallelic
																	str(self.correct_all),
																	str(self.wrong_all),
																	str(self.not_typed_all),
																	str(self.not_in_callset_all),
																	str(correct_biallelic),
																	str(wrong_biallelic),
																	str(not_typed_biallelic),
																	str(self.not_in_callset_biallelic)
														]) + '\n'
											)

	def print_matrix_to_file(self, txt_file):
		txt_file.write('\n###############################################################################################################################################\n')
		txt_file.write('#       Matrix for quality: '  + str(self.quality) +', allele frequency threshold: ' + str(self.allele_frequency) + ' and unique kmer count threshold: ' + str(self.unique_kmers) + '#\n')
		txt_file.write('  ###############################################################################################################################################\n\n')
		matrix_string = "\t0\t1\t2\t./.\n"
		for true in range(4):
			line = str(true) + '\t'
			for geno in range(4):
				line += str(self.confusion_matrix[true][geno]) + '\t'
			matrix_string += line + '\n'
		txt_file.write(matrix_string)


class GenotypeConcordanceComputer:
	def __init__(self, baseline_vcf, callset_vcf, use_qual):
		self._baseline_variants = {}
		self._callset_variants = defaultdict(list)
		self._duplicated_positions = defaultdict(lambda:False)
		self._total_baseline = 0
		self._total_callset = 0

		# read baseline variants
		for record in vcf.Reader(open(baseline_vcf, 'r')):
			# determine variant type
			pos, variant = extract_call(record)
			if pos in self._baseline_variants:
				# duplicated position, ignore it
				sys.stderr.write('Warning: position ' + str(pos.chrom) + ' ' + str(pos.position) + ' occurs more than once and will be skipped.\n')
				self._duplicated_positions[pos] = True
				continue
			self._total_baseline += 1
			self._baseline_variants[pos] = variant

		# read callset variants
		for record in vcf.Reader(open(callset_vcf, 'r')):
			if use_qual:
				pos, variant = extract_call(record, read_qual = True)
			else:
				pos, variant = extract_call(record, read_gl = True)
			self._callset_variants[pos].append(variant)
			self._total_callset += 1

		sys.stderr.write('Read ' + str(self._total_baseline) + ' variants from the baseline VCF (including duplicates).\n')
		sys.stderr.write('Read ' + str(self._total_callset) + ' variants from the callset VCF (including duplicates).\n')

	def nr_baseline_variants(self):
		return self._total_baseline

	def nr_callset_variants(self):
		return self._total_callset

	def print_statistics(self, varianttypes, quality, allele_freq, uk_count, txt_files, tsv_files):
		"""
		Compute genotype concordance statistics for all variants with
		a quality at least quality and for which the allele frequency 
		of the least frequent allele is at least allele_freq.
		"""

		# make sure there is a txt and tsv file for each variant type
		assert len(txt_files) == len(tsv_files) == len(varianttypes)

		# keep variant statistics for all variant types
		statistics = { var_type : GenotypingStatistics(quality, allele_freq, uk_count) for var_type in varianttypes }

		# check genotype predictions
		for pos, gt in self._baseline_variants.items():

			# if position was present multiple times in baseline, skip it
			if self._duplicated_positions[pos]:
				continue

			vartype = gt.get_type()
			statistics[vartype].total_baseline += 1

			callset_predictions = self._callset_variants[pos]

			if gt.get_binary_genotype() != -1:
				statistics[vartype].total_baseline_biallelic += 1

			if len(callset_predictions) > 1:
				print('Warning: multiple predictions for variant at position ' + str(pos.chrom) + ' ' + str(pos.position) + '. Using first.')

			if len(callset_predictions) > 0:
				callset_gt = callset_predictions[0]
				statistics[vartype].total_intersection += 1

				if callset_gt.get_quality() < quality or callset_gt.get_allele_frequency() < allele_freq or callset_gt.get_unique_kmers() < uk_count:
					statistics[vartype].not_typed_all += 1
					continue

				if (callset_gt.get_binary_genotype() != -1) and (gt.get_binary_genotype() != -1):
					statistics[vartype].confusion_matrix[gt.get_binary_genotype()][callset_gt.get_binary_genotype()] += 1

				if (callset_gt.get_genotype() == None) or (gt.get_genotype() == None):
					statistics[vartype].not_typed_all += 1
					continue

				# check if genotype predictions are identical
				if callset_gt.get_genotype() == gt.get_genotype():
					statistics[vartype].correct_all += 1
				else:
					print('wrong', vartype.name, pos.chrom, pos.position, callset_gt.get_genotype(), gt.get_genotype(), gt.get_binary_genotype())
					statistics[vartype].wrong_all += 1
			else:
				statistics[vartype].not_in_callset_all += 1
				if gt.get_binary_genotype() != -1:
					statistics[vartype].confusion_matrix[gt.get_binary_genotype()][3] += 1
					statistics[vartype].not_in_callset_biallelic += 1

		assert (statistics[vartype].correct_all + statistics[vartype].wrong_all + statistics[vartype].not_typed_all == statistics[vartype].total_intersection)
		# write statistics and confusion matrices to file
		for i,vartype in enumerate(varianttypes):
			statistics[vartype].print_stats_to_file(tsv_files[i])
			statistics[vartype].print_matrix_to_file(txt_files[i])

if __name__ == "__main__":

	# baseline: the file containing variants and true genotypes
	# callset: contains genotype predictions of genotyper	
	parser = argparse.ArgumentParser(prog='genotype-concordance.py', description=__doc__)
	parser.add_argument('baseline', metavar='BASELINE', help='baseline VCF (ground truth).')
	parser.add_argument('callset', metavar='CALLSET', help='callset VCF (genotyped variants).')
	parser.add_argument('prefix', metavar='OUTFILE', help='prefix of the output file name.')
	parser.add_argument('--qualities', default='0', metavar='GQ-THRESHOLDS', help='comma separated list of GQ-thresholds to consider (default: consider all variants regardless of quality).')
	parser.add_argument('--allele-frequencies', default='0', metavar='AF-THRESHOLDS', help='comma separated list of allele frequency thresholds to consider (default: consider all variants regardless of allele frequency.).')
	parser.add_argument('--unique-kmers', default='0', metavar='UK-THRESHOLDS', help='comma separated list of unique kmer counts for a variant (only works for VCFs with UK tag.).')
	parser.add_argument('--use-qual', default=False, action='store_true', help='use qualities in QUAL field instead of GQ fields (default: use GQ field).')
	args = parser.parse_args()

	quality_thresholds = [int(i) for i in args.qualities.split(',')]
	allele_freq_thresholds = [float(i) for i in args.allele_frequencies.split(',')]
	uk_thresholds = [int(i) for i in args.unique_kmers.split(',')]
	genotype_concordance = GenotypeConcordanceComputer(args.baseline, args.callset, args.use_qual)

	# prepare output files
	variantnames = []
	for vartype in VariantType:
		name = vartype.name.replace('_', '-')
		variantnames.append(name)

	tsv_files = [open(args.prefix + '_' + name + '.tsv', 'w') for name in variantnames]
	txt_files = [open(args.prefix + '_' + name + '.txt', 'w') for name in variantnames]

	# write headers to tsv files
	header = '\t'.join(['quality',
											'allele_frequency',
											'unique_kmers',
											'total_baseline',
											'total_baseline_biallelic',
											'total_intersection',
											'correct_all',
											'wrong_all',
											'not_typed_all',
											'correct_biallelic',
											'wrong_biallelic',
											'not_typed_biallelic',
											'nr_correct_all',
											'nr_wrong_all',
											'nr_not_typed_all',
											'nr_not_in_callset_all',
											'nr_correct_biallelic',
											'nr_wrong_biallelic',
											'nr_not_typed_biallelic',
											'nr_not_in_callset_biallelic'
								]) + '\n'
	for tsv_file in tsv_files:
		tsv_file.write(header)

	# compute statistics
	vartypes = [vartype for vartype in VariantType]
	genotype_concordance.print_statistics(vartypes, 0, 0.0, 0, txt_files, tsv_files)
	for uk_count in uk_thresholds:
		for quality in quality_thresholds:
			if uk_count == quality == 0:
				continue
			genotype_concordance.print_statistics(vartypes, quality, 0.0, uk_count, txt_files, tsv_files)

	for allele_freq in allele_freq_thresholds:
		if allele_freq == 0:
			continue
		genotype_concordance.print_statistics(vartypes, 0, allele_freq, 0, txt_files, tsv_files)

	# close files
	for i in range(len(tsv_files)):
		tsv_files[i].close()
		txt_files[i].close()