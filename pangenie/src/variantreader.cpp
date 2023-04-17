#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <regex>
#include <algorithm>
#include "variantreader.hpp"


using namespace std;

std::map< std::string, std::vector<Variant> > VariantReader::variants_per_chromosome;
std::map< std::string, std::vector<std::vector<std::string>>> VariantReader::variant_ids;
//FastaReader VariantReader::fasta_reader;

void parse_line(vector<DnaSequence>& result, string line, char sep) {
	string token;
	istringstream iss (line);
	while (getline(iss, token, sep)) {
		result.push_back(DnaSequence(token));
	}
}

void parse_line(vector<string>& result, string line, char sep) {
	string token;
	istringstream iss (line);
	while (getline(iss, token, sep)) {
		result.push_back(token);
	}
}

void VariantReader::insert_ids(string& chromosome, vector<DnaSequence>& alleles, vector<string>& variant_ids, bool reference_added) {
	vector<unsigned char> index = construct_index(alleles, reference_added);
	assert(index.size() < 256);
	// insert IDs in the lex. order of their corresponding alleles
	vector<string> sorted_ids;
	for (auto id : index) {
		sorted_ids.push_back(variant_ids[id]);
	}
	this->variant_ids[chromosome].push_back(sorted_ids);
}

string VariantReader::get_ids(string chromosome, vector<string>& alleles, size_t variant_index, bool reference_added) {
	vector<unsigned char> index = construct_index(alleles, reference_added);
	assert(index.size() < 256);
	vector<string> sorted_ids(index.size());
	for (unsigned char i = 0; i < index.size(); ++i) {
		sorted_ids[index[i]] = this->variant_ids.at(chromosome).at(variant_index)[i];
	}

	string result = "";
	for (unsigned char i = 0; i < sorted_ids.size(); ++i) {
		if (i > 0) result += ',';
		result += sorted_ids[i];
	}
	return result;
}

bool matches_pattern(string& sequence, regex& r) {
	size_t total_length = sequence.size();
	string control = "";
	auto begin = sequence.begin();
	size_t length = 1000;
	size_t letters_seen = 0;
	while (begin + letters_seen != sequence.end()) {
		size_t next_interval = min(length, total_length - letters_seen);
		if (!regex_match(begin + letters_seen, begin + letters_seen + next_interval, r)) {
			return false;
		}
		control += string(begin + letters_seen, begin + letters_seen + next_interval);
		letters_seen += next_interval;
	}
	assert(control == sequence);
	return true;
}

void parse_info_fields(vector<string>& result, string line) {
	vector<string> fields;
	parse_line(fields, line, ';');
	for (auto s : fields) {
		if (s.rfind("ID=", 0) == 0) {
			// ID field present
			parse_line(result, s.substr (3), ',');
		}
	}
}
VariantReader::VariantReader (const VariantReader &old_obj)
{
  kmer_size               =old_obj.kmer_size;
  nr_paths                =old_obj.nr_paths;
  nr_variants             =old_obj.nr_variants;
  add_reference           =old_obj.add_reference;
  samples                  =old_obj.samples;
  genotyping_outfile_open =false;
  phasing_outfile_open    =false;
}

VariantReader& VariantReader::operator=(const VariantReader& other)
{
    //cout<<"begin copy"<<endl;
    fasta_reader            =other.fasta_reader;
    //cout<<"Here"<<endl;
    kmer_size               =other.kmer_size;
    nr_paths                =other.nr_paths;
    nr_variants             =other.nr_variants;
    add_reference           =other.add_reference;
    samples                 =other.samples;
    genotyping_outfile_open =false;
    phasing_outfile_open    =false;
    return *this;
}

VariantReader::VariantReader(string filename, string reference_filename, size_t kmer_size, bool add_reference, vector<string> sample)
	:
     fasta_reader(reference_filename),
     kmer_size(kmer_size),
	 nr_variants(0),
	 add_reference(add_reference),
	 samples(sample),
	 genotyping_outfile_open(false),
	 phasing_outfile_open(false)
{
	if (filename.substr(filename.size()-3,3).compare(".gz") == 0) {
		throw runtime_error("VariantReader::VariantReader: Uncompressed VCF-file is required.");
	}
	ifstream file(filename);
	if (!file.good()) {
		throw runtime_error("VariantReader::VariantReader: input VCF file cannot be opened.");
	}
	string line;
	string previous_chrom("");
	size_t previous_end_pos = 0;
	map<unsigned int, string> fields = { {0, "#CHROM"}, {1, "POS"}, {2, "ID"}, {3, "REF"}, {4, "ALT"}, {5, "QUAL"}, {6, "FILTER"}, {7, "INFO"}, {8, "FORMAT"} };
	vector<Variant> variant_cluster;
	// read VCF-file line by line
	while (getline(file, line)) {
		if (line.size() == 0) continue;
		vector<string> tokens;
		parse_line(tokens, line, '\t');
		// header lines?
		if (tokens[0].substr(0,2) == "##") continue;
		if (tokens[0].at(0) == '#') {
			// check number of samples/paths given
			if (tokens.size() < 9) {
				throw runtime_error("VariantReader::VariantReader: not a proper VCF-file.");
			}
			if (tokens.size() < 10) {
				throw runtime_error("VariantReader::VariantReader: no haplotype paths given.");
			}
			// validate header line
			for (unsigned int i = 0; i < 9; ++i) {
				if (tokens[i] != fields[i]) {
					throw runtime_error("VariantReader::VariantReader: VCF header line is malformed.");
				}
			}
			this->nr_paths = (tokens.size() - 9)*2;
			// add one for reference path
			if (add_reference) this->nr_paths += 1;
			continue;
		}
		if (tokens.size() < 10) {
			throw runtime_error("VariantReader::VariantReader: malformed VCF-file, or no haplotype paths given in VCF.");
		}
		// get chromosome
		string current_chrom = tokens[0];
        string variantID=tokens[2];
        if(chromsomes.empty() || chromsomes.back()!= current_chrom)
            chromsomes.push_back(current_chrom);
		// get position
		size_t current_start_pos;
		stringstream sstream(tokens[1]);
		sstream >> current_start_pos;
		// VCF positions are 1-based
		current_start_pos -= 1;
		// if variant is contained in previous one, skip it
		if ((previous_chrom == current_chrom) && (current_start_pos < previous_end_pos)) {
			cerr << "VariantReader: skip variant at " << current_chrom << ":" << current_start_pos << " since it is contained in a previous one."  << endl;
			continue;
		}
		// if distance to next variant is larger than kmer_size, start a new cluster
		if ( (previous_chrom != current_chrom) || (current_start_pos - previous_end_pos) >= (kmer_size-1) ) {
			// merge all variants currently in cluster and store them
			add_variant_cluster(previous_chrom, &variant_cluster);
			variant_cluster.clear();
		}
		// get REF allele
		DnaSequence ref(tokens[3]);
		DnaSequence observed_allele;
		this->fasta_reader.get_subsequence(current_chrom, current_start_pos, current_start_pos + ref.size(), observed_allele);
		if (ref != observed_allele) {
			throw runtime_error("VariantReader::VariantReader: reference allele given in VCF does not match allele in reference fasta file at that position.");
		}
		size_t current_end_pos = current_start_pos + ref.size();
		// get ALT alleles
		vector<DnaSequence> alleles = {ref};
		// make sure alt alleles are given explicitly
		regex r("^[CAGTcagt,]+$");
		if (!matches_pattern(tokens[4], r)) {
//		if (!regex_match(tokens[4], r)) {
			// skip this position
			cerr << "VariantReader: skip variant at " << current_chrom << ":" << current_start_pos << " since alleles contain undefined nucleotides: " << tokens[4] << endl;
			continue;
		}
		parse_line(alleles, tokens[4], ',');

		// currently, number of alleles is limited to 256
		if (alleles.size() > 255) {
			throw runtime_error("VariantReader: number of alternative alleles is limited to 254 in current implementation. Make sure the VCF contains only alternative alleles covered by at least one of the haplotypes.");
		}

		// TODO: handle cases where variant is less than kmersize from start or end of the chromosome
		if ( (current_start_pos < (kmer_size*2) ) || ( (current_end_pos + (kmer_size*2)) > this->fasta_reader.get_size_of(current_chrom)) ) {
			cerr << "VariantReader: skip variant at " << current_chrom << ":" << current_start_pos << " since variant is less than 2 * kmer size from start or end of chromosome. " << endl;
			continue;
		}

		// store mapping of alleles to variant ids
		vector<string> var_ids;
		parse_info_fields(var_ids, tokens[7]);
		if (!var_ids.empty()) {
			insert_ids(current_chrom, alleles, var_ids, true);
		} else {
			this->variant_ids[current_chrom].push_back(vector<string>());
		}

		// make sure that there are at most 255 paths (including reference path in case it is requested)
		if (this->nr_paths > 255) {
			throw runtime_error("VariantReader: number of paths is limited to 254 in current implementation.");
		}


		// construct paths
		vector<unsigned char> paths = {};
		if (add_reference) paths.push_back((unsigned char) 0);
		unsigned char undefined_index = alleles.size();
		string undefined_allele = "N";
        bool phased=true;
		for (size_t i = 9; i < tokens.size(); ++i) {
			// make sure all genotypes are phased
            char gtChar='|';

			if (tokens[i].find('/') != string::npos) {
                gtChar='/';
                phased=false;
				//throw runtime_error("VariantReader::VariantReader: Found unphased genotype.");
			}
			vector<string> p ;
			parse_line(p, tokens[i], gtChar);
			if (p.size() != 2) {
				throw runtime_error("VariantReader::VariantReader: Found invalid genotype. Genotypes must be diploid (.|. if missing).");
			}
			for (string& s : p){
				// handle unknown genotypes '.'
				if (s == ".") {
					// add "NNN" allele to the list of alleles
					parse_line(alleles, undefined_allele, ',');
					paths.push_back(undefined_index);
					assert(undefined_index < 255);
					undefined_index += 1;
					undefined_allele += "N";
				} else {
					unsigned int p_index = atoi(s.c_str());
					if (p_index >= alleles.size()) {
						throw runtime_error("VariantReader::VariantReader: invalid genotype in VCF.");
					}
					assert(p_index < 255);
					paths.push_back( (unsigned char) p_index);
				}
			}
		}


		// determine left and right flanks
		DnaSequence left_flank;
		this->fasta_reader.get_subsequence(current_chrom, current_start_pos - kmer_size + 1, current_start_pos, left_flank);
		DnaSequence right_flank;
		this->fasta_reader.get_subsequence(current_chrom, current_end_pos, current_end_pos + kmer_size - 1, right_flank);
		// add Variant to variant_cluster
		Variant variant (left_flank, right_flank, current_chrom, current_start_pos, current_end_pos, alleles, paths,variantID,phased);
		variant_cluster.push_back(variant);
		previous_chrom = current_chrom;
		previous_end_pos = current_end_pos;
	}
	// add last cluster to list
	add_variant_cluster(previous_chrom, &variant_cluster);


    for(auto chrom: this->chromsomes)
    {
        this->variantsStatsPerSample[chrom]=std::vector<std::vector<unsigned short>>(size_of(chrom));
    }

	cerr << "Identified " << this->nr_variants << " variants in total from VCF-file." << endl;
}

void VariantReader::addVariantStat(unsigned int variantID, unsigned int sampleID,std::string chromosome,
                                   vector<unsigned char>& defined_alleles,std::vector<VariantStats> & singleton_stats)
{
    size_t stat_size=(defined_alleles.size()+2)* singleton_stats.size();
    if(variantsStatsPerSample[chromosome][variantID].size() == 0)
    {
        size_t size= stat_size* this->samples.size();
        variantsStatsPerSample[chromosome][variantID].resize(size);
    }
// UK
    for (size_t j = 0; j < singleton_stats.size(); ++j){
        size_t index= stat_size* sampleID;
        variantsStatsPerSample[chromosome][variantID][index++]=singleton_stats.at(j).nr_unique_kmers;
        variantsStatsPerSample[chromosome][variantID][index++]=singleton_stats[j].coverage;
        for (unsigned int a = 0; a < defined_alleles.size(); ++a) {
            variantsStatsPerSample[chromosome][variantID][index++]=singleton_stats.at(j).kmer_counts[a];
        }
    }

}

size_t VariantReader::get_kmer_size() const {
	return this->kmer_size;
}

void VariantReader::write_path_segments(std::string filename) const {
	ofstream outfile;
	outfile.open(filename);
	if (!outfile.good()) {
		stringstream ss;
		ss << "VariantReader::write_path_segments: File " << filename << " cannot be created. Note that the filename must not contain non-existing directories." << endl;
		throw runtime_error(ss.str());
	}
	// make sure to capture all chromosomes in the reference (including such for which no variants are given)
	vector<string> chromosome_names;
	this->fasta_reader.get_sequence_names(chromosome_names);
	for (auto element : chromosome_names) {
		size_t prev_end = 0;
		// check if chromosome was present in VCF and write allele sequences in this case
		auto it = this->variants_per_chromosome.find(element);
		if (it != this->variants_per_chromosome.end()) {
			for (auto variant : this->variants_per_chromosome.at(element)) {
				// generate reference unitig and write to file
				size_t start_pos = variant.get_start_position();
				outfile << ">" << element << "_reference_" << start_pos << "\n";
				string ref_segment;
				this->fasta_reader.get_subsequence(element, prev_end, start_pos, ref_segment);
				outfile << ref_segment << "\n";
				for (size_t allele = 0; allele < variant.nr_of_alleles(); ++allele) {
					// sequence name
					outfile << ">" << element << "_" << start_pos << "_" << allele << "\n";
					outfile << variant.get_allele_string(allele) << "\n";
				}
				prev_end = variant.get_end_position();
			}
		}
		// output reference sequence after last position on chromosome
		outfile << ">" << element << "_reference_end" << "\n";
		size_t chr_len = this->fasta_reader.get_size_of(element);
		string ref_segment;
		this->fasta_reader.get_subsequence(element, prev_end, chr_len, ref_segment);
		outfile << ref_segment << "\n";
	}
	outfile.close();
}

void VariantReader::get_chromosomes(vector<string>* result) const {

	// sort the chromosomes by size (decending).
	vector<pair<size_t,string>> chromosome_sizes;
	for (auto const& element : this->variants_per_chromosome) {
		chromosome_sizes.push_back(make_pair(element.second.size(), element.first));
	}
	sort(chromosome_sizes.rbegin(), chromosome_sizes.rend());

	// return chromosomes in sorted order
	for (auto const& element : chromosome_sizes) {
		result->push_back(element.second);
	}
}
std::vector<std::string> VariantReader::get_chromosomes_vcfSorted() const
{
    return chromsomes;
}

size_t VariantReader::size_of(string chromosome) const {
	if (this->variants_per_chromosome.find(chromosome) != this->variants_per_chromosome.end()) {
		return this->variants_per_chromosome.at(chromosome).size();
	} else {
		return 0;
	}
}

const Variant& VariantReader::get_variant(string chromosome, size_t index) const {
	if (index < size_of(chromosome)) {
		return this->variants_per_chromosome.at(chromosome)[index];
	} else {
		throw runtime_error("VariantReader::get_variant: index out of bounds.");
	}
}

void VariantReader::add_variant_cluster(string& chromosome, vector<Variant>* cluster) {
	if (!cluster->empty()) {
		// merge all variants in cluster
		Variant combined = cluster->at(0);
		for (size_t v = 1; v < cluster->size(); ++v) {
			combined.combine_variants(cluster->at(v));
		}
		combined.add_flanking_sequence();
		this->variants_per_chromosome[chromosome].push_back(combined);
		this->nr_variants += 1;
	}
}

const vector<Variant>& VariantReader::get_variants_on_chromosome(string chromosome) const {
	if (this->variants_per_chromosome.find(chromosome) != this->variants_per_chromosome.end()) {
		return this->variants_per_chromosome.at(chromosome);
	} else {
		throw runtime_error("VariantReader::get_variants_on_chromosome: chromosome " + chromosome + " not present in VCF.");
	}
}

string get_date() {
	time_t t = time(0);
	tm* now = localtime(&t);
	ostringstream oss;
	oss << (now->tm_year+1900) << setfill('0') << setw(2) << (now->tm_mon+1) << setw(2) << now->tm_mday;
	return string(oss.str());
}

void VariantReader::open_genotyping_outfile(string filename) {
    this->genotyping_outfile=std::vector<std::ofstream>(samples.size());
    for(unsigned i=0;i< samples.size();i++) {
        this->genotyping_outfile[i].open(filename+"_"+samples[i]+"_genotyping.vcf");
        if (!this->genotyping_outfile[i].is_open()) {
            throw runtime_error(
                    "VariantReader::open_genotyping_outfile: genotyping output file cannot be opened. Note that the filename must not contain non-existing directories.");
        }

        this->genotyping_outfile_open = true;

        // write VCF header lines
        this->genotyping_outfile[i] << "##fileformat=VCFv4.2" << endl;
        this->genotyping_outfile[i] << "##fileDate=" << get_date() << endl;
        // TODO output command line
        this->genotyping_outfile[i] << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" << endl;
        this->genotyping_outfile[i] << "##INFO=<ID=UK,Number=1,Type=Integer,Description=\"Total number of unique kmers.\">"
                                 << endl;
        this->genotyping_outfile[i]
                << "##INFO=<ID=AK,Number=R,Type=Integer,Description=\"Number of unique kmers per allele. Will be -1 for alleles not covered by any input haplotype path\">"
                << endl;
        this->genotyping_outfile[i]
                << "##INFO=<ID=MA,Number=1,Type=Integer,Description=\"Number of alleles missing in panel haplotypes.\">"
                << endl;
        this->genotyping_outfile[i] << "##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs.\">" << endl;
        this->genotyping_outfile[i] << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
        this->genotyping_outfile[i]
                << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality: phred scaled probability that the genotype is wrong.\">"
                << endl;
        this->genotyping_outfile[i]
                << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Comma-separated log10-scaled genotype likelihoods for absent, heterozygous, homozygous.\">"
                << endl;
        this->genotyping_outfile[i] << "##FORMAT=<ID=KC,Number=1,Type=Float,Description=\"Local kmer coverage.\">" << endl;
        this->genotyping_outfile[i] << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        this->genotyping_outfile[i] << "\t" << samples[i];
        this->genotyping_outfile[i] << endl;
    }
}

void VariantReader::open_phasing_outfile(string filename) { 
	this->phasing_outfile.open(filename);
	if (! this->phasing_outfile.is_open()) {
		throw runtime_error("VariantReader::open_phasing_outfile: phasing output file cannot be opened.");
	}

	this->phasing_outfile_open = true;

	// write VCF header lines
	this->phasing_outfile << "##fileformat=VCFv4.2" << endl;
	this->phasing_outfile << "##fileDate=" << get_date() << endl;
	// TODO output command line
	this->phasing_outfile << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" << endl;
	this->phasing_outfile << "##INFO=<ID=UK,Number=1,Type=Integer,Description=\"Total number of unique kmers.\">" << endl;
	this->phasing_outfile << "##INFO=<ID=AK,Number=R,Type=Integer,Description=\"Number of unique kmers per allele. Will be -1 for alleles not covered by any input haplotype path.\">" << endl;
	this->phasing_outfile << "##INFO=<ID=MA,Number=1,Type=Integer,Description=\"Number of alleles missing in panel haplotypes.\">" << endl;
	this->phasing_outfile << "##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs.\">" << endl;
	this->phasing_outfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	this->phasing_outfile << "##FORMAT=<ID=KC,Number=1,Type=Float,Description=\"Local kmer coverage.\">" << endl;
	this->phasing_outfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for(auto s :samples)
          this->phasing_outfile <<"\t"<<s;
        this->phasing_outfile << "\n";

}

void VariantReader::write_genotypes_of(string chromosome, string sample_name,const vector<GenotypingResult>& genotyping_result, bool ignore_imputed) {
	// outfile needs to be open
	if (!this->genotyping_outfile_open) {
		throw runtime_error("VariantReader::write_genotypes_of: output file needs to be opened before writing.");
	}

	if (this->variants_per_chromosome.find(chromosome) == this->variants_per_chromosome.end()) {
		cerr << "VariantReader::write_genotypes_of: no variants for given chromosome were written." << endl;
		return;
	}

	if (genotyping_result.size() != size_of(chromosome)) {
		throw runtime_error("VariantReader::write_genotypes_of: number of variants and number of computed genotypes differ.");
	}
    auto sample_it= std::find(samples.begin(),samples.end(), sample_name);
    if(sample_it == samples.end())
    {
        throw runtime_error("VariantReader::write_genotypes_of: unknown sample.");

    }
    unsigned  sample_index= sample_it-samples.begin();
	size_t counter = 0;
	for (size_t i = 0; i < size_of(chromosome); ++i) {
		Variant variant = this->variants_per_chromosome.at(chromosome)[i];

		// separate (possibly combined) variant into single variants and print a line for each
		vector<Variant> singleton_variants;
		vector<GenotypingResult> singleton_likelihoods;
//		vector<VariantStats> singleton_stats;
		variant.separate_variants(&singleton_variants, &genotyping_result.at(i), &singleton_likelihoods);
//		variant.variant_statistics(unique_kmers->at(i), singleton_stats);
      //  vector<VariantStats> singleton_stats= variantsStatsPerSample[chromosome][sample_name][i];
		for (size_t j = 0; j < singleton_variants.size(); ++j) {
			Variant v = singleton_variants[j];
			v.remove_flanking_sequence();
			this->genotyping_outfile[sample_index] << v.get_chromosome() << "\t"; // CHROM
			this->genotyping_outfile[sample_index] << (v.get_start_position() + 1) << "\t"; // POS
			this->genotyping_outfile[sample_index] << v.get_id() << "\t"; // ID
			this->genotyping_outfile[sample_index] << v.get_allele_string(0) << "\t"; // REF

			// get alternative allele
			size_t nr_alleles = v.nr_of_alleles();
			if (nr_alleles < 2) {
				ostringstream oss;
				oss << "VariantReader::write_genotypes_of: less than 2 alleles given for variant at position " << v.get_start_position() << endl;
				throw runtime_error(oss.str());
			}

			vector<string> alt_alleles;
			vector<unsigned char> defined_alleles = {0};
			for (size_t i = 1; i < nr_alleles; ++i) {
				DnaSequence allele = v.get_allele_sequence(i);
				// skip alleles that are undefined
				if (!v.is_undefined_allele(i)) {
					alt_alleles.push_back(allele.to_string());
					defined_alleles.push_back(i);
				}
			}

			string alt_string = "";
			for (unsigned char a = 0; a < alt_alleles.size(); ++a) {
				if (a > 0) alt_string += ',';
				alt_string += alt_alleles[a];
			}

			this->genotyping_outfile[sample_index] << alt_string << "\t"; // ALT
			this->genotyping_outfile[sample_index] << ".\t"; // QUAL
			this->genotyping_outfile[sample_index] << "PASS" << "\t"; // FILTER
			// output allele frequencies of all alleles
			ostringstream info;
			info << "AF="; // AF
			for (unsigned int a = 1; a < defined_alleles.size(); ++a) {
				if (a > 1) info << ",";
				info << setprecision(6) << v.allele_frequency(defined_alleles[a], this->add_reference);				
			}

			// keep only likelihoods for genotypes with defined alleles
			size_t nr_missing = v.nr_missing_alleles();
			GenotypingResult genotype_likelihoods = singleton_likelihoods.at(j);
			if (nr_missing > 0) genotype_likelihoods = singleton_likelihoods.at(j).get_specific_likelihoods(defined_alleles);
			nr_alleles = defined_alleles.size();
            size_t stat_size=(defined_alleles.size()+2)* singleton_variants.size();
            size_t indexBase= stat_size* sample_index + (defined_alleles.size()+2)*j;
            unsigned short nr_uniq_kmers=variantsStatsPerSample[chromosome][i][indexBase];


			info << ";UK=" << nr_uniq_kmers; // UK
			info << ";AK="; // AK
			for (unsigned int a = 0; a < nr_alleles; ++a) {
				if (a > 0) info << ",";
				info << variantsStatsPerSample[chromosome][i][indexBase+2+a];
			}
			info << ";MA=" << nr_missing;
	
			// if IDs were given in input, write them to output as well
			if (!this->variant_ids[v.get_chromosome()][counter].empty()) info << ";ID=" << get_ids(v.get_chromosome(), alt_alleles, counter, false);
			this->genotyping_outfile[sample_index] << info.str() << "\t"; // INFO
			this->genotyping_outfile[sample_index] << "GT:GQ:GL:KC" << "\t"; // FORMAT

			// determine computed genotype
			pair<int,int> genotype = genotype_likelihoods.get_likeliest_genotype();
			if (ignore_imputed && (nr_uniq_kmers == 0)) genotype = {-1,-1};
			if ( (genotype.first != -1) && (genotype.second != -1)) {

				// unique maximum and therefore a likeliest genotype exists
				this->genotyping_outfile[sample_index] << genotype.first << "/" << genotype.second << ":"; // GT

				// output genotype quality
				this->genotyping_outfile[sample_index] << genotype_likelihoods.get_genotype_quality(genotype.first, genotype.second) << ":"; // GQ
			} else {
				// genotype could not be determined 
				this->genotyping_outfile[sample_index] << "./.:0:"; // GT:GQ
			}

			// output genotype likelihoods
			vector<long double> likelihoods = genotype_likelihoods.get_all_likelihoods(nr_alleles);
			if (likelihoods.size() < 3) {
				ostringstream oss;
				oss << "VariantReader::write_genotypes_of: too few likelihoods (" << likelihoods.size() << ") computed for variant at position " << v.get_start_position() << endl;
				throw runtime_error(oss.str());
			}

			ostringstream oss;
			oss << log10(likelihoods[0]);
			for (size_t j = 1; j < likelihoods.size(); ++j) {
				oss << "," << setprecision(4) << log10(likelihoods[j]);
			}
			this->genotyping_outfile[sample_index] << oss.str(); // GL
			this->genotyping_outfile[sample_index] << ":" << variantsStatsPerSample[chromosome][i][indexBase+1] << "\n"; // KC
			counter += 1;
		}
	}
}

void VariantReader::write_phasing_of(string chromosome, const vector<GenotypingResult>& genotyping_result, vector<UniqueKmers*>* unique_kmers, bool ignore_imputed) {
	// outfile needs to be open
	if (! this->phasing_outfile_open) {
		throw runtime_error("VariantReader::write_phasing_of: output file needs to be opened before writing.");
	}
	if (genotyping_result.size() != size_of(chromosome)) {
		throw runtime_error("VariantReader::write_phasing_of: number of variants and number of computed phasings differ.");
	}

	size_t counter = 0;
	for (size_t i = 0; i < size_of(chromosome); ++i) {
		Variant variant = this->variants_per_chromosome.at(chromosome)[i];

		// separate (possibly combined) variant into single variants and print a line for each
		vector<Variant> singleton_variants;
		vector<GenotypingResult> singleton_likelihoods;
		variant.separate_variants(&singleton_variants, &genotyping_result.at(i), &singleton_likelihoods);
		vector<VariantStats> singleton_stats;
		variant.variant_statistics(unique_kmers->at(i), singleton_stats);

		for (size_t j = 0; j < singleton_variants.size(); ++j) {
			Variant v = singleton_variants[j];
			v.remove_flanking_sequence();
			this->phasing_outfile << v.get_chromosome() << "\t"; // CHROM
			this->phasing_outfile << (v.get_start_position() + 1) << "\t"; // POS
			this->phasing_outfile << v.get_id() << "\t"; // ID
			this->phasing_outfile << v.get_allele_string(0) << "\t"; // REF

			// get alternative alleles
			size_t nr_alleles = v.nr_of_alleles();
			if (nr_alleles < 2) {
				ostringstream oss;
				oss << "VariantReader::write_phasing_of: less than 2 alleles given for variant at position " << v.get_start_position() << endl;
				throw runtime_error(oss.str());
			}

			vector<string> alt_alleles;
			vector<unsigned char> defined_alleles = {0};
			for (size_t i = 1; i < nr_alleles; ++i) {
				DnaSequence allele = v.get_allele_sequence(i);
				// skip alleles that are undefined
				if (! v.is_undefined_allele(i)) {
					alt_alleles.push_back(allele.to_string());
					defined_alleles.push_back(i);
				}
			}

			string alt_string = "";
			for (unsigned char a = 0; a < alt_alleles.size(); ++a) {
				if (a > 0) alt_string += ',';
				alt_string += alt_alleles[a];
			}

			size_t nr_missing = v.nr_missing_alleles();
			GenotypingResult genotype_likelihoods = singleton_likelihoods.at(j);
			if (nr_missing > 0) genotype_likelihoods = singleton_likelihoods.at(j).get_specific_likelihoods(defined_alleles);

			this->phasing_outfile << alt_string << "\t"; // ALT
			this->phasing_outfile << ".\t"; // QUAL
			this->phasing_outfile << "PASS" << "\t"; // FILTER
			// output allele frequencies of all alleles
			ostringstream info;
			info << "AF=";
			for (unsigned int a = 1; a < defined_alleles.size(); ++a) {
				if (a > 1) info << ",";
				info << setprecision(6) << v.allele_frequency(defined_alleles[a], this->add_reference);				
			}
			info << ";UK=" << singleton_stats.at(j).nr_unique_kmers; // UK
			info << ";AK="; // AK
			for (unsigned int a = 0; a < nr_alleles; ++a) {
				if (a > 0) info << ",";
				info << singleton_stats.at(j).kmer_counts[a];
			}

			info << ";MA=" << nr_missing;

			// if IDs were given in input, write them to output as well
			if (!this->variant_ids[v.get_chromosome()][counter].empty()) info << ";ID=" << get_ids(v.get_chromosome(), alt_alleles, counter, false);

			this->phasing_outfile << info.str() << "\t"; // INFO
			this->phasing_outfile << "GT:KC" << "\t"; // FORMAT

			// determine phasing
			if (ignore_imputed && (singleton_stats.at(j).nr_unique_kmers == 0)){
				this->phasing_outfile << "./."; // GT (phased)
			} else {
				pair<unsigned char,unsigned char> haplotype = singleton_likelihoods.at(j).get_haplotype();
				// check if the haplotype allele is undefined
				bool hap1_undefined = v.is_undefined_allele(haplotype.first);
				bool hap2_undefined = v.is_undefined_allele(haplotype.second);
				string hap1 (1, haplotype.first);
				string hap2 (1, haplotype.second);
				if (hap1_undefined) hap1 = ".";
				if (hap2_undefined) hap2 = ".";
				this->phasing_outfile << (unsigned int) haplotype.first << "|" << (unsigned int) haplotype.second; // GT (phased)
			}
			this->phasing_outfile << ":" << singleton_stats[j].coverage << endl; // KC
			counter += 1;
		}
	}
}

void VariantReader::close_genotyping_outfile() {
    for(unsigned i=0;i< samples.size();i++) {
        this->genotyping_outfile[i].close();
        this->genotyping_outfile_open = false;
    }
}

void VariantReader::close_phasing_outfile() {
	this->phasing_outfile.close();
	this->phasing_outfile_open = false;
}

size_t VariantReader::nr_of_genomic_kmers() const {
	return this->fasta_reader.get_total_kmers(this->kmer_size);
}

size_t VariantReader::nr_of_paths() const {
	return this->nr_paths;
}

void VariantReader::get_left_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const {
	size_t cur_start = this->variants_per_chromosome.at(chromosome).at(index).get_start_position();
	size_t prev_end = 0;
	if (index > 0) prev_end = this->variants_per_chromosome.at(chromosome).at(index-1).get_end_position();
	size_t overhang_start = cur_start - length;
	if (overhang_start < prev_end) overhang_start = prev_end;
	size_t overhang_end = cur_start;
	this->fasta_reader.get_subsequence(chromosome, overhang_start, overhang_end, result);
}

void VariantReader::get_right_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const {
	size_t cur_end = this->variants_per_chromosome.at(chromosome).at(index).get_end_position();
	size_t next_start = this->fasta_reader.get_size_of(chromosome);
	if (index < (this->variants_per_chromosome.at(chromosome).size() - 1)) next_start = this->variants_per_chromosome.at(chromosome).at(index+1).get_start_position();
	size_t overhang_end = cur_end + length;
	if (overhang_end > next_start) overhang_end = next_start;
	size_t overhang_start = cur_end;
	this->fasta_reader.get_subsequence(chromosome, overhang_start, overhang_end, result);
}

void VariantReader::setSampleName(std::string sampleName)
{
  samples.clear();
  samples.push_back(sampleName);
}
