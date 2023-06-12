#include "uniquekmercomputer.hpp"
#include <jellyfish/mer_dna.hpp>
#include <iostream>
#include <cassert>
#include <map>

using namespace std;

void unique_kmers(DnaSequence& allele, unsigned char index, size_t kmer_size, map<jellyfish::mer_dna, vector<unsigned char>>& occurences) {
	//enumerate kmers
	map<jellyfish::mer_dna, size_t> counts;
	size_t extra_shifts = kmer_size;
	jellyfish::mer_dna::k(kmer_size);
	jellyfish::mer_dna current_kmer("");
	for (size_t i = 0; i < allele.size(); ++i) {
		char current_base = allele[i];
		if (extra_shifts == 0) {
			counts[current_kmer] += 1;
		}
		if (  ( current_base != 'A') && (current_base != 'C') && (current_base != 'G') && (current_base != 'T') ) {
			extra_shifts = kmer_size + 1;
		}
		current_kmer.shift_left(current_base);
		if (extra_shifts > 0) extra_shifts -= 1;
	}
	counts[current_kmer] += 1;

	// determine kmers unique to allele
	for (auto const& entry : counts) {
		if (entry.second == 1) occurences[entry.first].push_back(index);
	}
}

UniqueKmerComputer::UniqueKmerComputer (KmerCounter* genomic_kmers,  VariantReader* variants, string chromosome)
	:genomic_kmers(genomic_kmers),
	 variants(variants),
	 chromosome(chromosome)
{
	jellyfish::mer_dna::k(this->variants->get_kmer_size());
    compute_unique_kmers();
}

void UniqueKmerComputer::compute_unique_kmers() {
	size_t nr_variants = this->variants->size_of(this->chromosome);
    uniqKmers.resize(nr_variants);
#pragma omp parallel for shared(uniqKmers)
    for (size_t v = 0; v < nr_variants; ++v) {
        // set parameters of distributions
        size_t kmer_size = this->variants->get_kmer_size();



        const Variant& variant = this->variants->get_variant(this->chromosome, v);

        UniqueKmers* u = new UniqueKmers(variant.get_start_position(),variant.get_phase_status());
        size_t nr_alleles = variant.nr_of_alleles();
        // insert empty alleles (to also capture paths for which no unique kmers exist)
        assert(variant.nr_of_paths() < 65535);
        for (unsigned short p = 0; p < variant.nr_of_paths(); ++p) {
            unsigned char a = variant.get_allele_on_path(p);
            u->insert_empty_allele(a);
            u->insert_path(p,a);
        }


        for (unsigned char a = 0; a < nr_alleles; ++a) {
            // consider all alleles not undefined
            if (variant.is_undefined_allele(a)) {
                // skip kmers of alleles that are undefined
                u->set_undefined_allele(a);
                continue;
            }
            DnaSequence allele = variant.get_allele_sequence(a);
            unique_kmers(allele, a, kmer_size, u->occurences);
        }

        for (auto kmer = u->occurences.cbegin(), next_it = kmer; kmer != u->occurences.cend(); kmer = next_it)
        {
            bool ignoreKmer=false;
            size_t genomic_count = this->genomic_kmers->getKmerAbundance(kmer->first);
            size_t local_count = kmer->second.size();

            if ( (genomic_count - local_count) == 0 ) {
                // kmer unique to this region
                // determine read kmercount for this kmer

                // determine on which paths kmer occurs
                vector<size_t> paths;
                for (auto &allele : kmer->second) {
                    variant.get_paths_of_allele(allele, paths);
                }

                // skip kmer that does not occur on any path (uncovered allele)
                if (paths.size() == 0) {
                    ignoreKmer=true;
                }

                // skip kmer that occurs on all paths (they do not give any information about a genotype)
                if (paths.size() == variant.nr_of_paths()) {
                    ignoreKmer=true;
                }
            }
            else{
                ignoreKmer=true;
            }


            ++next_it;
            if (ignoreKmer)
            {
                u->occurences.erase(kmer);
            }
        }


        uniqKmers[v]=u;

        //delete u;
    }

}

void UniqueKmerComputer::compute_emissions(SamplesDatabase* database, EmissionProbabilities* result){
    size_t nr_variants = this->variants->size_of(this->chromosome);
    uint32_t numSamples=database->getNumSamples();
    vector<double> localCoverage(numSamples);
    vector<string> seqs;
    seqs.resize(100);
    cerr<<"Enter"<<endl;
    cerr<<"number of samples "<< numSamples<<endl;
//#pragma omp parallel for shared(result,uniqKmers) firstprivate(localCoverage,numSamples,seqs)
    for (size_t v = 0; v < nr_variants; ++v) {
        // set parameters of distributions
        size_t kmer_size = this->variants->get_kmer_size();

        compute_local_coverage(database,this->chromosome, v, 2 * kmer_size, localCoverage);

        const Variant &variant = this->variants->get_variant(this->chromosome, v);
        size_t nr_alleles = variant.nr_of_alleles();
        seqs.resize(nr_alleles);
        vector<unsigned char> defined_alleles;
        cerr<<variant.get_start_position()<<"-"<<variant.get_end_position()<<endl;
        for (unsigned char a = 0; a < nr_alleles; ++a) {
            // consider all alleles not undefined
            if (variant.is_undefined_allele(a)) {
                // skip kmers of alleles that are undefined
                continue;
            }
            DnaSequence allele = variant.get_allele_sequence(a);
            seqs[a]=allele.to_string();
            defined_alleles.push_back(a);
        }
        //cerr<<"point 1"<<endl;
        vector<unordered_map<string,uint32_t>> kmerCounts;
        database->getKmerCounts(seqs,kmerCounts);
        //cerr<<"point 2"<<endl;
        //cerr<<"kmer counts size "<<kmerCounts.size()<<endl;

        for(unsigned sampleID=0; sampleID< numSamples ; sampleID++) {
            cerr<<sampleID<<endl;
            UniqueKmers* sampleU= new UniqueKmers(*(uniqKmers[v]));
            string sampleName=database->getSampleName(sampleID);
            //bool debug= (sampleName=="SRR17029944");

            unsigned sampleKmerCoverage=database->getKmerCoverage(sampleID);
            sampleU->set_coverage(localCoverage[sampleID]);
            // if(debug) cerr<<"Local Coverage "<<localCoverage[sampleID]<<"\n";
            //		  debug=false;
            ProbabilityTable* probabilities=database->getSampleProbability(sampleID);
            size_t nr_kmers_used = 0;
            for (auto &kmer : sampleU->occurences) {
                if (nr_kmers_used > 300)
                    break;
                auto it =kmerCounts[sampleID].find(kmer.first.to_str());
                size_t read_kmercount =0;
                if(it != kmerCounts[sampleID].end() ) {
                    read_kmercount =it->second;
                }
                if (read_kmercount >
                    (2 * sampleKmerCoverage)) {
                    continue;
                }

                // determine probabilities
                CopyNumber cn = probabilities->get_probability(
                        localCoverage[sampleID], read_kmercount);

                long double p_cn0 = cn.get_probability_of(0);
                long double p_cn1 = cn.get_probability_of(1);
                long double p_cn2 = cn.get_probability_of(2);

                // skip kmers with only 0 probabilities
                if ((p_cn0 > 0) || (p_cn1 > 0) || (p_cn2 > 0)) {
                    nr_kmers_used += 1;
                    sampleU->insert_kmer(read_kmercount, kmer.second);
                    //  if(debug)
                    //	cerr<<"KMERFN"<<"\t"<<kmer.first<<"\t"<<read_kmercount<<"\n";
                }
            }
          //  cerr<<"point 4"<<endl;
            vector<Variant> singleton_variants;
            variant.separate_variants(&singleton_variants);
            vector<VariantStats> singleton_stats;
            variant.variant_statistics(sampleU, singleton_stats);
            variants->addVariantStat(v, sampleName,this->chromosome, defined_alleles,singleton_stats);
            result->compute(sampleU,v,sampleID);
            delete sampleU;
            //cerr<<"point 5"<<endl;
        }
    }
    cerr<<"point 6"<<endl;
}

void UniqueKmerComputer::compute_empty(vector<UniqueKmers*>* result) const {
	size_t nr_variants = this->variants->size_of(this->chromosome);
	for (size_t v = 0; v < nr_variants; ++v) {
		const Variant& variant = this->variants->get_variant(this->chromosome, v);
		UniqueKmers* u = new UniqueKmers(variant.get_start_position());
		size_t nr_alleles = variant.nr_of_alleles();

		// insert empty alleles and paths
		assert(variant.nr_of_paths() < 65535);
		for (unsigned short p = 0; p < variant.nr_of_paths(); ++p) {
			unsigned char a = variant.get_allele_on_path(p);
			u->insert_empty_allele(a);
			u->insert_path(p,a);
		}
		result->push_back(u);
	}
}

void UniqueKmerComputer::compute_local_coverage(SamplesDatabase* database,string chromosome, size_t var_index, size_t length,vector<double>& result) {
        DnaSequence left_overhang;
	DnaSequence right_overhang;
	
        uint32_t numSamples=database->getNumSamples();
        result.resize(numSamples);


        this->variants->get_left_overhang(chromosome, var_index, length, left_overhang);
	this->variants->get_right_overhang(chromosome, var_index, length, right_overhang);

	size_t kmer_size = this->variants->get_kmer_size();
	map <jellyfish::mer_dna, vector<unsigned char>> occurences;
	unique_kmers(left_overhang, 0, kmer_size, occurences);
	unique_kmers(right_overhang, 1, kmer_size, occurences);

        vector<string> seqs;
        seqs.push_back(left_overhang.to_string());
        seqs.push_back(right_overhang.to_string());
        vector<unordered_map<string,uint32_t>> kmerCounts;

        database->getKmerCounts(seqs,kmerCounts);

        for(unsigned sampleID=0; sampleID<numSamples ;sampleID++) {

          size_t sample_kmer_coverage = database->getKmerCoverage(sampleID);
	  size_t total_coverage = 0;
	  size_t total_kmers = 0;
	  string sampleName=database->getSampleName(sampleID);
	  //bool debug= (sampleName=="SRR17029944");
	  

          for (auto &kmer : occurences) {
            size_t genomic_count =
                this->genomic_kmers->getKmerAbundance(kmer.first);

            if (genomic_count == 1) {
              size_t read_count =
                  kmerCounts[sampleID][kmer.first.to_str()];
	      //	      if(debug)
	      //cerr<<"KMER"<<"\t"<<kmer.first.to_str()<<"\t"<<read_count<<"\n";
              // ignore too extreme counts
              if ((read_count < (sample_kmer_coverage / 4)) ||
                  (read_count > (sample_kmer_coverage * 4)))
                continue;

              total_coverage += read_count;
              total_kmers += 1;
            }
          }
          // in case no unique kmers were found, use constant kmer coverage
          if ((total_kmers > 0) && (total_coverage > 0)) {
             result[sampleID] = total_coverage / total_kmers;
          } else {
            result[sampleID] =  sample_kmer_coverage;
          }
        }
}
