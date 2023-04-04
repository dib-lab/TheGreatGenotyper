#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cassert>
#include "emissionprobabilitycomputer.hpp"

using namespace std;

EmissionProbabilityComputer::EmissionProbabilityComputer(UniqueKmers* uniquekmers, ProbabilityTable* probabilities)
	:uniquekmers(uniquekmers),
	 probabilities(probabilities),
	 all_zeros(true)
{
	vector<unsigned char> unique_alleles;
	uniquekmers->get_allele_ids(unique_alleles);
	unsigned char max_allele = *max_element(std::begin(unique_alleles), std::end(unique_alleles));
	this->state_to_prob = vector< vector<long double>>(max_allele+1, vector<long double>(max_allele+1));



	for (auto a1 : unique_alleles) {
		for (auto a2 : unique_alleles) {
			bool a1_is_undefined = uniquekmers->is_undefined_allele(a1);
			bool a2_is_undefined = uniquekmers->is_undefined_allele(a2);
			this->state_to_prob[a1][a2] = compute_emission_probability(a1, a2, a1_is_undefined, a2_is_undefined);
			if (this->state_to_prob[a1][a2] > 0) this->all_zeros = false;
		}

	}

//	if (this->all_zeros) cerr << "EmissionProbabilities at position " << uniquekmers->get_variant_position() << " are all zero. Set to uniform." << endl;
}

long double EmissionProbabilityComputer::get_emission_probability(unsigned char allele_id1, unsigned char allele_id2) const {
	if (this->all_zeros) return 1.0L;
	return this->state_to_prob[allele_id1][allele_id2];
}

long double EmissionProbabilityComputer::compute_emission_probability(unsigned char allele_id1, unsigned char allele_id2, bool a1_undefined, bool a2_undefined){
	long double result = 1.0L;
	for (size_t i = 0; i < this->uniquekmers->size(); ++i){
		unsigned int expected_kmer_count = this->uniquekmers->alleles.at(allele_id1).first.get_position(i) + this->uniquekmers->alleles.at(allele_id2).first.get_position(i);
		if (a1_undefined && a2_undefined) {
			// all kmers can have copy numbers 0-2
			result *= (1.0L / 3.0L) * (this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(0) + this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(1) + this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(2));
		} else if (a1_undefined || a2_undefined) {
			// two possible copy numbers
			assert (expected_kmer_count < 2);
			result *= 0.5L * (this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(expected_kmer_count) + this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(expected_kmer_count + 1));
		} else {
			// expected kmer count is known
			result *= this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(expected_kmer_count);
		}
	}
	return result;
}
EmissionProbabilities::EmissionProbabilities()
{

    nr_samples=0;

//	if (this->all_zeros) cerr << "EmissionProbabilities at position " << uniquekmers->get_variant_position() << " are all zero. Set to uniform." << endl;
}

EmissionProbabilities::EmissionProbabilities(SamplesDatabase* samples,unsigned  nr_variants)
{

    this->state_to_prob = std::vector<std::vector<long double> > (nr_variants);
    this->numAllelesPerVariant = std::vector<unsigned short> (nr_variants);
    nr_samples= samples->getNumSamples();
    all_zeros=std::vector<std::vector<bool> >(nr_variants,std::vector<bool>(nr_samples,true));
    probabilities.resize(nr_samples);
    for(unsigned i=0;i< nr_samples;i++)
        probabilities[i]=samples->getSampleProbability(i);

//	if (this->all_zeros) cerr << "EmissionProbabilities at position " << uniquekmers->get_variant_position() << " are all zero. Set to uniform." << endl;
}

long double EmissionProbabilities::get_emission_probability(unsigned variantID,unsigned sampleID, unsigned char allele_id1, unsigned char allele_id2) const {
    if (this->all_zeros[variantID][sampleID]) return 1.0L;
    unsigned short max_allele= this->numAllelesPerVariant[variantID];
    unsigned size= (max_allele*(max_allele+1))/2;
    if(allele_id1 > allele_id2)
        swap(allele_id1,allele_id2);

    unsigned index=((int)allele_id1*(int)max_allele) - (((int)allele_id1-1)*(int)allele_id1)/2 + ((int)allele_id2-(int)allele_id1);
    index+=(sampleID*size);
    return this->state_to_prob[variantID][index];
}
size_t EmissionProbabilities::getNumVariants(){
    return all_zeros.size();

}
void EmissionProbabilities::compute(UniqueKmers* uniq,unsigned variantID,unsigned sampleID)
{
    vector<unsigned char> unique_alleles;
    uniq->get_allele_ids(unique_alleles);
    sort(unique_alleles.begin(),unique_alleles.end());

    unsigned char max_allele = unique_alleles[unique_alleles.size()-1] +1;
    unsigned size= (max_allele*(max_allele+1))/2;
    if(state_to_prob[variantID].size()==0) {
        this->state_to_prob[variantID].resize(size*nr_samples);
    }

    numAllelesPerVariant[variantID]=max_allele;

    for (unsigned i=0;i< unique_alleles.size();i++) {
        for (unsigned j=i; j< unique_alleles.size(); j ++) {
            auto a1=unique_alleles[i];
            auto a2=unique_alleles[j];
            bool a1_is_undefined = uniq->is_undefined_allele(a1);
            bool a2_is_undefined = uniq->is_undefined_allele(a2);
            unsigned index=((int)a1*(int)max_allele) - (((int)a1-1)*(int)a1)/2 + ((int)a2-(int)a1);
            index+=(sampleID*size);
            this->state_to_prob[variantID][index] = compute_emission_probability(uniq,sampleID,a1, a2, a1_is_undefined, a2_is_undefined);
            if (this->state_to_prob[variantID][index] > 0) this->all_zeros[variantID][sampleID] = false;
        }
    }


}



long double EmissionProbabilities::compute_emission_probability(UniqueKmers* uniquekmers,unsigned  sample_id,unsigned char allele_id1, unsigned char allele_id2, bool a1_undefined, bool a2_undefined){
    long double result = 1.0L;
    for (size_t i = 0; i < uniquekmers->size(); ++i){
        unsigned int expected_kmer_count = uniquekmers->alleles.at(allele_id1).first.get_position(i) + uniquekmers->alleles.at(allele_id2).first.get_position(i);
        if (a1_undefined && a2_undefined) {
            // all kmers can have copy numbers 0-2
            result *= (1.0L / 3.0L) * (this->probabilities[sample_id]->get_probability(uniquekmers->local_coverage, uniquekmers->kmer_to_count[i]).get_probability_of(0) + this->probabilities[sample_id]->get_probability(uniquekmers->local_coverage, uniquekmers->kmer_to_count[i]).get_probability_of(1) + this->probabilities[sample_id]->get_probability(uniquekmers->local_coverage, uniquekmers->kmer_to_count[i]).get_probability_of(2));
        } else if (a1_undefined || a2_undefined) {
            // two possible copy numbers
            assert (expected_kmer_count < 2);
            result *= 0.5L * (this->probabilities[sample_id]->get_probability(uniquekmers->local_coverage, uniquekmers->kmer_to_count[i]).get_probability_of(expected_kmer_count) + this->probabilities[sample_id]->get_probability(uniquekmers->local_coverage, uniquekmers->kmer_to_count[i]).get_probability_of(expected_kmer_count + 1));
        } else {
            // expected kmer count is known
            result *= this->probabilities[sample_id]->get_probability(uniquekmers->local_coverage, uniquekmers->kmer_to_count[i]).get_probability_of(expected_kmer_count);
        }
    }
    return result;
}

void EmissionProbabilities::save(std::string filename){
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    unsigned nr_variants = this->state_to_prob.size();
    ofs.write(reinterpret_cast<const char*>(&nr_variants), sizeof(nr_variants));
    for (const auto& row : this->state_to_prob) {
        size_t size = row.size();
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
        ofs.write(reinterpret_cast<const char*>(row.data()), size * sizeof(long double));
    }
    ofs.write(reinterpret_cast<const char*>(numAllelesPerVariant.data()), nr_variants * sizeof(unsigned short));
    ofs.write(reinterpret_cast<const char*>(&nr_samples), sizeof(nr_samples));
    vector<char> tmp(nr_samples);
    for (const auto& row : this->all_zeros)
    {
        for(unsigned  i =0 ; i< row.size();i++)
            tmp[i]=(char)row[i];
        ofs.write(reinterpret_cast<const char*>(tmp.data()), nr_samples );
    }



}
void EmissionProbabilities::load(std::string filename){
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        throw runtime_error("Failed to open file: "+filename);
        return ;
    }
    unsigned nr_variants;
    ifs.read(reinterpret_cast<char*>(&nr_variants), sizeof(nr_variants));
    this->state_to_prob.resize(nr_variants);
    for (auto& row : this->state_to_prob) {
        size_t size;
        ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
        row.resize(size);
        ifs.read(reinterpret_cast<char*>(row.data()), size * sizeof(long double));
    }

    numAllelesPerVariant.resize(nr_variants);
    ifs.read(reinterpret_cast<char*>(numAllelesPerVariant.data()), nr_variants * sizeof(unsigned short));
    size_t loadedNr_samples;
    ifs.read(reinterpret_cast<char*>(&loadedNr_samples), sizeof(loadedNr_samples));
    if(nr_samples!=0 && loadedNr_samples!=nr_samples)
    {
        throw runtime_error("file is corrupted: "+filename);
    }
    nr_samples=loadedNr_samples;
    all_zeros.resize(nr_variants);
    vector<char> tmp(nr_samples);
    for (auto& row : this->all_zeros)
    {
        ifs.read(reinterpret_cast<char*>(tmp.data()), nr_samples);
        row.resize(nr_samples);
        for(unsigned  i =0 ; i< row.size();i++)
            row[i]=(bool)tmp[i];
    }


}
void EmissionProbabilities::destroy()
{
    this->state_to_prob = std::vector<std::vector<long double> > (1);
    this->numAllelesPerVariant = std::vector<unsigned short> (1);
    all_zeros=std::vector<std::vector<bool> >(1);

}
