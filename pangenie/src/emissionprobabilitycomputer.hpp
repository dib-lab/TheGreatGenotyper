#ifndef EMISSIONPROBABILITYCOMPUTER_H
#define EMISSIONPROBABILITYCOMPUTER_H

#include <vector>
#include <string>
#include <unordered_map>
#include "uniquekmers.hpp"
#include "copynumber.hpp"
#include "columnindexer.hpp"
#include "probabilitytable.hpp"
#include "SamplesDatabase.h"
#include "genotypingresult.hpp"
/** 
* Computes the emission probabilities for a variant position.
**/

typedef std::vector<std::vector<long double>> ProbabilityMatrix;

class EmissionProbabilityComputer {
public:
	/**
	* @param uniquekmers all unique kmers for this position
	 **/
	EmissionProbabilityComputer(UniqueKmers* uniquekmers, ProbabilityTable* probabilities);
	/** get emission probability for a state in the HMM **/
	long double get_emission_probability(unsigned char allele_id1, unsigned char allele_id2) const;

private:
	UniqueKmers* uniquekmers;
	ProbabilityTable* probabilities;
	bool all_zeros;
	ProbabilityMatrix state_to_prob;
	long double compute_emission_probability(unsigned char allele1, unsigned char allele2, bool allele1_undefined, bool allele2_undefined);
};

typedef std::vector<std::vector<std::vector<std::vector<long double> > > > ProbabilityMatrices;


class EmissionProbabilities {
public:
    /**
    * @param uniquekmers all unique kmers for this position
     **/
    EmissionProbabilities();
    EmissionProbabilities(SamplesDatabase* samples,unsigned  nr_variants);
    /** get emission probability for a state in the HMM **/
    long double get_emission_probability(unsigned variantID,unsigned sampleID, unsigned char allele_id1, unsigned char allele_id2,bool return_one_if_all_zeros=true) const;
    void compute(UniqueKmers* uniq,unsigned variantID,unsigned sampleID);
    void compute_most_likely_genotypes(std::vector<UniqueKmers*>* uniqKmers);
    size_t getNumVariants();
    void save(std::string filename);
    void load(std::string filename);
    void destroy();
    vector<vector<GenotypingResult> > result;
private:
    std::vector<ProbabilityTable*> probabilities;
    std::vector<std::vector<bool> > all_zeros;
    std::vector<std::vector<long double> >  state_to_prob;
    std::vector<unsigned short> numAllelesPerVariant;
    std::vector<std::vector<std::pair<unsigned char, unsigned char> > > most_likely_gts;
    std::vector<std::vector<long double > > gts_qual;
    unsigned nr_samples;
    long double compute_emission_probability(UniqueKmers* uniquekmers,unsigned  sample_id,unsigned char allele1, unsigned char allele2, bool allele1_undefined, bool allele2_undefined);

	friend class populationJointProbability;
};

# endif // EMISSIONPROBABILITYCOMPUTER_H
