#ifndef TRANSITIONPROBABILITYCOMPUTER_HPP
#define TRANSITIONPROBABILITYCOMPUTER_HPP

#include "uniquekmers.hpp"
#include "variantreader.hpp"
#include "emissionprobabilitycomputer.hpp"

/** 
* Computes the transition probabilities between variants.
**/

class TransitionProbabilityComputer {
public:
	TransitionProbabilityComputer(size_t from_variant, size_t to_variant, double recomb_rate, unsigned short nr_paths,bool phased =true, bool uniform = false, long double effective_N = 25000.0L);
	long double compute_transition_prob(unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4);
private:
	std::vector<long double> probabilities;
	bool uniform;
    bool phased;
	
};
class TransitionProbability{
public:
    TransitionProbability(){}
    TransitionProbability(VariantReader* variants,std::string chromsome);
    virtual long double get(unsigned from_variant, unsigned to_variant,unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4)=0;
    virtual long double get(unsigned from_variant, unsigned to_variant,unsigned nr_switches)=0;
    void save(std::string filename);
    void load(std::string filename);
protected:
    VariantReader* variants;
    std::string chromosome;
    std::vector< std::vector<long double> > probabilities;
    std::string type;
};

class LiStephens: public TransitionProbability{
public:
    LiStephens(VariantReader* variants,std::string chromsome,double recomb_rate, long double effective_N = 25000.0L);
    long double get(unsigned from_variant, unsigned to_variant,unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4);
    long double get(unsigned from_variant, unsigned to_variant,unsigned nr_switches);
};

class populationJointProbability: public TransitionProbability{
public:
    populationJointProbability(VariantReader* variants, std::string chromsome,std::vector<UniqueKmers*>* unique_kmers);
	populationJointProbability(VariantReader* variants, std::string chromsome, vector<EmissionProbabilities*> emissions,std::vector<UniqueKmers*>* unique_kmers);
	long double get(unsigned from_variant, unsigned to_variant,unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4);
    void normalize();
private:
	std::vector<UniqueKmers*>* unique_kmers;
};




#endif // TRANSITIONPROBABILITYCOMPUTER_HPP
