#ifndef TRANSITIONPROBABILITYCOMPUTER_HPP
#define TRANSITIONPROBABILITYCOMPUTER_HPP

#include "uniquekmers.hpp"
#include "variantreader.hpp"

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
    TransitionProbability(VariantReader* variants,std::string chromsome);
    void computeLiStephens(double recomb_rate, long double effective_N = 25000.0L);
    long double get(unsigned from_variant, unsigned to_variant,unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4);

private:
    VariantReader* variants;
    std::string chromosome;
    std::vector< std::vector<long double> > probabilities;
};




#endif // TRANSITIONPROBABILITYCOMPUTER_HPP
