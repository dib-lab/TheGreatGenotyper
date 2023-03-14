#include <cassert>
#include <math.h>
#include "transitionprobabilitycomputer.hpp"
#include <iostream>

using namespace std;

TransitionProbabilityComputer::TransitionProbabilityComputer(size_t from_variant, size_t to_variant, double recomb_rate, unsigned short nr_paths,bool phased, bool uniform, long double effective_N)
{
	assert(from_variant <= to_variant );
	this->uniform = uniform;
    this->phased = phased;
	// using same formula as in WhatsHap
//	long double distance = (to_variant - from_variant) * 0.000001 * ((long double) recomb_rate) * 4.0L * effective_N;
	long double distance = (to_variant - from_variant) * 0.000004L * ((long double) recomb_rate) * effective_N;
	// use Li-Stephans pair HMM transitions TODO: correct?
	long double recomb_prob = (1.0L - exp(-distance / (long double) nr_paths) )* (1.0L / (long double) nr_paths);
	long double no_recomb_prob = exp(-distance / (long double) nr_paths) + recomb_prob;
	this->probabilities = {no_recomb_prob*no_recomb_prob, no_recomb_prob*recomb_prob, recomb_prob*recomb_prob};
}

long double TransitionProbabilityComputer::compute_transition_prob(unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4){
	if (this->uniform) {
		return 1.0L;
	} else if (this->phased){
		// determine number of recombination events
		unsigned int nr_events = 0;
		if (path_id1 != path_id3) nr_events += 1;
		if (path_id2 != path_id4) nr_events += 1;
		return this->probabilities[nr_events];
	}
    else{
        unsigned int nr_events = 0;
        if (path_id1 == path_id3)
        {
            if(path_id2 != path_id4)
                nr_events+=1;
        }
        else if (path_id1 == path_id4)
        {
            if(path_id2 != path_id3)
                nr_events+=1;
        }
        else
        {
            nr_events+=1;
            if(path_id2 != path_id3 and path_id2 != path_id4  )
                nr_events+=1;
        }
        return this->probabilities[nr_events];
    }
}
