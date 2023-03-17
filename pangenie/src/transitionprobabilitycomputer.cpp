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
		//		cout<<"( "<<path_id1<<", "<<path_id2<<")"<<"\t( "<<path_id3<<", "<<path_id4<<")  -> "<<nr_events<<endl;
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
	//cout<<"( "<<path_id1<<", "<<path_id2<<")"<<"\t( "<<path_id3<<", "<<path_id4<<")  -> "<<nr_events<<endl;
        return this->probabilities[nr_events];
    }
}


TransitionProbability::TransitionProbability(VariantReader* variants,std::string chromosome):
variants(variants),
chromosome(chromosome)
{
    this->probabilities= std::vector<std::vector<long double> >(variants->size_of(chromosome));
}

void TransitionProbability::computeLiStephens(double recomb_rate, long double effective_N){
    size_t nr_variants = this->variants->size_of(this->chromosome);
#pragma omp parallel for
    for (size_t v = 0; v < nr_variants-1; ++v) {
        const Variant& curr_variant = this->variants->get_variant(this->chromosome, v);
        const Variant& next_variant = this->variants->get_variant(this->chromosome, v+1);

        const size_t curr_variant_pos = curr_variant.get_start_position();
        const size_t next_variant_pos = next_variant.get_start_position();
        size_t curr_nr_paths= (long double) curr_variant.nr_of_paths();
        size_t next_nr_paths= (long double) next_variant.nr_of_paths();

        bool phased= curr_variant.get_phase_status() && next_variant.get_phase_status();

        long double distance = (next_variant_pos - curr_variant_pos) * 0.000004L * ((long double) recomb_rate) * effective_N;
        // use Li-Stephans pair HMM transitions TODO: correct?
        long double recomb_prob = (1.0L - exp(-distance /  (long double)curr_nr_paths) )* (1.0L /  (long double)curr_nr_paths);
        long double no_recomb_prob = exp(-distance /  (long double)curr_nr_paths) + recomb_prob;

        std::vector<long double> props = {no_recomb_prob*no_recomb_prob, no_recomb_prob*recomb_prob, recomb_prob*recomb_prob};

        this->probabilities[v]= std::vector<long double>(curr_nr_paths*next_nr_paths*curr_nr_paths*next_nr_paths);
        size_t index= 0;
        unsigned int nr_events = 0;
        for(unsigned path1= 0; path1 < curr_nr_paths; path1++)
            for(unsigned path2= 0; path2 < curr_nr_paths; path2++)
                for(unsigned path3= 0; path3 < next_nr_paths; path3++)
                    for(unsigned path4= 0; path4 < next_nr_paths; path4++)
                    {
                        nr_events = 0;
                        if(phased) {
                            if (path1 != path3) nr_events += 1;
                            if (path2 != path4) nr_events += 1;
                            this->probabilities[v][index] = props[nr_events];
                        }
                        else{
                            if (path1 == path3)
                            {
                                if(path2 != path4)
                                    nr_events+=1;
                            }
                            else if (path1 == path4)
                            {
                                if(path2 != path3)
                                    nr_events+=1;
                            }
                            else
                            {
                                nr_events+=1;
                                if(path2 != path3 and path2 != path4  )
                                    nr_events+=1;
                            }
                        }
                        index++;
                    }

    }
}

long double TransitionProbability::get(unsigned from_variant, unsigned to_variant,unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4){
    if(to_variant < from_variant)
        swap(from_variant,to_variant);
    const Variant& curr_variant = this->variants->get_variant(this->chromosome, from_variant);
    const Variant& next_variant = this->variants->get_variant(this->chromosome, to_variant);

    size_t curr_nr_paths= (long double) curr_variant.nr_of_paths();
    size_t next_nr_paths= (long double) next_variant.nr_of_paths();

    size_t index = (path_id1 * curr_nr_paths * next_nr_paths * next_nr_paths) +
                   (path_id2 * next_nr_paths * next_nr_paths) +
                   (path_id3 * next_nr_paths) +
                    path_id4;

    return this->probabilities[from_variant][index];

}
