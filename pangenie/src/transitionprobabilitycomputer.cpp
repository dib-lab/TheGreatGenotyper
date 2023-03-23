#include <cassert>
#include <math.h>
#include "transitionprobabilitycomputer.hpp"
#include <iostream>
#include <fstream>

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
{}

LiStephens::LiStephens(VariantReader* variants,std::string chromosome,double recomb_rate, long double effective_N)
        :TransitionProbability(variants,chromosome )
{
    type="LiStepens";
    size_t nr_variants = this->variants->size_of(this->chromosome);
    this->probabilities= std::vector<std::vector<long double> >(nr_variants);
#pragma omp parallel for
    for (size_t v = 0; v < nr_variants-1; ++v) {
        const Variant& curr_variant = this->variants->get_variant(this->chromosome, v);
        const Variant& next_variant = this->variants->get_variant(this->chromosome, v+1);

        const size_t curr_variant_pos = curr_variant.get_start_position();
        const size_t next_variant_pos = next_variant.get_start_position();
        size_t curr_nr_paths= (long double) curr_variant.nr_of_paths();
        size_t next_nr_paths= (long double) next_variant.nr_of_paths();


        long double distance = (next_variant_pos - curr_variant_pos) * 0.000004L * ((long double) recomb_rate) * effective_N;
        // use Li-Stephans pair HMM transitions TODO: correct?
        long double recomb_prob = (1.0L - exp(-distance /  (long double)curr_nr_paths) )* (1.0L /  (long double)curr_nr_paths);
        long double no_recomb_prob = exp(-distance /  (long double)curr_nr_paths) + recomb_prob;

        this->probabilities[v] = {no_recomb_prob*no_recomb_prob, no_recomb_prob*recomb_prob, recomb_prob*recomb_prob};

        this->probabilities[v]= std::vector<long double>(curr_nr_paths*next_nr_paths*curr_nr_paths*next_nr_paths);
    }
}


long double LiStephens::get(unsigned from_variant, unsigned to_variant,unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4){
    if(to_variant < from_variant)
        swap(from_variant,to_variant);
    const Variant& curr_variant = this->variants->get_variant(this->chromosome, from_variant);
    const Variant& next_variant = this->variants->get_variant(this->chromosome, to_variant);
    bool phased =  curr_variant.get_phase_status() && next_variant.get_phase_status();
    unsigned int nr_events = 0;
    if (phased){

        // determine number of recombination events

        if (path_id1 != path_id3) nr_events += 1;
        if (path_id2 != path_id4) nr_events += 1;
        //		cout<<"( "<<path_id1<<", "<<path_id2<<")"<<"\t( "<<path_id3<<", "<<path_id4<<")  -> "<<nr_events<<endl;

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
    }
    return this->probabilities[from_variant][nr_events];

}

void TransitionProbability::save(std::string filename){
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    size_t size = this->probabilities.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto& row : this->probabilities) {
        size = row.size();
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
        ofs.write(reinterpret_cast<const char*>(row.data()), size * sizeof(long double));
    }
}
void TransitionProbability::load(std::string filename){
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        throw runtime_error("Failed to open file: "+filename);
        return ;
    }
    size_t size;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    this->probabilities.resize(size);
    for (auto& row : this->probabilities) {
        ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
        row.resize(size);
        ifs.read(reinterpret_cast<char*>(row.data()), size * sizeof(long double));
    }
}