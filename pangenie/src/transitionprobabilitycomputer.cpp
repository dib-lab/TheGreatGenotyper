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

    }
}


long double LiStephens::get(unsigned from_variant, unsigned to_variant,unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4){
    if(to_variant < from_variant) {
        swap(from_variant, to_variant);
        swap(path_id1, path_id3);
        swap(path_id2, path_id4);
    }
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

populationJointProbability::populationJointProbability(VariantReader* variants, std::string chromsome,std::vector<UniqueKmers*>* unique_kmers)
        : unique_kmers(unique_kmers)
{
    type="JointProbability";
    this->chromosome=chromsome;
    this->variants=variants;
}

populationJointProbability::populationJointProbability(VariantReader* variants, std::string chromsome, vector<EmissionProbabilities*> allemissions,std::vector<UniqueKmers*>* unique_kmers)
        : unique_kmers(unique_kmers)
{
    type="JointProbability";
    this->chromosome=chromsome;
    this->variants=variants;
    size_t nr_variants = this->variants->size_of(this->chromosome);
    this->probabilities = std::vector<std::vector<long double> >(nr_variants);

    int gts_per_state[100][100];
    #pragma omp parallel for firstprivate(gts_per_state)
    for (size_t v = 0; v < nr_variants-1; ++v) {
        const Variant& curr_variant = this->variants->get_variant(this->chromosome, v);
        const Variant& next_variant = this->variants->get_variant(this->chromosome, v+1);



        vector<unsigned char> curr_unique_alleles;
        (*unique_kmers)[v]->get_allele_ids(curr_unique_alleles);

        vector<unsigned char> next_unique_alleles;
        (*unique_kmers)[v+1]->get_allele_ids(next_unique_alleles);

        unsigned char curr_max_allele = (*unique_kmers)[v]->get_max_allele_id()+1;
        unsigned char next_max_allele = (*unique_kmers)[v+1]->get_max_allele_id()+1;

        const size_t curr_variant_pos = curr_variant.get_start_position();
        const size_t next_variant_pos = next_variant.get_start_position();
        long double recomb_rate=1.26;
        long double effective_N=25000.0L;
        auto curr_nr_paths = curr_variant.nr_of_paths();

//        long double distance = (next_variant_pos - curr_variant_pos) * 0.000004L * (recomb_rate) * effective_N;
//        // use Li-Stephans pair HMM transitions TODO: correct?
//        long double recomb_prob = (1.0L - exp(-distance /  (long double)curr_nr_paths) )* (1.0L /  (long double)curr_nr_paths);
//        long double no_recomb_prob = exp(-distance /  (long double)curr_nr_paths) + recomb_prob;
//

          long double distance = (next_variant_pos - curr_variant_pos);
          long double no_recomb_prob=exp( -0.00000001L*distance);
//        // use Li-Stephans pair HMM transitions TODO: correct?
//        long double recomb_prob = (1.0L - exp(-distance /  (long double)curr_nr_paths) )* (1.0L /  (long double)curr_nr_paths);
//        long double no_recomb_prob = exp(-distance /  (long double)curr_nr_paths) + recomb_prob;
//

        this->probabilities[v].resize(curr_max_allele* curr_max_allele* next_max_allele* next_max_allele);
        this->probabilities[v].assign(curr_max_allele* curr_max_allele* next_max_allele* next_max_allele,0.0);
        for(auto i: curr_unique_alleles)
            for(auto j: curr_unique_alleles)
            {
                gts_per_state[i][j]=0.0L;
            }
        for(auto emissions : allemissions) {
            for (unsigned i = 0; i < emissions->nr_samples; i++) {
                if (emissions->all_zeros[v][i] || emissions->all_zeros[v + 1][i])
                    continue;
// get likely and check quallity
//                long double curr_qual = emissions->gts_qual[v][i];
//                long double next_qual = emissions->gts_qual[v+1][i];
//                if(curr_qual < 0.7 || next_qual < 0.7)
//                    continue;

 //               pair<unsigned char,unsigned char> curr_gt = emissions->most_likely_gts[v][i];
//pair<unsigned char,unsigned char> next_gt = emissions->most_likely_gts[v+1][i];

                pair<unsigned char,unsigned char> curr_gt = emissions->result[i][v].get_likeliest_genotype();
                pair<unsigned char,unsigned char> next_gt = emissions->result[i][v+1].get_likeliest_genotype();


                size_t index = ((int)curr_gt.first * curr_max_allele * next_max_allele * next_max_allele) +
                               ((int)curr_gt.second * next_max_allele * next_max_allele) +
                               ((int)next_gt.first * next_max_allele) +
                               (int)next_gt.second;
                gts_per_state[(int)curr_gt.first][(int)curr_gt.second]+=1;
                this->probabilities[v][index] +=1;
            }
        }

//        for(auto i:curr_unique_alleles)
//            for(auto j : curr_unique_alleles)
//            {
//                cout<<int(i)<<"/"<<int(j)<<" : "<<gts_per_state[i][j]<<" , ";
//            }
//        cout<<"\n";

        for (auto c1 : curr_unique_alleles) {
            for (auto c2: curr_unique_alleles) {
                for (auto n1: next_unique_alleles) {
                    for (auto n2: next_unique_alleles) {
                        size_t index = ((int)c1 * curr_max_allele * next_max_allele * next_max_allele) +
                                       ((int)c2 * next_max_allele * next_max_allele) +
                                       ((int)n1 * next_max_allele) +
                                       (int)n2;
                        if(c1> c2)
                            swap(c1,c2);
                        if(n1 > n2)
                            swap(n1,n2);
                        size_t index_with_value = ((int)c1 * curr_max_allele * next_max_allele * next_max_allele) +
                                       ((int)c2 * next_max_allele * next_max_allele) +
                                       ((int)n1 * next_max_allele) +
                                       (int)n2;

                        this->probabilities[v][index] = this->probabilities[v][index_with_value];
                        if(gts_per_state[c1][c2] > 0.0L) {
                            this->probabilities[v][index] /= gts_per_state[c1][c2];
                        }
                        else{
                            this->probabilities[v][index]=0.0000000001L;
                        }
                    }
                }
            }
        }
//        for (auto c1 : curr_unique_alleles) {
//            for (auto c2 : curr_unique_alleles) {
//                for (auto n1 : next_unique_alleles) {
//                    for (auto n2 : next_unique_alleles) {
//                        long double jointPropSum=0;
//                        long double norm_a=0.00000000001L;
//                        long double norm_b=0.00000000001L;
//                        long double count=0;
//                        for(auto emissions : allemissions) {
//                            for (unsigned i = 0; i < emissions->nr_samples; i++) {
//                                if(emissions->all_zeros[v][i] || emissions->all_zeros[v+1][i])
//                                    continue;
//                                long double a=emissions->get_emission_probability(v,i,c1,c2,false);
//                                long double b=emissions->get_emission_probability(v+1,i,n1,n2,false);
//
//                                jointPropSum += (a*b);
//                                norm_a+= a*a;
//                                norm_b+= b*b;
//                            }
//                            count+=emissions->nr_samples;
//                        }
//
//                        jointPropSum /= (sqrt(norm_a) * sqrt(norm_b));
//                        size_t index = ((int)c1 * curr_max_allele * next_max_allele * next_max_allele) +
//                                       ((int)c2 * next_max_allele * next_max_allele) +
//                                       ((int)n1 * next_max_allele) +
//                                       (int)n2;
//
//                        this->probabilities[v][index]= jointPropSum * no_recomb_prob;
//                    }
//                }
//            }
//        }
//
//


    }

    normalize();
}
long double populationJointProbability::get(unsigned from_variant, unsigned to_variant,unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4)
{
    if(to_variant < from_variant) {
        swap(from_variant, to_variant);
        swap(path_id1, path_id3);
        swap(path_id2, path_id4);
    }
    const Variant& curr_variant = this->variants->get_variant(this->chromosome, from_variant);
    const Variant& next_variant = this->variants->get_variant(this->chromosome, to_variant);

    unsigned char from_allele1 =  curr_variant.get_allele_on_path(path_id1);
    unsigned char from_allele2 =  curr_variant.get_allele_on_path(path_id2);

    unsigned char to_allele1 =  next_variant.get_allele_on_path(path_id3);
    unsigned char to_allele2 =  next_variant.get_allele_on_path(path_id4);

    unsigned char from_max_allele = (*unique_kmers)[from_variant]->get_max_allele_id() +1;
    unsigned char to_max_allele = (*unique_kmers)[to_variant]->get_max_allele_id()+ 1;

    size_t index = ((int)from_allele1 * from_max_allele * to_max_allele * to_max_allele) +
                   ((int)from_allele2 * to_max_allele * to_max_allele) +
                   ((int)to_allele1 * to_max_allele) +
                   (int)to_allele2;
    return this->probabilities[from_variant][index];

}

void populationJointProbability::normalize()
{
    size_t nr_variants = this->variants->size_of(this->chromosome);
#pragma omp parallel for
    for (size_t v = 0; v < nr_variants-1; ++v) {
        const Variant& curr_variant = this->variants->get_variant(this->chromosome, v);
        const Variant& next_variant = this->variants->get_variant(this->chromosome, v+1);



        vector<unsigned char> curr_unique_alleles;
        (*unique_kmers)[v]->get_allele_ids(curr_unique_alleles);

        vector<unsigned char> next_unique_alleles;
        (*unique_kmers)[v+1]->get_allele_ids(next_unique_alleles);

        unsigned char curr_max_allele = (*unique_kmers)[v]->get_max_allele_id()+1;
        unsigned char next_max_allele = (*unique_kmers)[v+1]->get_max_allele_id()+1;

        this->probabilities[v].resize(curr_max_allele* curr_max_allele* next_max_allele* next_max_allele);
        for (auto c1 : curr_unique_alleles) {
            for (auto c2 : curr_unique_alleles) {
                long double sum=0;
                for (auto n1 : next_unique_alleles) {
                    for (auto n2 : next_unique_alleles) {
                        size_t index = ((int)c1 * curr_max_allele * next_max_allele * next_max_allele) +
                                       ((int)c2 * next_max_allele * next_max_allele) +
                                       ((int)n1 * next_max_allele) +
                                       (int)n2;
                        sum+=this->probabilities[v][index];
                    }
                }
                for (auto n1 : next_unique_alleles) {
                    for (auto n2 : next_unique_alleles) {
                        size_t index = ((int)c1 * curr_max_allele * next_max_allele * next_max_allele) +
                                       ((int)c2 * next_max_allele * next_max_allele) +
                                       ((int)n1 * next_max_allele) +
                                       (int)n2;
                        this->probabilities[v][index]/=sum;
                    }
                }
            }
        }




    }
}
