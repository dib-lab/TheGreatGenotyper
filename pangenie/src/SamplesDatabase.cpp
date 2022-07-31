//
// Created by Mostafa on 7/28/2022.
//

#include "SamplesDatabase.h"

#include "cli/load/load_graph.cpp"
#include "cli/stats.hpp"
#include "cli/query.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
SamplesDatabase::SamplesDatabase(){

}
SamplesDatabase::SamplesDatabase(string graph_path, string annotation_path, string descriptionFile,long double regularization){
  graph=mtg::cli::load_critical_dbg(graph_path);
  kSize=graph->get_k();
  auto config = std::make_unique<mtg::cli::Config>();
  config->query_counts=true;
  config->infbase=graph_path;
  config->infbase_annotators.push_back(annotation_path);
  anno_graph=mtg::cli::initialize_annotated_dbg(graph, *config);

  ifstream inputDescriptor(descriptionFile);
  uint32_t index=0;

  string sample_name,sampleLabel;
  int coverage;
  while(inputDescriptor>> sample_name>> coverage >> sampleLabel)
  {
    samples.push_back(SampleDescriptor(sample_name,sampleLabel,coverage,regularization));
    metaLabel_to_Index[sampleLabel]=index++;
  }
}

size_t SamplesDatabase::getKSize()
{
  return kSize;
}

size_t SamplesDatabase::getNumSamples()
{
  return samples.size();
}
int SamplesDatabase::getKmerCoverage(unsigned sampleIndex)
{
  return samples[sampleIndex].kmerCoverage;
}
string SamplesDatabase::getSampleName(unsigned sampleIndex)
{
  return samples[sampleIndex].SampleName;
}
ProbabilityTable* SamplesDatabase::getSampleProbability(unsigned sampleIndex)
{
  return &(samples[sampleIndex].probs);
}
void SamplesDatabase::getKmerCounts(vector<string>& seqs,vector<unordered_map<string,uint32_t>> &kmerCounts)
{
  kmerCounts.resize(samples.size());
  float discovery_fraction=0.1;
  float  presence_fraction=0.0;
  for(auto sequence: seqs)
  {
    if(sequence.size()<kSize)
      continue;

    LabelCountAbundancesVec result= anno_graph->get_kmer_counts(sequence,
                                                                samples.size(),
                                                                discovery_fraction,
                                                                 presence_fraction);
    size_t numKmers=sequence.size()+1-kSize;

    vector<string> kmers(numKmers);
    for(unsigned i=0; i < numKmers;i++)
    {
      kmers[i]=sequence.substr(i,kSize);
    }
    for(auto t: result)
    {
        string label=std::get<0>(t);
        auto counts=std::get<2>(t);
        if(numKmers != counts.size())
        {
          throw std::logic_error("number of kmers doesnt match query result");
        }
	auto it=metaLabel_to_Index.find(label);
	if(it == metaLabel_to_Index.end())
	{
	  throw std::logic_error("Sample "+label+" doesnt exist in graph description");
	}
        uint32_t sampleIndex= it->second;
        for(unsigned i = 0; i< numKmers; i++)
        {
          kmerCounts[sampleIndex][kmers[i]] = counts[i];
        }

    }

  }

}
