//
// Created by Mostafa on 7/28/2022.
//

#ifndef THEGREATGENOTYPER_SAMPLESDATABASE_H
#define THEGREATGENOTYPER_SAMPLESDATABASE_H

#include "graph/annotated_dbg.hpp"
#include <vector>
#include "probabilitytable.hpp"
#include <unordered_map>
using namespace std;

struct SampleDescriptor{
  string SampleName;
  string metagraphLabel;
  int kmerCoverage;
  ProbabilityTable probs;
  SampleDescriptor(string SampleName,string metagraphLabel,int kmerCoverage,long double regularization):
                                                                                                             SampleName(SampleName),metagraphLabel(metagraphLabel),kmerCoverage(kmerCoverage)
  {
    probs=ProbabilityTable(kmerCoverage / 4, kmerCoverage*4, 2*kmerCoverage, regularization);
  }
};

class SamplesDatabase {
public:
  SamplesDatabase();
  SamplesDatabase(string graphPath, string annotation_path, string descPath,long double regularization,bool log_scale);
  int getKmerCoverage(unsigned sampleIndex);
  size_t getKSize();
  size_t getNumSamples();
  string getSampleName(unsigned sampleIndex);
  vector<string> getSamplesName();
  void load_graph();
  void delete_graph();
  ProbabilityTable* getSampleProbability(unsigned sampleIndex);
  void getKmerCounts(vector<string>& seqs,vector<unordered_map<string,uint32_t>> & kmerCounts);
private:
  bool log_scale;
  size_t kSize;
  typedef std::vector<std::tuple<string, size_t, std::vector<size_t>>> LabelCountAbundancesVec;
  std::shared_ptr<mtg::graph::DeBruijnGraph> graph;
  std::unique_ptr<mtg::graph::AnnotatedDBG> anno_graph;
  vector<SampleDescriptor> samples;
  unordered_map<string,uint32_t> metaLabel_to_Index;
  string graph_path;
  string anno_path;
};

#endif // THEGREATGENOTYPER_SAMPLESDATABASE_H
