//
// Created by Mostafa on 7/26/2022.
//
#include <iostream>
#include <fstream>
#include <kmc_file.h>
#include "CLI11.hpp"
#include "histogram.hpp"


using namespace std;

int main(int argc, char *argv[]) {
  CLI::App app;
  string kmcInput;
  string outputfile;
  uint32_t max_count = 1000;
  bool largest_peak = true;

  app.add_option("-i,--input", kmcInput, "prefix for KMC files")->required();

  app.add_option("-m,--max-count", max_count, "Max Count for the kmers");

  app.add_option("-o,--output", outputfile,
                 "Output Path for the histogram file")
      ->required();

  CLI11_PARSE(app, argc, argv);

  uint32 _kmer_length;
  uint32 _mode;
  uint32 _counter_size;
  uint32 _lut_prefix_length;
  uint32 _signature_len;
  uint32 _min_count;
  uint64 _max_count;
  uint64 _total_kmers;
  CKMCFile kmer_data_base;
  if (!kmer_data_base.OpenForListing(kmcInput)) {
    cerr << "Cant open KMC DB " << kmcInput << endl;
    return -1;
  }
  kmer_data_base.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length,
                      _signature_len, _min_count, _max_count, _total_kmers);
  CKmerAPI kmer_object(_kmer_length);
  uint32 counter;
  std::string str;

  Histogram histogram(max_count);

  while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
    kmer_object.to_string(str);
    histogram.add_value(counter);
  }
  histogram.write_to_file(outputfile);

  // smooth the histogram
  histogram.smooth_histogram();
  // find peaks
  vector<size_t> peak_ids;
  vector<size_t> peak_values;
  histogram.find_peaks(peak_ids, peak_values);

  // identify the largest and second largest (if it exists)
  if (peak_ids.size() == 0) {
    cerr << "computeHistogram: no peak found in kmer-count histogram." << endl;
    return -1;
  }
  size_t kmer_coverage_estimate = -1;
  if (peak_ids.size() < 2) {
    cerr << "Histogram peak: " << peak_ids[0] << " (" << peak_values[0] << ")"
         << endl;
    kmer_coverage_estimate = peak_ids[0];
  } else {
    size_t largest, second, largest_id, second_id;
    if (peak_values[0] < peak_values[1]) {
      largest = peak_values[1];
      largest_id = peak_ids[1];
      second = peak_values[0];
      second_id = peak_ids[0];
    } else {
      largest = peak_values[0];
      largest_id = peak_ids[0];
      second = peak_values[1];
      second_id = peak_ids[1];
    }
    for (size_t i = 0; i < peak_values.size(); ++i) {
      if (peak_values[i] > largest) {
        second = largest;
        second_id = largest_id;
        largest = peak_values[i];
      } else if ((peak_values[i] > second) && (peak_values[i] != largest)) {
        second = peak_values[i];
        second_id = peak_ids[i];
      }
    }
    cerr << "Histogram peaks: " << largest_id << " (" << largest << "), "
         << second_id << " (" << second << ")" << endl;
    if (largest_peak) {
      kmer_coverage_estimate = largest_id;
    } else {
      kmer_coverage_estimate = second_id;
    }
  }

  // add expected abundance counts to end of hist file
  ofstream histofile;
  histofile.open(outputfile, ios::app);
  if (!histofile.good()) {
    stringstream ss;
    ss << "computeHistogram: File " << outputfile
       << " cannot be created. Note that the filename must not contain non-existing directories."
       << endl;
    return -1;
  }

  histofile << "parameters\t" << kmer_coverage_estimate / 2.0 << '\t'
            << kmer_coverage_estimate << endl;
  histofile.close();

  return 0;
}





