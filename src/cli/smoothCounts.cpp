#include "clean.hpp"

#include <ips4o.hpp>

#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/template_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/graph_extensions/node_weights.hpp"
#include "graph/graph_cleaning.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include <iostream>
#include "common/vector_map.hpp"
#include "seq_io/kmc_parser.hpp"
#include "common/utils/string_utils.hpp"
#include "histogram.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using mtg::seq_io::FastaWriter;
using mtg::seq_io::ExtendedFastaWriter;


int smoothCounts(Config *config) {
    assert(config);
    assert(config->outfbase.size());
    cerr<<"Starting SMoothing counts"<<endl;
   // logger->trace();
    const auto &files = config->fnames;
    uint32_t max_count = 1000;
    assert(files.size() == 1);

    string fasta_file=files[0];

    u_int32_t kSize=config->k;

    ExtendedFastaWriter<uint32_t> writer(config->outfbase,
                                         "kmer_counts",
                                         kSize );
    vector<uint32_t> counts;
    counts.reserve(1000000);
    seq_io::read_extended_fasta_file_critical<uint32_t>(fasta_file, "kmer_counts",
                                                [&](size_t k, const seq_io::kseq_t *read_stream, const uint32_t *kmer_counts) {
                                                    if (k != config->k) {
                                                        logger->error("File {} contains counts for k-mers of "
                                                                                   "length {} but graph is constructed with k={}",
                                                                      fasta_file, k, config->k);
                                                        exit(1);
                                                    }
                                                    assert(read_stream->seq.l >= k && "sequences can't be shorter than k-mers");
                                                    const uint32_t number_of_kmers = read_stream->seq.l - k + 1;


                                                    counts.resize(number_of_kmers);
                                                    uint32_t sum_kmer_count = 0;
                                                    for(unsigned i = 0 ;i < number_of_kmers; i++)
                                                    {
                                                        sum_kmer_count +=  *(kmer_counts+i);
                                                    }
                                                    const uint32_t avg_kmer_count = sum_kmer_count/ number_of_kmers;

                                                    for(unsigned i = 0 ;i < number_of_kmers; i++)
                                                    {
                                                        counts[i]= avg_kmer_count;
                                                    }

                                                    writer.write(read_stream->seq.s,counts);
                                                }
    );
    logger->trace("Smoothing finished");
    cerr<<"Smoothing Finished"<<endl;
    return 0;
}

} // namespace cli
} // namespace mtg
