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


int addCount(Config *config) {
    assert(config);
    assert(config->outfbase.size());
    logger->trace("Starting Add Count");
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
    Histogram histogram(max_count);
    //read the fasta
    seq_io::read_fasta_file_critical(fasta_file,
                                     [&](seq_io::kseq_t *stream) {
                                         string seq= stream->seq.s;
                                         counts.resize(seq.size()-kSize+1);
                                         string header= stream->comment.s;

                                         auto tokens=utils::split_string(header," ",true);
                                         unsigned i=0;
                                         for(;i < tokens.size();i++)
                                         {
                                            // std::cout<<t<<std::endl;
                                             if(utils::starts_with(tokens[i],"ab:Z:"))
                                             {

                                                // std::cout<<stream->name.s<<" : ";
                                                 break;
                                             }
                                         }

                                         for(unsigned j=0;j<seq.size()-kSize+1;j++)
                                         {
                                             if(j==0)
                                             {
                                                 counts[j]=std::atoi(tokens[i].substr(5,tokens[i].size()-5).c_str());
                                                 i++;
                                             }
                                             else
                                             {
                                                 counts[j] = std::atoi(tokens[i++].c_str());
                                             }
                                             histogram.add_value(counts[j]);
                                            // std::cout<<counts[j]<<" ";

                                         }
                                         //std::cout<<std::endl;
                                         writer.write(seq,counts);
                                     },
                                     false);

    //write the sequencues with extended
    string outputfile=config->outfbase+".hist";
    bool largest_peak = true;
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

    config->min_count = std::max(1u, config->min_count);

    if (!config->to_fasta) {
        logger->error("Clean graph can be serialized only in form "
                      "of its contigs or unitigs. Add flag --to-fasta");
        exit(1);
    }

    Timer timer;
    logger->trace("Loading graph...");

    auto graph = load_critical_dbg(files.at(0));
    // try loading k-mer counts
    auto node_weights = graph->load_extension<graph::NodeWeights>(files.at(0));

    if (node_weights) {
        if (auto *dbg_succ = dynamic_cast<graph::DBGSuccinct*>(graph.get()))
            dbg_succ->reset_mask();

        if (!node_weights->is_compatible(*graph)) {
            logger->error("Loaded k-mer counts not compatible with graph '{}'", files.at(0));
            exit(1);
        }
    }

    if (graph->get_mode() == graph::DeBruijnGraph::PRIMARY) {
        logger->warn("Cleaning primary graphs is not recommended."
                     " Consider building a canonical graph for cleaning instead.");
    }

    if (config->min_count > 1
            || config->max_count < std::numeric_limits<unsigned int>::max()
            || config->min_unitig_median_kmer_abundance != 1
            || config->count_slice_quantiles[0] != 0
            || config->count_slice_quantiles[1] != 1) {
        if (!node_weights) {
            logger->error("Cannot load k-mer counts from file '{}'", files.at(0));
            exit(1);
        }

        if (config->min_unitig_median_kmer_abundance == 0) {
            // skip zero k-mer counts for dummy k-mers in DBGSuccinct
            const auto _graph = dynamic_cast<graph::DBGSuccinct*>(graph.get())
                    ? std::make_shared<graph::MaskedDeBruijnGraph>(graph,
                        [&](auto i) { return (*node_weights)[i] > 0; }, true)
                    : graph;

            uint64_t cutoff
                = estimate_min_kmer_abundance(*_graph, *node_weights,
                                              config->num_singleton_kmers);

            if (cutoff != static_cast<uint64_t>(-1)) {
                config->min_unitig_median_kmer_abundance = cutoff;
            } else {
                if (config->fallback_abundance_cutoff == -1) {
                    logger->error("Cannot estimate expected minimum k-mer abundance"
                                  " and fallback is disabled (--fallback -1). Terminating.");
                    std::exit(129);
                }
                logger->warn("Cannot estimate expected minimum k-mer abundance."
                             " Using fallback value: {}", config->fallback_abundance_cutoff);
                config->min_unitig_median_kmer_abundance = config->fallback_abundance_cutoff;
            }
        }

        if (config->min_count > 1
                || config->max_count < std::numeric_limits<unsigned int>::max()) {
            const auto &weights = *graph->get_extension<graph::NodeWeights>();

            graph = std::make_shared<graph::MaskedDeBruijnGraph>(graph,
                [&](auto i) { return weights[i] >= config->min_count
                                    && weights[i] <= config->max_count; },
                true,
                graph->get_mode()
            );
            graph->add_extension(node_weights);

            assert(node_weights->is_compatible(*graph));
        }
    }

    logger->trace("Graph loaded in {} sec", timer.elapsed());

    if (dynamic_cast<const graph::MaskedDeBruijnGraph *>(graph.get())) {
        logger->trace("Extracting sequences from subgraph...");
    } else {
        logger->trace("Extracting sequences from graph...");
    }

    timer.reset();

    auto call_clean_contigs = [&](auto callback, size_t num_threads) {
        if (config->min_unitig_median_kmer_abundance != 1) {
            assert(node_weights);
            if (!node_weights->is_compatible(*graph)) {
                logger->error("k-mer counts are not compatible with the subgraph");
                exit(1);
            }

            logger->info("Threshold for median k-mer abundance in unitigs: {}",
                         config->min_unitig_median_kmer_abundance);

            graph->call_unitigs([&](const std::string &unitig, const auto &path) {
                if (!is_unreliable_unitig(path,
                                          *node_weights,
                                          config->min_unitig_median_kmer_abundance))
                    callback(unitig, path);
            }, num_threads, config->min_tip_size, graph->get_mode() == graph::DeBruijnGraph::CANONICAL);

        } else if (config->unitigs || config->min_tip_size > 1 || config->smoothing_window > 1) {
            graph->call_unitigs(callback,
                                num_threads,
                                config->min_tip_size,
                                graph->get_mode() == graph::DeBruijnGraph::CANONICAL);

        } else {
            graph->call_sequences(callback, num_threads,
                                  graph->get_mode() == graph::DeBruijnGraph::CANONICAL);
        }
    };

    auto dump_contigs_to_fasta = [&](const std::string &outfbase, auto call_contigs) {
        std::mutex seq_mutex;
        if (node_weights) {
            if (!node_weights->is_compatible(*graph)) {
                logger->error("k-mer counts are not compatible with the subgraph");
                exit(1);
            }

            ExtendedFastaWriter<uint32_t> writer(outfbase,
                                                 "kmer_counts",
                                                 graph->get_k(),
                                                 config->header,
                                                 config->enumerate_out_sequences,
                                                 get_num_threads() > 1);
            call_contigs([&](const std::string &contig, const auto &path) {
                std::vector<uint32_t> kmer_counts;
                kmer_counts.reserve(path.size());
                for (auto node : path) {
                    kmer_counts.push_back((*node_weights)[node]);
                }
                if(config->log_counts)
                {
                    for(unsigned i =0;i< kmer_counts.size(); i++) {
                        double tmp = log2(kmer_counts[i]);
                        double base = int(tmp);
                        if (tmp - base > 0.5)
                            base++;
                        kmer_counts[i] = base;
                    }
                }
                // smooth k-mer counts in the unitig
                utils::smooth_vector(config->smoothing_window, &kmer_counts);

                std::lock_guard<std::mutex> lock(seq_mutex);
                writer.write(contig, kmer_counts);
            }, get_num_threads());

        } else {
            FastaWriter writer(outfbase, config->header,
                               config->enumerate_out_sequences,
                               get_num_threads() > 1);

            call_contigs([&](const std::string &contig, const auto &) {
                std::lock_guard<std::mutex> lock(seq_mutex);
                writer.write(contig);
            }, get_num_threads());
        }
    };

    assert(config->count_slice_quantiles.size() >= 2);

    if (config->count_slice_quantiles[0] == 0
            && config->count_slice_quantiles[1] == 1) {
        dump_contigs_to_fasta(config->outfbase, call_clean_contigs);

    } else {
        if (!node_weights) {
            logger->error("Need k-mer counts for binning k-mers by abundance");
            exit(1);
        }
        assert(node_weights->is_compatible(*graph));

        auto &weights = node_weights->get_data();

        assert(graph->max_index() + 1 == weights.size());

        // compute clean count histogram
        tsl::hopscotch_map<uint64_t, uint64_t> count_hist;

        if (config->min_unitig_median_kmer_abundance != 1 || config->min_tip_size > 1) {
            // cleaning required
            sdsl::bit_vector removed_nodes(weights.size(), 1);

            call_clean_contigs([&](const std::string&, const auto &path) {
                for (auto i : path) {
                    assert(weights[i]);
                    count_hist[weights[i]]++;
                    removed_nodes[i] = 0;
                }
            }, get_num_threads());

            call_ones(removed_nodes, [&weights](auto i) { weights[i] = 0; });

        } else if (auto dbg_succ = std::dynamic_pointer_cast<graph::DBGSuccinct>(graph)) {
            // use entire graph without dummy BOSS edges
            graph->call_nodes([&](auto i) {
                if (uint64_t count = weights[i])
                    count_hist[count]++;
            });
        } else {
            // use entire graph
            graph->call_nodes([&](auto i) {
                assert(weights[i]);
                count_hist[weights[i]]++;
            });
        }
        // must not have any zero weights
        assert(!count_hist.count(0));

        std::vector<std::pair<uint64_t, uint64_t>> count_hist_v(count_hist.begin(),
                                                                count_hist.end());
        count_hist.clear();

        ips4o::parallel::sort(count_hist_v.begin(), count_hist_v.end(),
                              utils::LessFirst(), get_num_threads());

        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 1; i < config->count_slice_quantiles.size(); ++i) {
            // extract sequences for k-mer counts bin |i|
            assert(config->count_slice_quantiles[i - 1] < config->count_slice_quantiles[i]);

            auto filebase = utils::remove_suffix(config->outfbase, ".gz", ".fasta")
                    + fmt::format(".{}.{}", config->count_slice_quantiles[i - 1],
                                            config->count_slice_quantiles[i]);

            if (!count_hist_v.size()) {
                dump_contigs_to_fasta(filebase, [](auto, auto) {});
                continue;
            }

            uint64_t min_count = config->count_slice_quantiles[i - 1] > 0
                ? utils::get_quantile(count_hist_v, config->count_slice_quantiles[i - 1])
                : 1;
            uint64_t max_count = config->count_slice_quantiles[i] < 1
                ? utils::get_quantile(count_hist_v, config->count_slice_quantiles[i])
                : std::numeric_limits<uint64_t>::max();

            logger->info("k-mer count thresholds:\n"
                         "min (including): {}\n"
                         "max (excluding): {}", min_count, max_count);

            assert(node_weights->is_compatible(*graph));

            graph::MaskedDeBruijnGraph graph_slice(graph,
                [&](auto i) { return weights[i] >= min_count && weights[i] < max_count; },
                true,
                graph->get_mode()
            );

            // the outer for loop is parallelized, so don't start another thread
            // pool here
            dump_contigs_to_fasta(filebase, [&](auto dump_sequence, auto /* num_threads */) {
                graph_slice.call_sequences(dump_sequence, 1,
                                           graph_slice.get_mode() == graph::DeBruijnGraph::CANONICAL);
            });
        }
    }

    logger->trace("Graph cleaning finished in {} sec", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg
