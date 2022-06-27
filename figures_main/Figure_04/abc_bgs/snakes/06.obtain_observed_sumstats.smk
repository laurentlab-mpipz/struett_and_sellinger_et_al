
"""
snake to obtain the summary statistics from an infered tree sequence,
 e.g. the CEU cluster of A. thaliana
 the sequence length has to be of the same size as the simulations
"""

import tskit
import tqdm

rule obtain_observed_sumstats:
    output: "data/observed/stats.{observed_tree_sequence}.table.gzip"
    input: config["observed_tree_sequence"]
    params: configuration = config, verbose = False
    threads: 1
    run:
        def get_sample_set_from_tree_sequence(tree_sequence, nsam):
            """obtain a sample set of given size"""
            if tree_sequence.num_samples <= nsam:
                print(f"sample size larger than samples in tree sequence")
                nsam = tree_sequence.num_samples
                print(f"set sample size to {nsam}\n\n")

            tables = tree_sequence.tables

            samples = np.where(tables.nodes.flags == tskit.NODE_IS_SAMPLE)[0]

            return sorted(np.random.choice(samples, size=nsam, replace=False))

        def main(tree_sequence):
            vp = print_utils(params.verbose)
            
            # create tsabc instance
            tsabc = tsabc2()

            # read the tree sequence
            if not isinstance(tree_sequence, tskit.trees.TreeSequence):
                tree_sequence = tskit.load(tree_sequence)

            # define samples to observe from
            sample_set_to_observe = get_sample_set_from_tree_sequence(tree_sequence,
                params.configuration["sample_size_observations"])

            # reduce to sample size
            ts = tree_sequence.simplify(samples=sample_set_to_observe)

            # get sfs
            sfs_stat_mode, sfs_stat_polarity = params.configuration['sfs_stat_properties']
            polarised = lambda x: x == "unfolded"
            vector_sfs = ts.allele_frequency_spectrum(sample_sets=None,
                windows=None, mode=sfs_stat_mode, span_normalise=False,
                polarised=polarised(sfs_stat_polarity))
            del polarised

            # params for binning
            chromosome_length = int(float(params.configuration['chromosome_length']))

            # get binned ld
            num_considered_mutatations = int(float(
                params.configuration['ld_num_of_considered_mutations']))
            breaks = np.array([float(x) for x in params.configuration['ld_binning_breaks']])
            breaks = sorted(breaks[breaks <= chromosome_length])
            vector_ld = tsabc.get_binned_ld_from_ts(ts, num_considered_mutatations,
                chromosome_length, breaks, np.random.randint(1, 4294967296))
            del breaks

            # get binned tl_true
            breaks = tsabc.get_breaks_for_binnning_tl_true(
                float(params.configuration['tl_true_binning_scale_factor']),
                np.mean([float(params.configuration['pod_population_size_recent'][0]),
                    float(params.configuration['pod_population_size_ancient'][0])]),
                int(float(params.configuration['tl_true_number_of_bins'])))
            tl_true, num_pairs, tl_table = tsabc.get_list_tmrcas_and_lengths(ts,
                int(float(params.configuration['tl_true_max_number_of_pairs'])))
            vector_tl_true = tsabc.get_binned_tmrcas_from_true_tl_array(tl_true,
                breaks, chromosome_length, num_pairs,
                float(params.configuration['tl_true_minimal_length_bin_size']))

            # get tm_true
            tmrcas = [np.array(tl_true_of_pair)[:,0] for tl_true_of_pair in tl_table]
            dtmrcas = tsabc.discretize_list_tmrca_lists(tmrcas, breaks)
            tms = tsabc.transition_matrices_from_tiscretized_tmrcas(dtmrcas, len(breaks))
            vector_tm_true = np.mean(tms, axis=0).flatten('C')
            del tmrcas, dtmrcas, tms, breaks

            # get tm_win
            winpis = tsabc.get_list_of_pairwise_diversities(ts,
                int(float(params.configuration['tl_true_max_number_of_pairs'])),
                params.configuration['tm_win_properties'][2],
                int(float(params.configuration['tm_win_properties'][0])), 2,
                int(float(params.configuration['tm_win_properties'][1])),
                np.random.randint(1, 4294967296) ).to_numpy()[:,2:].transpose()
            breaks = tsabc.get_breaks_for_binning_tm_win(
                float(params.configuration['tm_win_binning_scale_factor']),
                int(float(params.configuration['tm_win_properties'][0])),
                float(params.configuration['mutation_rate']),
                np.mean([float(params.configuration['pod_population_size_recent'][0]),
                    float(params.configuration['pod_population_size_ancient'][0])]),
                int(float(params.configuration['tm_win_properties'][3])) )
            dwinpis = np.array(tsabc.discretize_list_tmrca_lists(winpis, breaks) )
            tms = tsabc.transition_matrices_from_tiscretized_tmrcas(dwinpis, len(breaks))
            vector_tm_win = np.mean(tms, axis=0).flatten('C')
            del winpis, dwinpis, tms, breaks

            # execution time
            exec_time = vp.time_measure()

            # concatenate into dataframe
            df_which_strains = tsabc.vector_to_DataFrame(sample_set_to_observe, 'sample', None)
            df_exec = tsabc.vector_to_DataFrame([exec_time], 'exec_time', None)
            df_sfs = tsabc.vector_to_DataFrame(vector_sfs, 'sfs', None)
            df_ld = tsabc.vector_to_DataFrame(vector_ld, 'ld', None)
            df_tl_true = tsabc.vector_to_DataFrame(vector_tl_true, 'tl_true', None)
            df_tm_true = tsabc.vector_to_DataFrame(vector_tm_true, 'tm_true', None)
            df_tm_win = tsabc.vector_to_DataFrame(vector_tm_win, 'tm_win', None)
            df = pd.concat([df_which_strains, df_exec, df_sfs,
                df_ld, df_tl_true, df_tm_true, df_tm_win], axis=1)

            return df

        # load tree sequence
        ts = tskit.load(input[int(wildcards.observed_tree_sequence)])

        # sample n times
        sumstat_list = []
        error_counter = 0
        for i in tqdm.tqdm(range(params.configuration["nsample_observations"])):
            try:
                sumstat_list.append(main(ts))
            except:
                error_counter += 1
        print(f"successful sumstat calculations: {params.configuration['nsample_observations'] - error_counter}")
        
        df = pd.concat(sumstat_list)

        # save as binary
        df.to_pickle(str(output))




