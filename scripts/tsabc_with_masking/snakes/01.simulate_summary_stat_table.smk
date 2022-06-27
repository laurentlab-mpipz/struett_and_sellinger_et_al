
"""
snake to generate the abc table including parameters and summary statistics
"""

import msprime
import multiprocessing

rule draw_prior_param_and_simulate_summary_stat_table:
    output: "data/simulated/params_and_sumstats_{ident}_.gzip"
    input:
        rand_table = "data/rands/rands.gzip",
        mask_file = "resources/exons_for_slim.txt"
    params: configuration = config, verbose = False
    threads: config["nthreads_per_batch"]
    resources:
        mem_mb = 4000 * config["nthreads_per_batch"]
    run:
        import warnings

        # restrict threads of this job
        os.environ['MKL_NUM_THREADS']='1' 
        os.environ['OPENBLAS_NUM_THREADS']='1'
        os.environ['OMP_NUM_THREADS']='1'


        os.environ["NUM_INTER_THREADS"]="1"
        os.environ["NUM_INTRA_THREADS"]="1"

        os.environ["XLA_FLAGS"] = ("--xla_cpu_multi_thread_eigen=false "
                                   "intra_op_parallelism_threads=1")

        global main_simulator

        def main_simulator(identifier_of_simulation):
            vp = print_utils(params.verbose)
            ident = int(identifier_of_simulation)

            # read the seeds table
            df = pd.read_pickle(input.rand_table)
            df_seeds = df.loc[df['id'] == ident][['rand_prior',
                "rand_simulation"]]
            rand_prior = df_seeds.iloc[0]['rand_prior']
            rand_simulation = df_seeds.iloc[0]['rand_simulation']
            del df, df_seeds

            # create tsabc instance
            tsabc = tsabc2()
            #tsabc.test()

            # draw prior
            np.random.seed(rand_prior)
            population_size_recent = tsabc.draw_parameter_from_prior_distribution(
                params.configuration['population_size_recent'])
            population_size_ancient = population_size_recent
            population_size_change_time = 10
            selfing_rate_recent = tsabc.draw_parameter_from_prior_distribution(
                params.configuration['selfing_rate_recent'])
            selfing_rate_ancient = tsabc.draw_parameter_from_prior_distribution(
                params.configuration['selfing_rate_ancient'])
            selfing_rate_change_time = tsabc.draw_parameter_from_prior_distribution(
                params.configuration['selfing_rate_change_time'])
            vector_prior_parameter = np.array([population_size_recent,
                selfing_rate_recent, selfing_rate_ancient,selfing_rate_change_time])
            vp.vprint(f"parameters\t{' '.join(vector_prior_parameter.astype(str))}")

            # simulate
            recombination_rate = float(params.configuration['recombination_rate'])
            mutation_rate = float(params.configuration['mutation_rate'])
            chromosome_length = int(float(params.configuration['chromosome_length']))
            sample_size = int(float(params.configuration['sample_size_simulation']))
            generations_dtwf = int(float(params.configuration['generations_dtwf']))

            df_list = []
            for independent_region_no in range(params.configuration["num_independent_regions"]):
                print("region: ", independent_region_no+1, "of", params.configuration["num_independent_regions"])

                rand_simulation = np.random.randint(*params.configuration["randint_range_for_simulations"])
                ts_unmutated = tsabc.run_transition_to_selfing_model(
                    population_size_recent, population_size_change_time, population_size_ancient,
                    selfing_rate_recent, selfing_rate_change_time, selfing_rate_ancient,
                    recombination_rate, sample_size, chromosome_length, generations_dtwf,
                    rand_simulation)

                ts = msprime.mutate(ts_unmutated, rate = mutation_rate, keep = False)
                del ts_unmutated

                before_sites = ts.num_sites
                mask_all = pd.read_csv(input.mask_file, sep="\t", header=None)
                mask_chr = mask_all.loc[mask_all[0]==independent_region_no+1]
                region_array = np.array(list(zip(mask_chr[1], mask_chr[2])))
                # check for disjointness
                change_last_pos = False
                range_l = 0
                range_r = ts.sequence_length
                last_l, last_r = 0, 0
                filtered_regions = []
                for i, (l, r) in enumerate(region_array, start=1):
                    ISOK = True
                    if l > r:
                        warnings.warn("Bad interval: right <= left")
                        ISOK = False
                    elif l < last_r:
                        warnings.warn(f"Intervals must be disjoint. {last_r, l, r}")
                        ISOK = False
                    elif l < range_l or r > range_r:
                        if i == len(region_array):
                            change_last_pos = True
                        else:
                            warnings.warn(f"Intervals must be within {range_l} and {range_r}: {l, r}")
                    last_r = r

                    if ISOK:
                        filtered_regions.append([l, r])

                # another disjoint test
                bad_regions = []
                i_, j_ = -2, -1
                for ix, (i, j) in enumerate(filtered_regions):
                    if not i > j_:
                        bad_regions.append(ix)
                    i_, j_ = i, j

                if len(bad_regions):
                    for index in sorted(bad_regions, reverse=True):
                        del filtered_regions[index]

                plus_filtered_regions = []
                minus_filtered_regions = []
                odd_regions = []
                # remove genes on reverse strand
                for i, j in filtered_regions:
                    if i < j:
                        plus_filtered_regions.append([i, j])
                    elif i > j:
                        minus_filtered_regions.append([i, j])
                    else:
                        odd_regions.append([i, j])

                ts = ts.delete_intervals(plus_filtered_regions)

                print("before|after:", before_sites, ts.num_sites)

                # get sfs
                sfs_stat_mode, sfs_stat_polarity = params.configuration['sfs_stat_properties']
                polarised = lambda x: x == "unfolded"
                vector_sfs = ts.allele_frequency_spectrum(sample_sets=None,
                    windows=None, mode=sfs_stat_mode, span_normalise=False,
                    polarised=polarised(sfs_stat_polarity))
                del polarised
                vp.vprint(f"sfs\t{' '.join(vector_sfs.astype(int).astype(str))}")

                # get binned ld
                num_considered_mutatations = int(float(
                    params.configuration['ld_num_of_considered_mutations']))
                breaks = np.array([float(x) for x in params.configuration['ld_binning_breaks']])
                breaks = sorted(breaks[breaks <= chromosome_length])
                vector_ld = tsabc.get_binned_ld_from_ts(ts, num_considered_mutatations,
                    chromosome_length, breaks, np.random.randint(1, 4294967296))
                del breaks
                vp.vprint(f"ld\t{' '.join(vector_ld.astype(str))}")

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
                vp.vprint(f"vector_tl_true\t{' '.join(vector_tl_true.astype(str))}")

                # get tm_true
                tmrcas = [np.array(tl_true_of_pair)[:,1] for tl_true_of_pair in tl_table]
                dtmrcas = tsabc.discretize_list_tmrca_lists(tmrcas, breaks)
                tms = tsabc.transition_matrices_from_tiscretized_tmrcas(dtmrcas, len(breaks))
                vector_tm_true = np.mean(tms, axis=0).flatten('C')
                del tmrcas, dtmrcas, tms, breaks
                vp.vprint(f"vector_tm_true\t{' '.join(vector_tm_true.astype(str))}")
                vp.vprint(f"zero probability content: {round(100*(vector_tm_true==0).sum()/len(vector_tm_true), 2)}%")

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
                vp.vprint(f"vector_tm_win\t{' '.join(vector_tm_win.astype(str))}")
                vp.vprint(f"zero probability content: {round(100*(vector_tm_win==0).sum()/len(vector_tm_win), 2)}%")

                # execution time
                exec_time = vp.time_measure()

                # concatenate and save to pickle
                df_ident = tsabc.vector_to_DataFrame([ident], 'ident', ident)
                df_exec = tsabc.vector_to_DataFrame([exec_time], 'exec_time', ident)
                df_rand = tsabc.vector_to_DataFrame([rand_prior, rand_simulation],
                    'rand', ident)
                df_prior_parameter = tsabc.vector_to_DataFrame(vector_prior_parameter,
                    'param', ident)
                df_sfs = tsabc.vector_to_DataFrame(vector_sfs, 'sfs', ident)
                df_ld = tsabc.vector_to_DataFrame(vector_ld, 'ld', ident)
                df_tl_true = tsabc.vector_to_DataFrame(vector_tl_true, 'tl_true', ident)
                df_tm_true = tsabc.vector_to_DataFrame(vector_tm_true, 'tm_true', ident)
                df_tm_win = tsabc.vector_to_DataFrame(vector_tm_win, 'tm_win', ident)
                df = pd.concat([df_ident, df_exec, df_rand, df_prior_parameter, df_sfs,
                    df_ld, df_tl_true, df_tm_true, df_tm_win], axis=1)
                del df_ident, df_exec, df_rand, df_prior_parameter, df_sfs, df_ld, df_tl_true, df_tm_true, df_tm_win

                print("\n____________")
                print(f"num_mutations:\t{ts.num_mutations}")
                print(f"num_sites:\t{ts.num_sites}")
                print(f"sum_sfs:\t{vector_sfs.sum()}")
                print("exec_time", "{:.2f} sec".format(exec_time), sep="\t")

                df_list.append(df)
                del df

            # make average of the sumstats
            df_average = pd.DataFrame(pd.concat(df_list).mean()).transpose()

            return df_average

        # conditional multithreading
        if threads > 1:
            with multiprocessing.Pool(threads) as p:
                df = pd.concat(p.map(main_simulator, batch_dict_simulations[float(wildcards.ident)]))
        else:
            df = pd.concat([main_simulator(i) for i in batch_dict_simulations[float(wildcards.ident)]])


        # save as binary
        df.to_pickle(str(output))

rule aggregate_simulated_data:
    output: "data/table/params_and_sumstats.table.gzip"
    input:
        simstats = expand("data/simulated/params_and_sumstats_{ident}_.gzip",
            ident=batch_dict_simulations.keys())
    threads: 1
    run:
        df = pd.concat([pd.read_pickle(infile) for infile in input.simstats])
        df.to_pickle(str(output))
