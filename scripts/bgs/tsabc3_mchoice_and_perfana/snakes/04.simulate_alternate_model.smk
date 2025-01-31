
"""
snake to generate the abc table including parameters and summary statistics
"""

import msprime
import os

rule draw_prior_param_and_simulate_summary_stat_table_alternate_model:
    output: "data/simulated_alt_model/params_and_sumstats_{ident}_.gzip"
    input: rand_table = "data/rands/rands_alternate.gzip"
    params: configuration = config, verbose = False
    threads: 1
    resources:
        mem_mb=500
    run:
        # restrict threads of this job
        os.environ['MKL_NUM_THREADS']='1' 
        os.environ['OPENBLAS_NUM_THREADS']='1'
        os.environ['OMP_NUM_THREADS']='1'


        os.environ["NUM_INTER_THREADS"]="1"
        os.environ["NUM_INTRA_THREADS"]="1"

        os.environ["XLA_FLAGS"] = ("--xla_cpu_multi_thread_eigen=false "
                                   "intra_op_parallelism_threads=1")

        # simulator
        def main(identifier_of_simulation):
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
                params.configuration['alt_model_population_size_recent'])
            population_size_ancient = tsabc.draw_parameter_from_prior_distribution(
                params.configuration['alt_model_population_size_ancient'])
            population_size_change_time = tsabc.draw_parameter_from_prior_distribution(
                params.configuration['alt_model_population_size_change_time'])
            selfing_rate_recent = tsabc.draw_parameter_from_prior_distribution(
                params.configuration['alt_model_selfing_rate'])
            selfing_rate_ancient = selfing_rate_recent
            selfing_rate_change_time = 10
            vector_prior_parameter = np.array([population_size_recent,
                population_size_ancient,
                population_size_change_time,
                selfing_rate_recent])
            vp.vprint(f"parameters\t{' '.join(vector_prior_parameter.astype(str))}")

            # simulate
            recombination_rate = float(params.configuration['recombination_rate'])
            mutation_rate = float(params.configuration['mutation_rate'])
            chromosome_length = int(float(params.configuration['chromosome_length']))
            sample_size = int(float(params.configuration['sample_size_simulation']))
            generations_dtwf = int(float(params.configuration['generations_dtwf']))

            df_list = []
            for independent_region_no in range(params.configuration["num_independent_regions"]):
                print("region: ", independent_region_no, "of", params.configuration["num_independent_regions"])

                rand_simulation = np.random.randint(*params.configuration["randint_range_for_simulations"])

                ts_unmutated = tsabc.run_transition_to_selfing_model(
                    population_size_recent, population_size_change_time,
                    population_size_ancient, selfing_rate_recent,
                    selfing_rate_change_time, selfing_rate_ancient, recombination_rate,
                    sample_size, chromosome_length, generations_dtwf, rand_simulation)

                ts = msprime.mutate(ts_unmutated, rate = mutation_rate, keep = False)
                del ts_unmutated

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
                tmrcas = [np.array(tl_true_of_pair)[:,0] for tl_true_of_pair in tl_table]
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

                print("\n____________\nthis has been great, thank you")
                print(f"num_mutations:\t{ts.num_mutations}")
                print(f"num_sites:\t{ts.num_sites}")
                print(f"sum_sfs:\t{vector_sfs.sum()}")
                print("exec_time", "{:.2f} sec".format(exec_time), sep="\t")

                df_list.append(df)
                del df

            # make average of the sumstats
            df_average = pd.DataFrame(pd.concat(df_list).mean()).transpose()

            return df_average

        # execution
        df = pd.concat([main(my_ident) for my_ident in batch_dict_simulations_alt[float(wildcards.ident)]])

        # save as binary
        df.to_pickle(str(output))

rule aggregate_simulated_data_alternate_model:
    output: "data/table/params_and_sumstats.alternate_model.table.gzip"
    input:
        simstats = expand("data/simulated_alt_model/params_and_sumstats_{ident}_.gzip",
            ident=batch_dict_simulations_alt.keys())
    threads: 1
    run:
        df = pd.concat([pd.read_pickle(infile) for infile in input.simstats])
        df.to_pickle(str(output))
