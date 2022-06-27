
import numpy as np
import pandas as pd
import pyslim

localrules: collect_sumstats_to_pod, rename_collected

# expanding vars; used in input function
num_pod_simulations = int(float(config["pod_number_of_simulations_per"]))

rule sumstat_on_pyslim:
    output:
        df = "results/stats/pod_{pod_no}/tsigma_{tsigma}.gzip"
    input:
        sumstat_on_pyslim_input_func
    params:
        recombrate = float(config["recombination_rate"]),
        pop_size_outcrossing = int(float(set(config["pod_population_size_recent"]).pop())),
        nsam = int(float(config["sample_size_simulation"])),
        configuration = config,
        verbose = False
    log:
        std = "results/log/stats/pod_{pod_no}/tsigma_{tsigma}.log",
        err = "results/log/stats/pod_{pod_no}/tsigma_{tsigma}.err"
    threads: 1
    run:
        treeseqs = (pyslim.load(infile) for infile in input)  # load all treeseqs

        # recaptitate; only for the use of the summarizing function; no mutating
        recaps = [ts.recapitate(recombination_rate=params.recombrate,
            Ne=params.pop_size_outcrossing) for ts in treeseqs]

        tss_1hapInd = (sample_1perInd(recap) for recap in recaps)  # sample a single haplotype per individual

        # subsample
        tss_1perInd_nsam = [random_sample_from_treeseq(ts_1hapInd,
            sample_size=params.nsam) for ts_1hapInd in tss_1hapInd]

        # some tree stats
        with open(log.err, 'w') as loge:
            for tx, ts in enumerate(tss_1perInd_nsam, start=1):
                print("_"*80, file=loge)
                print(f"Tree no: {tx}/{len(tss_1perInd_nsam)}", file=loge)
                print(ts, file=loge)
                print("="*80, end="\n"*2, file=loge)

        # function to calculate the stats
        def calc_stat(ts, ident):
            vp = print_utils(params.verbose)

            # create tsabc instance
            tsabc = tsabc2()
            #tsabc.test()

            # value preparations
            chromosome_length = int(ts.sequence_length)

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
            dwinpis = tsabc.discretize_list_tmrca_lists(winpis, breaks)
            tms = tsabc.transition_matrices_from_tiscretized_tmrcas(dwinpis, len(breaks))
            vector_tm_win = np.mean(tms, axis=0).flatten('C')
            del winpis, dwinpis, tms, breaks
            vp.vprint(f"vector_tm_win\t{' '.join(vector_tm_win.astype(str))}")
            vp.vprint(f"zero probability content: {round(100*(vector_tm_win==0).sum()/len(vector_tm_win), 2)}%")

            # execution time
            exec_time = vp.time_measure()

            # pod_iloc_from tsigma
            pod_tsigma_s = np.array([int(i) for i in params.configuration["pod_selfing_rate_change_time"]])
            pod_iloc = np.where(pod_tsigma_s == int(float(wildcards.tsigma)))[0]

            # make prior vecotr
            vector_pior = np.array([
                int(float(params.configuration["pod_population_size_recent"][int(float(pod_iloc))])),
                float(params.configuration["pod_selfing_rate_recent"][int(float(pod_iloc))]),
                float(params.configuration["pod_selfing_rate_ancient"][int(float(pod_iloc))]),
                int(float(params.configuration["pod_selfing_rate_change_time"][int(float(pod_iloc))])),
            ])

            # concatenate and save to pickle
            df_params = tsabc.vector_to_DataFrame(vector_pior, 'param', ident)
            df_sfs = tsabc.vector_to_DataFrame(vector_sfs, 'sfs', ident)
            df_ld = tsabc.vector_to_DataFrame(vector_ld, 'ld', ident)
            df_tl_true = tsabc.vector_to_DataFrame(vector_tl_true, 'tl_true', ident)
            df_tm_true = tsabc.vector_to_DataFrame(vector_tm_true, 'tm_true', ident)
            df_tm_win = tsabc.vector_to_DataFrame(vector_tm_win, 'tm_win', ident)
            df = pd.concat([df_params, df_sfs, df_ld,
                df_tl_true, df_tm_true, df_tm_win], axis=1)
            del df_sfs, df_tl_true, df_tm_true, df_tm_win, df_ld

            if params.verbose:
                print("\n____________\nthis has been great, thank you")
                print(f"num_mutations:\t{ts.num_mutations}")
                print(f"num_sites:\t{ts.num_sites}")
                print(f"sum_sfs:\t{vector_sfs.sum()}")
                print("exec_time", "{:.2f} sec".format(exec_time), sep="\t")

            return df

        df = pd.concat([calc_stat(ts, ix) for ix, ts in 
            enumerate(tss_1perInd_nsam)])
        ps = pd.DataFrame(df.mean(axis=0)).transpose()

        ps.to_pickle(output.df, compression="gzip")

rule collect_sumstats_to_pod:
    output:
        df = temp("results/stats/prename/all_pod_indices.tsigma_{tsigma}.gzip")
    input: 
        expand("results/stats/pod_{pod_no}/tsigma_{{tsigma}}.gzip",
            pod_no=range(num_pod_simulations))
    params:
    log:
        std = "results/log/collect.tsigma_{tsigma}.log"
    threads: 1
    run:
        df = pd.concat([pd.read_pickle(f, compression="gzip")
            for f in input], axis=0)

        # get all 6 params separately
        param_0 = df.param_0.copy()
        param_1 = df.param_0.copy()
        param_2 = [10 for i in range(len(df))]
        param_3 = df.param_1.copy()
        param_4 = df.param_2.copy()
        param_5 = df.param_3.copy()

        df["param_0"] = param_0
        df["param_1"] = param_1
        df["param_2"] = param_2
        df["param_3"] = param_3
        df["param_4"] = param_4
        df["param_5"] = param_5

        # re-order the columns
        odered_colnames = []
        _ = [odered_colnames.append(i) for i in df.columns if "param" in i]
        _ = [odered_colnames.append(i) for i in df.columns if "sfs" in i]
        _ = [odered_colnames.append(i) for i in df.columns if "ld" in i]
        _ = [odered_colnames.append(i) for i in df.columns if "tl_true" in i]
        _ = [odered_colnames.append(i) for i in df.columns if "tm_true" in i]
        _ = [odered_colnames.append(i) for i in df.columns if "tm_win" in i]

        assert len(odered_colnames) <= len(df.columns), "weired colsort"
        if len(odered_colnames) < len(df.columns):
            # append the missing values to ordered_colnames
            _ = [odered_colnames.append(el) for el in df.columns
                if el not in odered_colnames]
        assert len(odered_colnames) == len(df.columns), "missing columns after sorting"

        df = df[odered_colnames]
        
        df.to_pickle(output.df, compression=None)

rule rename_collected:
    output:
        df = "results/stats/params_and_sumstats_pod_{pod_iloc}.table.gzip"
    input:
        df = lambda wc: f"results/stats/prename/all_pod_indices.tsigma_{tsigma_from_pod_iloc(wc)}.gzip"
    log: 
        err = "results/log/stats/params_and_sumstats_pod_{pod_iloc}.table.gzip"
    params:
    conda: "config/env.yml"
    threads: 1
    shell:
        r"""
        cp {input.df} {output.df} 2>{log.err} 
        """


