
import pandas as pd

# expanding vars; used in input function
num_pod_simulations = int(float(config["pod_number_of_simulations_per"]))

rule calc_more_stats:
    output:
        df = temp("results/more_stats/noag/pod_{pod_no}.tsigma_{tsigma}.gzip")
    input:
        tss = sumstat_on_pyslim_input_func
    log:
        std = "results/log/more_stats/noag/pod_{pod_no}.tsigma_{tsigma}.log"
    params:
        recombrate = float(config["recombination_rate"]),
        pop_size_outcrossing = int(float(set(config["pod_population_size_recent"]).pop())),
        nsam = int(float(config["sample_size_simulation"])),
    threads: 1
    run:
        # calculate the vector of stats
        stats_collection = [more_stats(pyslim.load(ts), params, log)
            for ts in input.tss]

        # get average through all the different loci
        df = pd.concat(stats_collection).mean(axis=0)

        df.to_pickle(output.df, compression="gzip")


rule calc_more_stats_fpop:
    output:
        df = temp(
            "results/more_stats_fpop/noag/pod_{pod_no}.tsigma_{tsigma}.gzip")
    input:
        tss = sumstat_on_pyslim_input_func
    log:
        std = "results/log/more_stats_fpop/noag/pod_{pod_no}.tsigma_{tsigma}.log"
    params:
        recombrate = float(config["recombination_rate"]),
        pop_size_outcrossing = int(float(set(
            config["pod_population_size_recent"]).pop())),
        nsam = int(float(config["sample_size_simulation"])),
    threads: 1
    run:
        # calculate the vector of stats
        stats_collection = [more_stats_fpop(pyslim.load(ts), params, log)
            for ts in input.tss]

        # get average through all the different loci
        df = pd.concat(stats_collection).mean(axis=0)

        df.to_pickle(output.df, compression="gzip")


rule collect_more_stats:
    output:
        df = "results/more_stats/tsigma_{tsigma}.gzip"
    input:
        expand("results/more_stats/noag/pod_{pod_no}.tsigma_{{tsigma}}.gzip",
            pod_no=range(num_pod_simulations))
    log: "results/log/more_stats/tsigma_{tsigma}.log"
    params:
    threads: 1
    run:
        df = pd.concat([pd.read_pickle(infile, compression="gzip")
            for infile in input], axis=1).transpose()

        df.to_pickle(output.df, compression="gzip")


rule collect_more_stats_fpop:
    output:
        df = "results/more_stats_fpop/tsigma_{tsigma}.gzip"
    input:
        expand("results/more_stats_fpop/noag/pod_{pod_no}.tsigma_{{tsigma}}.gzip",
            pod_no=range(num_pod_simulations))
    log: "results/log/more_stats_fpop/tsigma_{tsigma}.log"
    params:
    threads: 1
    run:
        df = pd.concat([pd.read_pickle(infile, compression="gzip")
            for infile in input], axis=1).transpose()

        df.to_pickle(output.df, compression="gzip")


rule collected_stats_to_csv:
    output:
        csv = "results/more_stats/more_stats.collected.csv"
    input:
        expand("results/more_stats/tsigma_{tsigma}.gzip",
            tsigma=grep_all_tsigma())
    log: "results/log/more_stats/more_stats.collected.csv"
    params:
    threads: 1
    run:
        listed_df = []
        for infile in input:
            df = pd.read_pickle(infile, compression="gzip")
            df["tsigma"] = int(infile.split("tsigma_")[1].split(".")[0])
            listed_df.append(df)
            del df

        df = pd.concat(listed_df, axis=0)
        df.to_csv(output.csv, index=False)


rule collected_stats_to_csv_fpop:
    output:
        csv = "results/more_stats_fpop/more_stats.collected.csv"
    input:
        expand("results/more_stats_fpop/tsigma_{tsigma}.gzip",
            tsigma=grep_all_tsigma())
    log: "results/log/more_stats_fpop/more_stats.collected.csv"
    params:
    threads: 1
    run:
        listed_df = []
        for infile in input:
            df = pd.read_pickle(infile, compression="gzip")
            df["tsigma"] = int(infile.split("tsigma_")[1].split(".")[0])
            listed_df.append(df)
            del df

        df = pd.concat(listed_df, axis=0)
        df.to_csv(output.csv, index=False)


