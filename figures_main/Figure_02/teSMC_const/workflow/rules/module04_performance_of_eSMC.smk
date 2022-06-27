
import pandas as pd

localrules: mse_eSMC, aggregate_mse_eSMC, rmse_eSMC, aggregate_rmse_eSMC, calc_rmse_rb_eSMC

rule mse_eSMC:
    output:
        demgr = temp("results/mse/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.{inference_model}.pop.csv"),
        sigma = temp("results/mse/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.{inference_model}.sigma.csv")
    input:
        full_data = "results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.{inference_model}.RDS"
    group: "aggregate_mse"
    script:
        "../scripts/mse.R"

rule aggregate_mse_eSMC:
    output:
        table = "results/mse/table.csv"
    input:
        fn = expand("results/mse/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.{inference_model}.{param}.csv",
            param=["pop", "sigma"],
            demography=range(len(config["population_sizes_backward_in_time"])),
            replication=range(config["technical_replicates"]),
            sample_size=config["sample_size"],
            mutation_rate=config["mutation_rate"],
            recombination_rate=config["recombination_rate"],
            chromosome_length=config["chromosome_length"],
            inference_model=config["teSMC"]["inference_models"]
        )
    group: "aggregate_mse"
    run:
        def wildcards_from_filename(filename):
            demography, rep, nsam, mu, r, L, infmodel, param = re.findall("results/mse/demographic_scenario_(.+)/" +
                "rep_(.+)\.nsam_(.+)\.mu_(.+)\.r_(.+)\.length_(.+)\.(.+)\.(.+)\.csv", filename)[0]

            return {
                "demography" : [int(demography)],
                "rep" : [int(rep)],
                "nsam" : [int(nsam)],
                "mu" : [float(mu)],
                "r" : [float(r)],
                "L" : [int(float(L))],
                "infmodel" : [infmodel],
                "param" : [param]
            }

        def main():
            df_dict = {}
            for fn in tqdm.tqdm(input.fn):
                # extract wildcards from name
                wc = wildcards_from_filename(fn)
                wc["mse"] = [pd.read_csv(fn)["V1"][0]]
                wc["filename"] = [fn]

                for k, v in wc.items():
                    if k in df_dict:
                        df_dict[k].extend(v)
                    else:
                        df_dict[k] = v

            df = pd.DataFrame.from_dict(df_dict)
            df.to_csv(output.table)

        main()

rule rmse_eSMC:
    output:
        rmse = temp("results/rmse/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.csv")
    input:
        allinf = "results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.csv"
    # group: "aggregate_mse"
    threads: 1
    script:
        "../scripts/rmse.R"

rule aggregate_rmse_eSMC:
    output:
        table = "results/rmse/table_point_estimates.csv"
    input:
        fn = expand("results/rmse/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.csv",
            demography=range(len(config["population_sizes_backward_in_time"])),
            replication=range(config["technical_replicates"]),
            sample_size=config["sample_size"],
            mutation_rate=config["mutation_rate"],
            recombination_rate=config["recombination_rate"],
            chromosome_length=config["chromosome_length"]
        )
    group: "aggregate_mse"
    run:
        def main():
            df_list = []
            for fn in tqdm.tqdm(input.fn):
                df = pd.read_csv(fn, index_col=0)
                df["file"] = fn
                df_list.append(df)
                del df

            df = pd.concat(df_list)

            df.to_csv(output.table)

        main()

rule calc_rmse_rb_eSMC:
    output:
        table = "results/rmse/table.csv",
        plot = "results/figures/rb_rmse.pdf"
    input:
        pesti = "results/rmse/table_point_estimates.csv"
    group: "aggregate_mse"
    script:
        "../scripts/rmse_calc.R"

rule aggregate_lh:
    output:
        table = "results/lh/lhd.csv"
    input:
        fn = expand("results/esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.{inference_model}.log",
            demography=range(len(config["population_sizes_backward_in_time"])),
            replication=range(config["technical_replicates"]),
            sample_size=config["sample_size"],
            mutation_rate=config["mutation_rate"],
            recombination_rate=config["recombination_rate"],
            chromosome_length=config["chromosome_length"],
            inference_model=config["teSMC"]["inference_models"]
        )
    group: "lh"
    run:
        non_count = 0
        lh = []

        for infile in tqdm.tqdm(input.fn, " lh"):
            new_lh = None
            old_lh = None
            with open(infile, "r") as fin:
                for l in fin:
                    if "new likelihood" in l.lower():
                        new_lh = float(l.rstrip().split()[-1][0:-1])
                    elif "old likelihood" in l.lower():
                        old_lh = float(l.rstrip().split()[-1][0:-1])

            demographic_scenario = infile.split("/")[2].split("_")[2]
            wc = infile.split("/")[3].split(".")
            rep, nsam, mu, r, L, infmodel = (
                wc[0].split("_")[1], wc[1].split("_")[1], wc[2].split("_")[1],
                wc[3].split("_")[1], wc[4].split("_")[1], wc[5]
                )

            if new_lh is None:
                non_count += 1
                this_lh = old_lh
            else:
                this_lh = new_lh

            lh.append([this_lh, demographic_scenario, rep, nsam, mu, r, L,
                infmodel, infile])

        print(f"we found {non_count} failed optimizations")
        df = pd.DataFrame(lh, columns=["lh", "demographic_scenario", "rep",
            "nsam", "mu", "r", "L", "infmodel", "infile"])

        df.to_csv(output.table)

rule point_estimates_for_oneTrans:
    output: 
        csv = "results/esmc/point_estim_for_OneTrans.csv"
    input:
        rds = "results/esmc/table.rds"
    script:
        "../scripts/pointEstimsFor_OneTrans.R"














