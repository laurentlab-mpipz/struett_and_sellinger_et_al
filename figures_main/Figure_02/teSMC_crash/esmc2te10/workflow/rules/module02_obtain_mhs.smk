
"""
output mhs from tree sequences;
output tmrca on seqlen (what Thibaut calls ARG on n=2) as input for eSMC
"""

localrules: write_mhs_from_tree_sequence


rule write_mhs_from_tree_sequence:
    output:
        mhs = "results/mhs/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.mhs"
    input:
        ts = expand("results/tree_sequences/demographic_scenario_{{demography}}/"
            "rep_{{replication}}.nsam_{{sample_size}}.mu_{{mutation_rate}}."
            "r_{{recombination_rate}}.length_{{chromosome_length}}.ts_{chr_number}",
            chr_number=range(config["chromosome_number"]))
    params:
        chr_identifier = "{replication}"
    group: "create_input_files_for_analysis"
    threads: 1
    priority: 50
    run:
        # create new file:
        with open(output.mhs, "w") as outf:
            pass

        # loop through independent chromosomes
        for chr_num, ints in enumerate(input.ts, start=1):
            ts = tskit.load(ints)

            previous_site_position = 1
            nmultiallelics = 0
            print("we consider haplotypes being A|T for 0|1")
            with open(output.mhs, "a") as outf:
                for variant in tqdm.tqdm(ts.variants(), total=ts.num_sites):
                    # num_called_sites since the last heterozygous site
                    num_called_sites = int(round(variant.position, 0)) - int(round(
                        previous_site_position, 0))
                    previous_site_position = int(round(variant.position, 0))

                    # do not allow for multiallelic sites
                    if not num_called_sites:
                        nmultiallelics += 1
                        continue

                    # haplotypes
                    haplotypes = "".join(["A" if i == 0 else "T" if i == 1 else "-"
                        for i in variant.genotypes])

                    print(chr_num, int(round(variant.position, 0)),
                        num_called_sites, haplotypes, sep="\t", end="\n", file=outf)
            print(f"we found {nmultiallelics} multiallelic sites")
            del ts

rule provide_esmc_input:
    output:
        ps_esmc = "results/ps_esmc/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.pars"
    input:
        ts = "results/tree_sequences/demographic_scenario_{demography}/"
            "rep_{replication}.nsam_{sample_size}.mu_{mutation_rate}."
            "r_{recombination_rate}.length_{chromosome_length}.ts"
    params:
        config = config
    group: "create_input_files_for_analysis"
    threads: 1
    run:
        ts = tskit.load(input.ts)
        N_0 = params.config["population_sizes_backward_in_time"][int(float(
            wildcards.demography))][0][0]

        print(f"pop size at time 0: {N_0}")

        with open(output.ps_esmc, "w") as f:
            f.write(f"num sites{ts.num_sites}\n")
            for tree in ts.trees():
                samples = [i for i in ts.samples()]
                assert len(samples) == int(float(wildcards.sample_size)), "sample size is not correct"
                u, v = samples[0:2]
                f.write(f"{str(tree.span)}time{str(tree.tmrca(u, v)/(4*float(N_0)))}\n")
