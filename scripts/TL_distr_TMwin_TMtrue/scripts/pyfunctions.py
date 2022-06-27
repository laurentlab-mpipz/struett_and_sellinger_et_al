import sys
import numpy as np
import msprime
import tskit
import itertools
import pandas as pd
from heapq import nlargest
import time

class print_utils:
    def __init__(self, verbose: bool):
        """print utilities"""
        self.verbose = verbose
        self.start = time.time()

    def vprint(self, content: str):
        if self.verbose:
            print(content)
        else: print(".")

    def time_measure(self):
        """returns elapsed time since creation of this instance"""
        self.end = time.time()
        return self.end - self.start

    def get_string_from_yaml_configuration_file_for_Rscript(self, my_list):
        """returns stringified version of the sumstat composition entry of the
        configuration yaml"""
        return 'NEXT'.join([' '.join(l) for l in my_list])

class tsabc2:
    def __init__(self):
        """container class for function of the tsabc pipeline"""
        pass

    def test(self):
        """empty test function"""
        try:
            print("executing tsabc2.test().. successful")
        except:
            sys.exit("could not execute tsabc2.test()")

    def draw_integer(self, my_lower, my_upper, my_distribution_type):
        """draw an integer from a closed interval (inclusive)"""
        my_lower = int(float(my_lower))
        my_upper = int(float(my_upper))

        if my_distribution_type == "uniform":
            return np.random.randint(my_lower, my_upper+1)
        elif my_distribution_type == "exponential":
            sys.exit("""exponential (geometric) distributions are not implemented, yet""")
        else: sys.exit("""only accepting 'uniform' or 'exponential' as distribution type
technically an exponential integer refers to a geometric distribution""")

    def draw_float(self, my_lower, my_upper, my_distribution_type):
        """draw a float from a half-open interval (inclusive)"""
        my_lower = float(my_lower)
        my_upper = float(my_upper)

        if my_distribution_type == "uniform":
            return np.random.uniform(my_lower, my_upper)
        elif my_distribution_type == "exponential":
            return()
        else: sys.exit("""only accepting 'uniform' or 'exponential' as distribution type
technically an exponential integer refers to a geometric distribution""")

    def draw_parameter_from_prior_distribution(self, prior_distribution_definition):
        """this method draws a single parameter value from a defined prior
        distribution and returns this value"""
        assert len(prior_distribution_definition) == 4, "need a tuple containting 4 values: lower, upper, valuetype, distributiontype"

        if prior_distribution_definition[2] == "int":
            return self.draw_integer(*prior_distribution_definition[0:2],
                prior_distribution_definition[3])
        elif prior_distribution_definition[2] == "float":
            return self.draw_float(*prior_distribution_definition[0:2],
                prior_distribution_definition[3])
        else: sys.exit("only accepting 'int' or 'float' as value type")

    def run_transition_to_selfing_model(self, p_s_r, t_p, p_s_a, s_r, t_s, s_a, r, nsample, L, n_dtwf, seed):
        """building the chronology of events and apply to msprime"""
        # transform the model into a time-row matrix
        tables = np.array(self.get_time_sequence_for_transition_to_selfing_model(
            p_s_r, t_p, p_s_a, s_r, t_s, s_a, r))
        times, pop_sizes, sigmas, rec_rates = tables

        # get effective pop sizes, r per time
        matrix = np.array(self.get_N_r_over_time(times, pop_sizes, sigmas, rec_rates))

        # empty workspace
        del times, pop_sizes, sigmas, rec_rates, tables

        # general demographic events; here no pop size, rec size changes
        demography_general = []
        demography_general.append(msprime.SimulationModelChange(n_dtwf,
            "smc_prime"))

        # run the first time segment
        events_phase_0 = [event for event in demography_general
            if event.time < matrix[1,0]]
        recomb_map_phase_0 = msprime.RecombinationMap.uniform_map(
            L, matrix[0,2], num_loci=L)
        print("seed {}: pop_size {} to {} at {}; selfing {} to {} at {}".format(
            seed, p_s_r, p_s_a, t_p, s_r, s_a, t_s))
        print("seed {}: beginning phase {}..".format(seed, 0))
        ts = msprime.simulate(
            sample_size=nsample,
            Ne=matrix[0,1],
            recombination_map=recomb_map_phase_0,
            demographic_events=events_phase_0,
            end_time=matrix[1,0],
            model="dtwf")

        # run the simulation for the remaining consecutive time segments
        for i in range(1, np.size(matrix,0)):
            print("seed {}: beginning phase {}..".format(seed, i))

            # end time definition
            if i < np.size(matrix,0)-1: end_time = matrix[i+1,0]
            else: end_time = float("inf")

            # event filtering
            events_of_this_phase = [event for event in demography_general
                if event.time >= matrix[i,0] and event.time < end_time]

            # fix the end time
            if end_time == float("inf"): end_time = None

            recomb_map__this_phase = msprime.RecombinationMap.uniform_map(
                L, matrix[i,2], num_loci=L)

            ts = msprime.simulate(
                Ne=matrix[i,1],
                recombination_map=recomb_map__this_phase,
                demographic_events=events_of_this_phase,
                end_time=end_time,
                from_ts=ts)

        ts = ts.simplify()
        print("seed {}: trees created..".format(seed))

        return ts

    def get_time_sequence_for_transition_to_selfing_model(self, population_size_recent, t_population_size_change, population_size_ancient, sigma_recent, t_sigma_change, sigma_ancient, recombination_rate):
        """This functions translates the 7 parameters of our abc model into the time series of parameters, that can be applied to the run-the-model function"""
        # times; if both times are equal, will reduce
        times = list(set([0, t_population_size_change, t_sigma_change]))
        times.sort()

        # population sizes and sigma-values over time
        if (t_population_size_change > t_sigma_change):
            pop_sizes = [population_size_recent, population_size_recent,
                        population_size_ancient]
            sigmas = [sigma_recent, sigma_ancient, sigma_ancient]
        elif (t_population_size_change < t_sigma_change):
            pop_sizes = [population_size_recent, population_size_ancient,
                        population_size_ancient]
            sigmas = [sigma_recent, sigma_recent, sigma_ancient]
        elif (len(times) == 2):
            pop_sizes = [population_size_recent, population_size_ancient]
            sigmas = [sigma_recent, sigma_ancient]
        else:
            assert False, "unexpected marginal case"

        # recombination rate assumed to be constant
        rec_rates = [recombination_rate] * len(times)

        return times, pop_sizes, sigmas, rec_rates

    def get_N_r_over_time(self, times: list, pop_sizes: list, sigmas: list, r: list):
        """
        takes real pop sizes and recombination rates and ouputs the effective
        values given a selfing rate
        """
        lam = lambda l: all([True if len(l[i-1]) == len(l[i]) else False
            for i in range(1,len(l))])

        assert lam([times, pop_sizes, sigmas, r]), "time-event-lists-lengths differ"
        assert times[0]==0, "simulation should start at time = 0"

        return(np.column_stack((times, [self.get_N_r_given_sigma(pop_sizes[i], r[i],
            sigmas[i]) for i in range(len(times))])))

    def get_N_r_given_sigma(self, N, r, sigma):
        F_is = sigma / ( 2 - sigma )
        N_new = N * ( 1 - 0.5 * sigma )
        r_new = r * ( 1 - F_is )
        return N_new, r_new

    def get_binned_ld_from_ts(self, ts, num_max_random_mutations, chrom_length, breaks, seed):
        """calculate binned linkage disequilibrium for a subset of the mutations"""
        assert isinstance(ts, tskit.trees.TreeSequence)
        ld_table = self.get_ld(ts, int(float(num_max_random_mutations)), seed)
        ld_binned = self.get_bins_from_ld(ld_table, breaks, chrom_length)
        # ld_distances = ld_binned[:,0]
        ld_vector = ld_binned[:,1]
        return ld_vector

    def get_ld(self, ts_mutated, num_sites, seed):
        np.random.seed(seed)

        m_ids = [mut.site for mut in ts_mutated.mutations()]
        lf = lambda my_num, my_ids: my_num if my_num<=len(my_ids) else len(my_ids)
        num_sites = lf(num_sites, m_ids)
        m_ids_sampled = np.random.choice(m_ids, num_sites, replace = False)

        tables = ts_mutated.dump_tables()

        # mutations
        old_mutations = tables.mutations.copy()
        tables.mutations.clear()
        [tables.mutations.add_row(*old_mutations[i]) for i in m_ids_sampled]

        # sites
        old_sites = tables.sites.copy()
        tables.sites.clear()
        [(tables.sites.add_row(*old_sites[i])) for i in m_ids_sampled]

        # rename the mutations.site (reference to id of site.id)
        tables.mutations.set_columns(
            site=np.array([i for i in range(len(tables.mutations))], dtype=np.int32),
            node=tables.mutations.node,
            derived_state=tables.mutations.derived_state,
            derived_state_offset=tables.mutations.derived_state_offset,
            parent=tables.mutations.parent,
            metadata=tables.mutations.metadata,
            metadata_offset=tables.mutations.metadata_offset)

        tables.sort()

        ts = tables.tree_sequence()
        ts = ts.simplify(reduce_to_site_topology=True)

        ldc = tskit.LdCalculator(ts)
        m_m = ldc.r2_matrix() # took 5 min for 500 mutations on test data

        sites_id = ts.tables.mutations.asdict()["site"]
        pos_dic = ts.tables.sites.asdict()["position"]

        return(np.array([[abs(pos_dic[i] - pos_dic[j]), m_m[i, j]]
            for i, j in itertools.combinations(sites_id, 2)]))

    def get_bins_from_ld(self, ld_table, breaks, chrom_length):
        breaks = [b for b in breaks if b < chrom_length]
        breaks.append(chrom_length)
        df = pd.DataFrame(ld_table)
        df['bins'] = pd.cut(df.iloc[:,0], bins = breaks, include_lowest = True)
        df = df.dropna()
        df = df.groupby(['bins']).mean()
        return(df.to_numpy())

    def get_breaks_for_binnning_tl_true(self, tl_true_binning_scale_factor, max_population_size, num_bins):
        """get the breaks, lower bound 0, upper bound inf returns a np.array with the breaks"""
        breaks = np.array([tl_true_binning_scale_factor*max_population_size*(-1)*np.log(1-i/(num_bins)) for i in range(num_bins)])
        breaks = np.append(breaks, np.inf)

        return breaks

    def get_list_tmrcas_and_lengths(self, ts, num_pairs_max):
        """from ts calculate list of tmrca and lengths for pairwise comparison"""
        assert isinstance(ts, tskit.trees.TreeSequence)

        pair_ts_list = self.get_list_of_tupled_ts(ts, num_pairs_max, 2,
            np.random.randint(0, 4294967296))
        tl_table = [[[tree.interval[1]-tree.interval[0],
            tree.time(tree.root)] if len(tree.roots) == 1
            else warnings.warn("Unrooted tree")
            for tree in tseq.trees()]
            for tseq in pair_ts_list]

        num_pairs = len(tl_table)

        # reformat
        tl = tl_table[0]
        for matrix in tl_table[1:]:
            for row in matrix:
                tl.append(list(row))
        tl = np.array(tl)

        return (tl, num_pairs, tl_table)

    def get_list_of_tupled_ts(self, ts, num_pairs_max, tuple_size, seed):
        all_pairs = itertools.combinations(ts.samples(), tuple_size)
        my_pairs = self.sample_from_iterable(all_pairs, num_pairs_max, seed)
        return [ts.simplify(samples=pair)
            for pair in my_pairs]

    def sample_from_iterable(self, iterable, samplesize, seed):
        """sample from iterable without replacement"""
        np.random.seed(seed)
        return (x for _, x in nlargest(samplesize, ((np.random.uniform(), x)
            for x in iterable)))

    def get_binned_tmrcas_from_true_tl_array(self, tl_true_arr, breaks_age, chrom_length, num_pairs, min_len_break):
        """returns flattened 2d histogram based on given breakpoints
        the breaks for the age follow the PSMC approach,
        while the length uses a simple log10 distribution of breakpoints of same number of bins"""
        LEN = tl_true_arr[:,0]
        AGE = tl_true_arr[:,1]

        # define breakpoints for length to the same dimension of the age breaks
        breaks_len = 10**np.arange(0, np.log10(chrom_length), (np.log10(chrom_length))/(len(breaks_age) - 1))[1:]
        breaks_len = np.sort(np.append([0, np.inf], breaks_len[breaks_len > min_len_break]))

        bins = np.histogram2d(LEN, AGE, [breaks_len, breaks_age], density=False)[0]/num_pairs
        vector_tl_true_binned = bins.flatten('C')

        # verbose information
        if (False):
            print(f"tl_true_vector dimension:\t{len(vector_tl_true_binned)}")
            print(f"tl_true_vector non-zeros:\t{round(100*(vector_tl_true_binned != 0).sum()/len(vector_tl_true_binned), 2)}%")

        return vector_tl_true_binned

    def discretize_list_tmrca_lists(self, tmrcas, breaks):
        """discretize tmrcasby given breaks"""
        return [pd.cut(tmrca, breaks, labels = [i for i in range(len(breaks)-1)], include_lowest=True) for tmrca in tmrcas]

    def transition_matrices_from_tiscretized_tmrcas(self, dtmrcas, nstates):
        """calculate transition matrixes from a list of pandas.cut objects and return them as a list of 2d numpy arrays"""
        return np.array([self.get_transition_matrix(dtmrca, nstates) for dtmrca in dtmrcas])

    def get_transition_matrix(self, transitions, nstates):
        """from a sequence of states calculate the transition matrix"""

        # convert state names into ranks
        # t = transitions.categories
        t = transitions

        n = nstates #number of states

        M = [[0]*n for _ in range(n)]

        for (i,j) in zip(t,t[1:]):
            M[i][j] += 1

        #now convert to probabilities:
        for row in M:
            s = sum(row)
            if s > 0:
                # row[:] = [np.log(s)-np.log(f) for f in row] # log likelihood
                row[:] = [f/s for f in row] # probability

        return np.array(M)

    def get_list_of_pairwise_diversities(self, ts, num_pairs_max, stat_mode, window, tuple_size, window_step, seed):
        """get a list of pairwise diversities from a tree sequence: pandas.DataFrame"""
        assert stat_mode in ["site", "node", "branch"]
        assert window % window_step == 0, "sorry for this"

        span_normalise = False

        pairs = np.array(list(itertools.combinations(ts.samples(), tuple_size)))
        pair_sample = list(self.sample_from_iterable(pairs, num_pairs_max, seed))

        # window definition
        windows = np.arange(0, ts.sequence_length, window)
        if windows[-1]<ts.sequence_length: windows = np.append(windows,
            ts.sequence_length)

        tables = ts.dump_tables()
        ts = tables.tree_sequence()

        div_matrix = np.array([ts.diversity(
            sample_sets=pair,
            windows=windows,
            mode=stat_mode,
            span_normalise=span_normalise) for pair in pair_sample]).transpose()

        div_matrix = pd.DataFrame(div_matrix)

        start_pos = windows[0:len(windows)-1]
        div_matrix.insert(loc=0, column='END_POS', value = np.delete(windows, 0))
        div_matrix.insert(loc=0, column='START_POS', value = start_pos)

        if window - window_step > 0:
            for i in range(int(window/window_step)-1):

                # increase window breakpoints by steps
                windows = windows + window_step
                if windows[-1] > ts.sequence_length:
                    windows = np.delete(windows, -1)
                if windows[-1] < ts.sequence_length:
                    windows = np.append(windows, ts.sequence_length)
                if windows[0] != 0:
                    windows = np.append(0, windows)
                while windows[2] - windows[1] < window:
                    windows = np.delete(windows, 1)

                div = np.array([ts.diversity(
                    sample_sets=pair,
                    windows=windows,
                    mode=stat_mode,
                    span_normalise=span_normalise) for pair in pair_sample]).transpose()
                div = pd.DataFrame(div)
                start_pos = windows[0:len(windows)-1]
                div.insert(loc=0, column='END_POS', value = np.delete(windows, 0))
                div.insert(loc=0, column='START_POS', value = start_pos)

                # remove the smaller windows at the beginning and the end
                div = div.drop([0, len(div.index)-1])

                div_matrix = div_matrix.append(div)
        else:
            if (False): print("window_step equals window_size")

        div_matrix.sort_values(by=['START_POS'], inplace=True)
        return div_matrix

    def sample_from_iterable(self, iterable, samplesize, seed):
        """sample from iterable without replacement"""
        np.random.seed(seed)
        return (x for _, x in nlargest(samplesize, ((np.random.uniform(), x)
            for x in iterable)))

    def get_breaks_for_binning_tm_win(self, discretizing_scale, window_pi_size, mutation_rate, max_pop_size, num_bins):
        """calculate the breakpoints to discretize the diversity for
        calculating the transition matrix of a sliding window pi"""
        prelim_breaks = np.array([-1*discretizing_scale*window_pi_size*mutation_rate*
            max_pop_size*np.log(1-i/(num_bins)) for i in range(num_bins)])

        # remove double interstitious entries
        breaks = [0]
        for b in prelim_breaks:
            if int(b) != int(breaks[-1]):
                breaks.append(b)
        breaks.append(np.inf)

        return np.array(breaks)

    def vector_to_DataFrame(self, my_vector, my_vectors_prefix, index):
        """from a vector return an indexed 1-row pandas dataframe with
        numerated columns' names and specified preifx (e.g. sfs_12)"""
        return pd.DataFrame(data=[my_vector],
            columns=[f"{my_vectors_prefix}_{i}"
                for i in range(len(my_vector))], index=[index])
