
# helper functions for snake

def check_parameter_constraints(config):
    """test the properties of the parameters from the config file
    """
    try:
        # population_sizes_backward_in_time
        pop_size_function = config["population_sizes_backward_in_time"]
        sigma_function = config["selfing_rates_backward_in_time"]

        assert len(pop_size_function) == len(sigma_function)

        for pop_size_fn_inner, sigma_function_inner in zip(pop_size_function,
            sigma_function):

            assert pop_size_fn_inner[0][1] == 0, "first time segment should start at present (=0)"
            assert sigma_function_inner[0][1] == 0, "first time segment should start at present (=0)"

            for b in pop_size_fn_inner:
                assert len(b) == 2, "expecting size, time tuple for pop size function"

            for b in sigma_function_inner:
                assert len(b) == 2, "expection rate, time tuple for sigma function"

            for p, t in pop_size_fn_inner:
                assert isinstance(p, int)
                assert isinstance(t, int)

            for s, t in sigma_function_inner:
                assert isinstance(s, float)
                assert isinstance(t, int)

        return True
    except:
        return False