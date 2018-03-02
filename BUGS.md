# Fri 02 Mar 2018 02:22:41 PM CET:

    pmraic.pl --maui --phyml=/cluster/bin/phyml data.dat

This did not listen to `--phyml`, but tried to locate the binary locally (failing)

# Fri 02 Mar 2018 02:23:40 PM CET

Handle `--bash=6 -cpu=4` to run forked runs (6X4)

# Fri 02 Mar 2018 02:24:34 PM CET

Add SLURM support?

# Fri 02 Mar 2018 02:25:17 PM CET

Check degrees of freedom. Compare PartitionFinder:

    # number of free parameters in substitution model, listed as "model+base_frequencies"
    # and the model string for PhyML as the second of the tuple.
    _base_models = {
        "JC"    :   (0+0, "-m 000000 -f '0.25, 0.25, 0.25, 0.25'"),
        "K80"   :   (1+0, "-m 010010 -f '0.25, 0.25, 0.25, 0.25'"),
        "TrNef" :   (2+0, "-m 010020 -f '0.25, 0.25, 0.25, 0.25'"),
        "K81"   :   (2+0, "-m 012210 -f '0.25, 0.25, 0.25, 0.25'"),
        "TVMef" :   (4+0, "-m 012314 -f '0.25, 0.25, 0.25, 0.25'"),
        "TIMef" :   (3+0, "-m 012230 -f '0.25, 0.25, 0.25, 0.25'"),
        "SYM"   :   (5+0, "-m 012345 -f '0.25, 0.25, 0.25, 0.25'"),
        "F81"   :   (0+3, "-m 000000 -f e"),
        "HKY"   :   (1+3, "-m 010010 -f e"),
        "TrN"   :   (2+3, "-m 010020 -f e"),  
        "K81uf" :   (2+3, "-m 012210 -f e"),
        "TVM"   :   (4+3, "-m 012314 -f e"),
        "TIM"   :   (3+3, "-m 012230 -f e"),
        "GTR"   :   (5+3, "-m 012345 -f e")

# Fri 02 Mar 2018 02:25:57 PM CET

From `summarize_files` function description:

    BUG: run_HoH is not initialized correclty if running script using --summarize
    #                Need to set model name, df, and nbranches, before calling pmraic_print
    #                First, make sure initialize_HoH is called, then iterate over all files to be summarized, and:
    #                @models = check_models_arg(@$models);
    #                @run_HoH{@models} = @HoH{@models}; # hash slice to create %run_HoH from a subset of %HoH

