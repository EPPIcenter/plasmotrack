# print() shows flagged Infections by class, the 10 worst Infections and the order pairs

    Code
      print(conv)
    Output
      <pt_convergence> 2 Runs, 2 Draws each
        run  draws  iterations  log_posterior   gap
        a        2  1-2                 -3.00  0.00
        b        2  1-2                 -3.00  0.00
      Log-posterior gap: 0.00 (reported, not flagged)
      Quantities: 86, 28 flagged (stuck, R-hat > 1.01 or ess_bulk < 200)
        group                     quantities  flagged  stuck  constant  max_rhat  min_ess_bulk
        log_posterior                      1        0      0         1         -             -
        global                             3        0      0         3         -             -
        expected_false_positives          12        0      0        12         -             -
        expected_false_negatives          12        0      0        12         -             -
        infection_duration                12        0      0        12         -             -
        allele_frequency                   2        0      0         2         -             -
        edge                              12       12     12         0       Inf             -
        order_pair                        12        0      0        12         -             -
        p_latent_parent                   12       10     10         2       Inf             -
        functional                         8        6      6         2       Inf             -
      Infections: 12, 9 flagged (spread > 0.2)
        class               infections
        same_top                     2
        root_flip                    0
        latent_vs_observed           6
        direction_flip               0
        linked_parents               0
        different_parent             1
      10 worst Infections by spread:
        infection  spread  class               worst_rhat  top_parent          
        11          1.000  latent_vs_observed         Inf  a: 1, b: latent     
        12          1.000  different_parent           Inf  a: 2, b: 3          
        10          0.900  latent_vs_observed         Inf  a: 1, b: latent     
        9           0.800  latent_vs_observed         Inf  a: 1, b: latent     
        8           0.700  latent_vs_observed         Inf  a: 1, b: latent     
        7           0.600  latent_vs_observed         Inf  a: 1, b: latent     
        6           0.500  latent_vs_observed         Inf  a: 1, b: latent     
        5           0.400  same_top                   Inf  a: latent, b: latent
        4           0.300  same_top                   Inf  a: latent, b: latent
        3           0.200  -                          Inf  a: latent, b: latent
      Order pairs: 12, 0 flagged (spread > 0.2)

