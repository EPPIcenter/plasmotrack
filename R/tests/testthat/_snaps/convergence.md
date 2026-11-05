# print() shows the Runs, the Draws used, the gap and flagged counts per group

    Code
      print(pt_convergence(a))
    Output
      <pt_convergence> 2 Runs, 4 Draws each
      Runs cut to the shortest: a 4, b 6 Draws
        run  draws  iterations  log_posterior   gap
        a        4  1-4                 -7.00  0.00
        b        4  1-4                 -9.00  2.00
      Log-posterior gap: 2.00 (reported, not flagged)
      Quantities: 17, 2 flagged (stuck, R-hat > 1.01 or ess_bulk < 200)
        group                     quantities  flagged  stuck  constant  max_rhat  min_ess_bulk
        log_posterior                      1        1      0         0      1.89             -
        global                             3        1      0         2      1.08             -
        expected_false_positives           1        0      0         1         -             -
        expected_false_negatives           1        0      0         1         -             -
        infection_duration                 1        0      0         1         -             -
        allele_frequency                   2        0      0         2         -             -
        edge                               0        0      0         0         -             -
        order_pair                         0        0      0         0         -             -
        p_latent_parent                    1        0      0         1         -             -
        functional                         7        0      0         7         -             -
      Infections: 1, 0 flagged (spread > 0.2)
      Infections by spread:
        infection  spread  class  worst_rhat  top_parent          
        1           0.000  -               -  a: latent, b: latent
      Order pairs: 0, 0 flagged (spread > 0.2)

---

    Code
      print(pt_convergence(a, iterations = c(2, 4)))
    Output
      <pt_convergence> 2 Runs, 3 Draws each, iterations 2-4
        run  draws  iterations  log_posterior   gap
        a        3  2-4                 -7.00  0.00
        b        3  2-4                 -9.00  2.00
      Log-posterior gap: 2.00 (reported, not flagged)
      Quantities: 17, 2 flagged (stuck, R-hat > 1.01 or ess_bulk < 200)
        group                     quantities  flagged  stuck  constant  max_rhat  min_ess_bulk
        log_posterior                      1        1      0         0       Inf             -
        global                             3        1      0         2      1.22             -
        expected_false_positives           1        0      0         1         -             -
        expected_false_negatives           1        0      0         1         -             -
        infection_duration                 1        0      0         1         -             -
        allele_frequency                   2        0      0         2         -             -
        edge                               0        0      0         0         -             -
        order_pair                         0        0      0         0         -             -
        p_latent_parent                    1        0      0         1         -             -
        functional                         7        0      0         7         -             -
      Infections: 1, 0 flagged (spread > 0.2)
      Infections by spread:
        infection  spread  class  worst_rhat  top_parent          
        1           0.000  -               -  a: latent, b: latent
      Order pairs: 0, 0 flagged (spread > 0.2)

# print() on the golden Runs

    Code
      print(pt_convergence(a))
    Output
      <pt_convergence> 2 Runs, 10 Draws each
        run               draws  iterations  log_posterior   gap
        single-chain         10  25-70                5.90  3.04
        replica-exchange     10  25-70                8.94  0.00
      Log-posterior gap: 3.04 (reported, not flagged)
      Quantities: 49, 47 flagged (stuck, R-hat > 1.01 or ess_bulk < 200)
        group                     quantities  flagged  stuck  constant  max_rhat  min_ess_bulk
        log_posterior                      1        1      0         0      1.37            10
        global                             3        3      0         0      1.81            10
        expected_false_positives           3        3      0         0      1.06            10
        expected_false_negatives           3        3      0         0      1.13            10
        infection_duration                 3        3      0         0      1.21            10
        allele_frequency                  20       20      0         0      2.66            10
        edge                               1        1      0         0      1.06            10
        order_pair                         1        0      0         1         -             -
        p_latent_parent                    3        2      0         1         1            10
        functional                        11       11      0         0      1.32            10
      Infections: 3, 0 flagged (spread > 0.2)
      Infections by spread:
        infection  spread  class  worst_rhat  top_parent                                    
        3           0.096  -             1.1  single-chain: 1, replica-exchange: 1          
        2           0.002  -            1.06  single-chain: latent, replica-exchange: latent
        1           0.000  -            1.21  single-chain: latent, replica-exchange: latent
      Order pairs: 1, 0 flagged (spread > 0.2)

