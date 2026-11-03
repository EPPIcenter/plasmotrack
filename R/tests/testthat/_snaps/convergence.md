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
        p_latent_parent                    1        0      0         1         -             -
        functional                         7        0      0         7         -             -

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
        p_latent_parent                    1        0      0         1         -             -
        functional                         7        0      0         7         -             -

# print() on the golden Runs

    Code
      print(pt_convergence(a))
    Output
      <pt_convergence> 2 Runs, 10 Draws each
        run               draws  iterations  log_posterior   gap
        single-chain         10  25-70               11.17  0.00
        replica-exchange     10  25-70                8.39  2.78
      Log-posterior gap: 2.78 (reported, not flagged)
      Quantities: 48, 47 flagged (stuck, R-hat > 1.01 or ess_bulk < 200)
        group                     quantities  flagged  stuck  constant  max_rhat  min_ess_bulk
        log_posterior                      1        1      0         0      1.35            10
        global                             3        3      0         0      1.11            10
        expected_false_positives           3        3      0         0      1.18            10
        expected_false_negatives           3        3      0         0      1.08            10
        infection_duration                 3        3      0         0      1.29            10
        allele_frequency                  20       20      0         0      2.04            10
        edge                               1        1      0         0      1.11            10
        p_latent_parent                    3        2      0         1      1.16            10
        functional                        11       11      0         0      1.63            10

