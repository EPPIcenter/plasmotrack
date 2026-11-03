# print() shows one line per Run, the epi match and the engine

    Code
      print(a)
    Output
      <pt_analysis> 2 Runs, 3 Infections, 2 Loci
        run   draws  iterations  windows  thin
        one      10  25-70             3     5
        four     10  25-70             3     5
      Epi: 2 of 3 Infections matched; columns district, travel
      Engine: DuckDB, memory_limit 1.0 GiB, threads 1, temp_directory <temp>

---

    Code
      print(b)
    Output
      <pt_analysis> 1 Run, 3 Infections, 2 Loci
        run           draws  iterations  windows  thin
        single-chain     10  25-70             3     5
      Epi: none
      Engine: closed

