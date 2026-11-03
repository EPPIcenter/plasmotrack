# an input whose SHA-256 differs from the manifests' errors

    Code
      pt_read(golden_store(), input = wrong)
    Condition
      Error in `pt_read()`:
      ! `input` is not the input of Run "single-chain": its SHA-256 differs.
      x '<dir>/input.json': a488f8cf8924
      i Run manifests: cd84e24b6713

# a missing or mismatched input in cli_arguments gives NA times and one message

    Code
      a <- pt_read(golden_store())
    Message
      No input JSON found; observation times are NA and `pt_candidates()` errors.
      x Run single-chain holds no copy of its input
      x <input>, named in Run single-chain's cli_arguments: not found
      i Give the input with `pt_read(..., input = )`.

# pt_candidates() errors when pt_read() found no input

    Code
      pt_candidates(a)
    Condition
      Error in `pt_candidates()`:
      ! This Analysis has no input JSON, which holds the candidate pairs.
      i Give it with `pt_read(..., input = )`.

# the Store's copy may be gzipped, and a copy with the wrong hash gives NA times and one message

    Code
      b <- pt_read(other)
    Message
      No input JSON found; observation times are NA and `pt_candidates()` errors.
      x Run single-chain's copy input.json: its SHA-256 differs from the manifest's
      x <input>, named in Run single-chain's cli_arguments: not found
      i Give the input with `pt_read(..., input = )`.

# print() says where the input came from, or that it was not found

    Code
      print(a)
    Output
      <pt_analysis> 1 Run, 3 Infections, 2 Loci
        run  status    draws  iterations  windows  thin
        one  finished     10  25-70             3     5
      Epi: none
      Input: <golden input>
      Engine: DuckDB, memory_limit 1.0 GiB, threads 1, temp_directory <temp>
    Code
      print(b)
    Output
      <pt_analysis> 1 Run, 3 Infections, 2 Loci
        run  status    draws  iterations  windows  thin
        one  finished     10  25-70             3     5
      Epi: none
      Input: not found; observation times are NA (give pt_read() `input =`)
      Engine: DuckDB, memory_limit 1.0 GiB, threads 1, temp_directory <temp>

