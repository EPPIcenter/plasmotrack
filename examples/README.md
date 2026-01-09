# Example Files

This directory contains example input files for testing Plasmotrack.

## Files

- `example_input.json` - Example input JSON file with 2 loci and 3 nodes
- `example_symptomatic_idp.txt` - Example Infection Duration Probability distribution for symptomatic infections
- `example_asymptomatic_idp.txt` - Example Infection Duration Probability distribution for asymptomatic infections

## Quick Test Run

To run a quick test with these example files:

```bash
# Create output directory
mkdir -p examples/output

# Run with minimal parameters (quick test)
./build/release/tools/cli/Model/transmission_networks_model \
    --input examples/example_input.json \
    --output-dir examples/output \
    --symptomatic-idp examples/example_symptomatic_idp.txt \
    --asymptomatic-idp examples/example_asymptomatic_idp.txt \
    --burnin 100 \
    --sample 200 \
    --thin 10
```

**Note:** This is a minimal example for testing purposes. For real analyses, use appropriate burnin, sample, and thin values as described in the main README.
