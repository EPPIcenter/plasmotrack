# Plasmotrack: Transmission Network Inference

[![DOI](https://zenodo.org/badge/231125107.svg)](https://doi.org/10.5281/zenodo.18202135)

A C++ library and toolset for modeling transmission networks using genetic data from plasmodium falciparum. The core model implementation is located in `src/impl/model/Model/Model.h`.

**Version:** 1.0.0

## Citation

*Citation information will be added upon preprint release. Please check back for the full citation once the manuscript is published.*

If you use Plasmotrack in your research, please cite the associated manuscript (details to be added).

## Quick Start

### Manual Build

```bash
# Configure and build (release)
cmake --preset=release
cmake --build --preset=release

# Run tests
./build/release/test/transmission_networks_tests

# Run CLI tool
./build/release/tools/cli/Model/transmission_networks_model
```

## Table of Contents

- [Citation](#citation)
- [Quick Start](#quick-start)
- [Prerequisites](#prerequisites)
- [Repository Setup](#repository-setup)
- [Conan Dependency Management](#conan-dependency-management)
- [CMake Build Process](#cmake-build-process)
- [Building the Model](#building-the-model)
- [Running Models](#running-models)
  - [Running the CLI Tool](#running-the-cli-tool)
  - [Input Format](#input-format)
  - [IDP File Format](#idp-file-format)
  - [Output Structure](#output-structure)
  - [Interpreting Results](#interpreting-results)
- [Model Overview](#model-overview)
- [Troubleshooting](#troubleshooting)
- [Development Workflow](#development-workflow)

## Prerequisites

### System Requirements

- **CMake**: Version 4.2.1 or higher
- **C++ Compiler**: C++20 compatible compiler (GCC 10+, Clang 12+, or MSVC 2019+)
- **Python**: Version 3.12 or higher (for Python tools)
- **Build Tools**: Ninja (recommended) or Make

### Required Tools

- **Conan**: Version 2.0.5 or higher (for dependency management)
- **CMake**: For build configuration
- **Git**: For cloning the repository

### Optional Dependencies

- **OpenMP**: For parallel computation support

### Installing Conan

Conan 2.x can be installed via pip:

```bash
pip install conan
```

Verify the installation:

```bash
conan --version
```

## Repository Setup

### Cloning the Repository

```bash
git clone <repository-url>
cd transmission_nets
```

### Repository Structure

Key directories:

- `src/` - Core library source code
  - `src/impl/model/Model/` - Main model implementation (`Model.h`, `Model.cpp`)
  - `src/core/` - Core utilities and data structures
  - `src/model/` - Model components and processes
- `test/` - Test suite (GoogleTest)
- `tools/cli/Model/` - Command-line interface tool
- `python/` - Python simulation tools
- `cmake/` - CMake modules and utilities
- `lib/` - Third-party libraries (spline, imgui)

## Conan Dependency Management

This project uses Conan 2.x for dependency management. Dependencies are automatically installed during the CMake configuration phase via the `conan_provider.cmake` mechanism.

### Dependencies

The following dependencies are managed by Conan (defined in `conandata.yml`):

- `eigen/3.4.0` - Linear algebra library
- `fmt/9.1.0` - Formatting library
- `nlohmann_json/3.11.3` - JSON library
- `boost/1.83.0` - Boost C++ libraries
- `zlib/1.3.1` - Compression library

### Automatic Dependency Installation

The project includes `conan_provider.cmake` which can automatically install dependencies when CMake's `find_package()` is called. However, this requires CMake 3.24+ for the dependency provider mechanism. 

**Note**: CMake 4.2.1+ is required for the dependency provider mechanism used by Conan.

To enable automatic dependency installation, add this line near the top of your CMakeLists.txt (before any `find_package()` calls):

```cmake
include(conan_provider.cmake)
```

The `conanfile.py` defines:
- Package type: `application`
- Generator: `CMakeDeps`
- Layout: `cmake_layout` (Conan 2.x standard layout)

### Manual Conan Setup (if needed)

If automatic dependency installation fails, you can manually install dependencies:

```bash
conan install . --output-folder=build/conan --build=missing
```

Then configure CMake to use the Conan-generated files:

```bash
cmake -DCMAKE_TOOLCHAIN_FILE=build/conan/build/Release/generators/conan_toolchain.cmake -S . -B build
```

## CMake Build Process

> **📖 For comprehensive build instructions, see [BUILD.md](BUILD.md)**

### Build Directory Structure

The project uses CMake presets for standardized build configurations. All builds are out-of-source and located in the `build/` directory:

- `build/debug/` - Debug build with symbols
- `build/release/` - Optimized release build
- `build/release-coverage/` - Release build with coverage instrumentation
- `build/debug-clang/` - Debug build with Clang compiler
- `build/release-clang/` - Release build with Clang compiler

### CMake Configuration Using Presets

The project uses CMake presets for consistent builds. This is the recommended approach.

#### Configure Using Presets

For a Debug build:

```bash
cmake --preset=debug
```

For a Release build:

```bash
cmake --preset=release
```

For a Release build with coverage:

```bash
cmake --preset=release-coverage
```

For Clang builds:

```bash
cmake --preset=debug-clang
cmake --preset=release-clang
```

#### Build Using Presets

After configuring, build using:

```bash
cmake --build --preset=debug
cmake --build --preset=release
```

Or combine configure and build:

```bash
cmake --preset=release && cmake --build --preset=release
```

#### Manual Configuration (Alternative)

If you prefer manual configuration without presets:

```bash
mkdir -p build/release
cd build/release
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ../..
```

### CMakeLists.txt Structure

The project has a hierarchical CMake structure with centralized dependency management:

1. **Root `CMakeLists.txt`**:
   - Sets C++20 standard and compiler options
   - Integrates Conan for dependency management
   - Includes `cmake/Dependencies.cmake` to find all dependencies
   - Defines project-wide variables and paths
   - Includes subdirectories (src, test, tools)

2. **`cmake/Dependencies.cmake`**:
   - Centralized dependency management
   - Finds all required dependencies (Boost, Eigen3, nlohmann_json, fmt, ZLIB)
   - Finds optional dependencies (OpenMP)
   - Sets up imported targets for use throughout the project

3. **`src/CMakeLists.txt`**:
   - Builds the `transmission_networks` static library
   - Organizes source files by component (core, model, impl)
   - Links dependencies using imported targets from centralized function
   - No duplicate `find_package()` calls

4. **`test/CMakeLists.txt`**:
   - Downloads and builds GoogleTest
   - Creates test executable `transmission_networks_tests`
   - Links against the `transmission_networks` library and dependencies
   - Uses imported targets, no `find_package()` calls

5. **`tools/cli/Model/CMakeLists.txt`**:
   - Builds the `transmission_networks_model` CLI executable
   - Links against the `transmission_networks` library and dependencies
   - Uses imported targets, no `find_package()` calls

### Build Output Locations

After building with presets, outputs are located in:

- **Libraries**: `build/<preset>/lib/libtransmission_networks.a` (or `.so` on Linux)
- **Executables**: 
  - `build/<preset>/test/transmission_networks_tests` - Test suite
  - `build/<preset>/tools/cli/Model/transmission_networks_model` - CLI tool
- **Install targets**: `bin/` and `lib/` in the project root (if installed)

## Building the Model

### Step-by-Step Build Instructions

1. **Configure CMake using presets** (from project root):

```bash
cmake --preset=release
```

During configuration, Conan will automatically:
- Detect your system profile
- Install missing dependencies
- Generate CMake configuration files

2. **Build the library and executables**:

```bash
cmake --build --preset=release
```

Or build specific targets:

```bash
# Build only the library
cmake --build --preset=release --target transmission_networks

# Build only the tests
cmake --build --preset=release --target transmission_networks_tests

# Build only the CLI tool
cmake --build --preset=release --target transmission_networks_model
```

### Build Variants

#### Debug Build

Optimized for debugging with symbols:

```bash
cmake --preset=debug
cmake --build --preset=debug
```

#### Release Build

Optimized for performance:

```bash
cmake --preset=release
cmake --build --preset=release
```

#### Coverage Build

For code coverage analysis:

```bash
cmake --preset=release-coverage
cmake --build --preset=release-coverage
```

### Compiler Flags

The project uses the following compiler flags (defined in root `CMakeLists.txt`):

- **Common (all builds)**: `-Wall -Wextra -Wno-unused-function -Wno-unused-parameter`
- **GCC-specific**: `-Wno-maybe-uninitialized` (GCC only)
- **Conditional flags**:
  - `-Werror`: Only if `TRANSMISSION_NETWORKS_WERROR=ON` (default: OFF)
  - `-fno-omit-frame-pointer`: Only for Debug and RelWithDebInfo builds
  - `-g`: Added automatically by CMake for Debug builds
- **Release builds**:
  - **Default optimization**: `-O2 -DNDEBUG` (both GCC and Clang)
  - **Aggressive optimization**: `-O3` available via `TRANSMISSION_NETWORKS_AGGRESSIVE_OPTIMIZATION=ON`
  - **Native CPU optimization**: `-march=native` available via `TRANSMISSION_NETWORKS_NATIVE_OPTIMIZATION=ON`
- **OpenMP**: `-fopenmp` added automatically if OpenMP is found (linked per-target)

## Running Models

### Running the CLI Tool

The `transmission_networks_model` executable is the main CLI tool for running transmission network inference models using MCMC (Markov Chain Monte Carlo) with replica exchange.

#### Basic Usage

```bash
./build/release/tools/cli/Model/transmission_networks_model \
    --input data/input.json \
    --output-dir results/ \
    --symptomatic-idp data/symptomatic_idp.txt \
    --asymptomatic-idp data/asymptomatic_idp.txt
```

#### Command-Line Options

**Required Options:**
- `--input, -i <file>`: Input JSON file containing network data (see [Input Format](#input-format) below)
- `--output-dir, -o <directory>`: Output directory where results will be written
- `--symptomatic-idp <file>`: File path to Symptomatic IDP (Infection Duration Probability) distribution
- `--asymptomatic-idp <file>`: File path to Asymptomatic IDP (Infection Duration Probability) distribution

**MCMC Parameters:**
- `--burnin, -b <int>`: Number of burn-in steps (default: 5000)
- `--sample, -s <int>`: Total number of sampling steps (default: 10000)
- `--thin, -t <int>`: Thinning interval - number of steps between logged samples (default: 1000)

**Replica Exchange Parameters:**
- `--numchains, -n <int>`: Number of chains to run in replica exchange algorithm (default: 1)
- `--numcores, -c <int>`: Number of CPU cores to use (default: 1)
- `--gradient, -g <float>`: Lower temperature of gradient to use in replica exchange (default: 0.0)

**Other Options:**
- `--seed <long>`: Random seed for reproducibility. Use -1 to generate a random seed (default: -1)
- `--hotload, -h`: Hotload (resume) parameters from the output directory
- `--null-model`: Run the null model (ignores genetic data, uses only temporal constraints)
- `--version`: Display version information
- `--help`: Display help message

#### Example Usage

**Basic run with default parameters:**
```bash
./build/release/tools/cli/Model/transmission_networks_model \
    --input data/network.json \
    --output-dir results/run1/ \
    --symptomatic-idp data/symptomatic_idp.txt \
    --asymptomatic-idp data/asymptomatic_idp.txt
```

**Extended run with custom MCMC parameters:**
```bash
./build/release/tools/cli/Model/transmission_networks_model \
    --input data/network.json \
    --output-dir results/run2/ \
    --symptomatic-idp data/symptomatic_idp.txt \
    --asymptomatic-idp data/asymptomatic_idp.txt \
    --burnin 10000 \
    --sample 50000 \
    --thin 100 \
    --seed 12345
```

**Parallel run with replica exchange:**
```bash
./build/release/tools/cli/Model/transmission_networks_model \
    --input data/network.json \
    --output-dir results/run3/ \
    --symptomatic-idp data/symptomatic_idp.txt \
    --asymptomatic-idp data/asymptomatic_idp.txt \
    --numchains 4 \
    --numcores 4 \
    --gradient 0.1
```

**Resume from previous run (hotload):**
```bash
./build/release/tools/cli/Model/transmission_networks_model \
    --input data/network.json \
    --output-dir results/run1/ \
    --symptomatic-idp data/symptomatic_idp.txt \
    --asymptomatic-idp data/asymptomatic_idp.txt \
    --hotload
```

**Null model (no genetics):**
```bash
./build/release/tools/cli/Model/transmission_networks_model \
    --input data/network.json \
    --output-dir results/null_model/ \
    --symptomatic-idp data/symptomatic_idp.txt \
    --asymptomatic-idp data/asymptomatic_idp.txt \
    --null-model
```

### Input Format

The input file must be a JSON file (optionally gzip-compressed with `.gz` extension) containing network data.

#### JSON Schema

The input JSON file must contain two main sections: `loci` and `nodes`. 

```json
{
  "loci": [
    {
      "locus": "AS1",
      "num_alleles": 8,
      "allele_freqs": [0.0013, 0.0017, 0.001, 0.7533, 0.2027, 0.0389, 0.0006, 0.0005]
    },
    {
      "locus": "AS11",
      "num_alleles": 12,
      "allele_freqs": [0.0006, 0.0001, 0.0001, 0.0199, 0.1868, 0.0365, 0.4397, 0.2626, 0.0352, 0.0081, 0.007, 0.0034]
    }
  ],
  "nodes": [
    {
      "id": "1",
      "observation_time": 234.0737,
      "symptomatic": true,
      "observed_genotype": [
        {
          "locus": "AS1",
          "genotype": "00001000"
        },
        {
          "locus": "AS11",
          "genotype": "000000010000"
        },
        {
          "locus": "AS12",
          "genotype": "00000101"
        }
      ],
      "allowed_parents": ["0"]
    },
    {
      "id": "2",
      "observation_time": 306.2297,
      "symptomatic": true,
      "observed_genotype": [
        {
          "locus": "AS1",
          "genotype": "00010000"
        }
      ],
      "allowed_parents": ["1", "0"]
    }
  ]
}
```

#### Field Descriptions

**Top-level fields:**
- `loci` (required): Array of locus definitions
- `nodes` (required): Array of infection event/node definitions

**Locus object:**
- `locus` (required): String identifier for the locus (e.g., "AS1", "AS11")
- `num_alleles` (required): Integer number of alleles at this locus
- `allele_freqs` (optional): Array of allele frequencies (floats, should sum to approximately 1.0). If not provided, frequencies are estimated from the observed data. Each frequency corresponds to one allele at the locus.

**Node object:**
- `id` (required): String identifier for the node/infection event (e.g., "1", "node1")
- `observation_time` (required): Float representing the time when the infection was observed (must be non-negative). Units are arbitrary but should be consistent across all nodes.
- `symptomatic` (optional): Boolean indicating if the infection is symptomatic (default: `true`). Used to select the appropriate infection duration probability distribution.
- `observed_genotype` (required): Array of observed genetic data objects
- `allowed_parents` (required): Array of node ID strings that can be parents of this node. Use `"0"` to indicate the source population (external infection). This field constrains the possible transmission network structure based on epidemiological or other prior knowledge.

**Observed genotype object:**
- `locus` (required): String matching a locus label from the `loci` array
- `genotype` (required): String of binary characters ('0' or '1') representing presence/absence of each allele. 
  - Length must match `num_alleles` for the corresponding locus
  - '1' indicates the allele is present, '0' indicates absence
  - Missing data can be represented as all zeros (`"000..."`) or an empty string

**Additional Notes:**
- **Genotype format**: Genotypes are strings, not arrays. Each character position corresponds to one allele.
- **Missing data**: Missing genotypes are handled as latent variables and will be inferred during MCMC sampling.
- **Network constraints**: The `allowed_parents` field is critical for constraining the transmission network. Nodes can only have parents from this list, which allows incorporation of temporal, spatial, or other epidemiological constraints.

### IDP File Format

IDP (Infection Duration Probability) files specify the probability distribution for the time between infection and detection. These are simple text files containing comma-separated probability values, one per line.

**Format:**
- Each line contains comma-separated probability values
- Values should be non-negative and typically sum to 1.0 (though normalization may be applied)
- Example for a 10-day distribution:

```
0.05
0.1
0.15
0.2
0.2
0.15
0.1
0.05
0.0
0.0
```

This represents probabilities for days 0-9, where day 0 has probability 0.05, day 1 has 0.1, etc.

**Creating IDP files:**
You can create IDP files from epidemiological data or use theoretical distributions. The length of the distribution should cover the expected range of infection durations in your study.

### Output Structure

The model writes output to the specified output directory with the following structure:

```
output_dir/
├── parameters/
│   ├── mean_coi.csv.gz                    # Mean complexity of infection
│   ├── mean_strains_tx.csv.gz             # Mean number of strains transmitted
│   ├── parent_set_size_prob.csv.gz        # Parent set size probability
│   ├── infection_order.csv.gz             # Infection event ordering
│   ├── eps_pos/                           # False positive rates per node
│   │   ├── node1.csv.gz
│   │   ├── node2.csv.gz
│   │   └── ...
│   ├── eps_neg/                           # False negative rates per node
│   │   ├── node1.csv.gz
│   │   ├── node2.csv.gz
│   │   └── ...
│   ├── infection_duration/                # Infection duration per node
│   │   ├── node1.csv.gz
│   │   ├── node2.csv.gz
│   │   └── ...
│   ├── allele_frequencies/                # Allele frequencies per locus
│   │   ├── L1.csv.gz
│   │   ├── L2.csv.gz
│   │   └── ...
│   ├── genotypes/                         # Latent genotypes per node
│   │   ├── node1/
│   │   │   ├── L1.csv.gz
│   │   │   ├── L2.csv.gz
│   │   │   └── ...
│   │   └── ...
│   └── latent_parents/                    # Latent parent genotypes
│       └── ...
├── stats/
│   └── likelihood.csv.gz                 # Log-likelihood values
└── parent_sets/                          # Parent set distributions
    ├── node1_ps.csv.gz
    ├── node2_ps.csv.gz
    └── ...
```

#### Output File Formats

**Parameter files** (`.csv.gz`):
- Compressed CSV files with one value per line
- Each line represents a sampled value from the MCMC chain
- Files are gzip-compressed to save space

**Parent set files** (`*_ps.csv.gz`):
- CSV format with columns: `parent_set,prob,iter`
- Each row represents a parent set configuration and its probability
- Used to infer transmission network structure

**Likelihood file** (`likelihood.csv.gz`):
- Contains log-likelihood values for each iteration
- Format: `llik,<value>`

### Interpreting Results

1. **Transmission Networks**: Check `parent_sets/` directory for inferred parent-child relationships
2. **Model Parameters**: Review `parameters/` for estimated epidemiological parameters
3. **Convergence**: Monitor `stats/likelihood.csv.gz` to assess MCMC convergence
4. **Genotype Inference**: Examine `parameters/genotypes/` for inferred latent genotypes

### Running Tests

Execute the test suite:

```bash
./build/release/test/transmission_networks_tests
```

Run with verbose output:

```bash
./build/release/test/transmission_networks_tests --gtest_output=xml:test_results.xml
```

### Signal Handling

The model handles interrupt signals gracefully:
- **SIGINT** (Ctrl+C), **SIGTERM**: Gracefully finalizes output and exits
- **SIGUSR1**, **SIGUSR2**: Finalizes current output without exiting
- **SIGQUIT**, **SIGABRT**: Attempts graceful shutdown

When interrupted, the model will:
1. Complete the current iteration
2. Write all buffered output to disk
3. Exit cleanly (output remains consistent)

## Model Overview

### Algorithm

The transmission network inference model uses **Bayesian MCMC (Markov Chain Monte Carlo)** with **replica exchange** (parallel tempering) to infer:

1. **Transmission Network Structure**: Parent-child relationships between infection events
2. **Model Parameters**: Epidemiological parameters (mean COI, transmission probabilities, etc.)
3. **Latent Variables**: Unobserved genotypes, infection times, and parent assignments

### Key Components

- **Observation Process**: Models genetic data observation with false positive/negative rates
- **Transmission Process**: Models how infections are transmitted through the network
- **Source Transmission**: Models transmission from source population
- **Node Transmission**: Models transmission between nodes
- **Prior Distributions**: Bayesian priors for model parameters

### Replica Exchange

When using multiple chains (`--numchains > 1`), the model employs replica exchange (parallel tempering):
- Multiple chains run at different "temperatures" (probability scales)
- Chains periodically exchange states to improve mixing
- Higher temperature chains explore the parameter space more freely
- Lower temperature chains focus on high-probability regions
- Temperature gradient is controlled by `--gradient` parameter

### Null Model

The `--null-model` flag runs a simplified version that:
- Ignores all genetic/genotype data
- Uses only temporal constraints (observation times)
- Useful for baseline comparisons and testing temporal-only inference

### Performance Tips

1. **Start with single chain** (`--numchains 1`) for faster initial runs
2. **Use multiple cores** (`--numcores`) when running multiple chains
3. **Adjust thinning** (`--thin`): Lower values = more samples but larger output files
4. **Monitor convergence**: Check likelihood values stabilize before trusting results
5. **Use hotload** (`--hotload`) to resume long runs that were interrupted

## Troubleshooting

### Common Build Issues

#### Conan Not Found

**Error**: `Conan install failed` or `conan: command not found`

**Solution**:
1. Install Conan: `pip install conan`
2. Verify installation: `conan --version`
3. Ensure Conan is in your PATH

#### CMake Version Too Old

**Error**: `CMake 4.2.1 or higher is required`

**Solution**: Upgrade CMake to version 4.2.1 or higher. On Linux:
```bash
# Using package manager
sudo apt-get install cmake  # Ubuntu/Debian
sudo yum install cmake      # CentOS/RHEL

# Or download from cmake.org
```

#### C++20 Compiler Not Found

**Error**: `C++20 standard requested but compiler does not support it`

**Solution**:
1. Use GCC 10+, Clang 12+, or MSVC 2019+
2. Specify compiler explicitly:
   ```bash
   cmake -DCMAKE_CXX_COMPILER=g++-10 ..
   ```

#### Dependency Resolution Failures

**Error**: `Package 'eigen/3.4.0' not found` or similar

**Solution**:
1. Update Conan remotes:
   ```bash
   conan remote add conancenter https://center.conan.io
   ```
2. Clear Conan cache and retry:
   ```bash
   conan remove "*" -c
   ```
3. Manually install dependencies:
   ```bash
   conan install . --build=missing
   ```

#### OpenMP Not Found

**Error**: `OpenMP not found` (warning, not fatal)

**Solution**: OpenMP is optional but recommended. Install:
- **Linux**: `sudo apt-get install libomp-dev` (Ubuntu/Debian)
- **macOS**: `brew install libomp`
- **Windows**: Included with Visual Studio

### CMake Configuration Errors

#### Missing find_package Targets

**Error**: `Could not find a package configuration file provided by "eigen3"`

**Solution**: This usually means Conan dependencies weren't installed. 

1. **If using CMake 3.24+**: The `conan_provider.cmake` should handle this automatically. Ensure it's included in CMakeLists.txt:
   ```cmake
   include(conan_provider.cmake)
   ```

2. **If using CMake 4.2.1+**: The dependency provider should work automatically. If not, manually run Conan install:
   ```bash
   conan install . --output-folder=build/conan --build=missing
   cmake -DCMAKE_TOOLCHAIN_FILE=build/conan/build/Release/generators/conan_toolchain.cmake -S . -B build
   ```

3. **Alternative**: Use the manual Conan setup method described in the Conan section above

#### Build Type Not Specified

**Warning**: `CMAKE_BUILD_TYPE is not set`

**Solution**: Always specify build type:
```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Compiler Compatibility Issues

#### GCC Version Too Old

**Error**: C++20 features not recognized

**Solution**: Use GCC 10 or newer. Check version:
```bash
g++ --version
```

#### Clang Compatibility

If using Clang, ensure it's configured correctly:
```bash
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..
```

### Linker Errors

#### Undefined References

**Error**: `undefined reference to 'boost::...'`

**Solution**: Ensure all dependencies are properly linked. Check `src/CMakeLists.txt` for required libraries.

#### Library Not Found

**Error**: `cannot find -ltransmission_networks`

**Solution**: Build the library first:
```bash
cmake --build --preset=release --target transmission_networks
```

## Development Workflow

### Recommended Build Configuration

For development, use Debug builds for easier debugging:

```bash
cmake --preset=debug
cmake --build --preset=debug
```

For performance testing, use Release builds:

```bash
cmake --preset=release
cmake --build --preset=release
```

### Using CMake Presets

The project includes `CMakePresets.json` with predefined configurations. List available presets:

```bash
cmake --list-presets
```

Configure and build using presets:

```bash
cmake --preset=release
cmake --build --preset=release
```

### IDE Integration

#### CLion

1. Open the project in CLion
2. CLion will automatically detect CMakeLists.txt
3. Select build configuration from the toolbar
4. Build and run from the IDE

#### VS Code

1. Install the CMake Tools extension
2. Open the project folder
3. Select a kit (compiler) and build type
4. Use the CMake Tools panel to configure and build

#### Generate compile_commands.json

The project sets `CMAKE_EXPORT_COMPILE_COMMANDS=ON`, so `compile_commands.json` is automatically generated in the build directory (e.g., `build/release/compile_commands.json`). This enables:
- Better IDE code completion
- clangd support
- Other language server features

### Continuous Integration

For CI/CD, typical workflow:

```bash
# Install dependencies
pip install conan
conan --version

# Configure and build using presets
cmake --preset=release
cmake --build --preset=release

# Run tests
./build/release/test/transmission_networks_tests
```

### Project Structure Reference

- **Model Implementation**: `src/impl/model/Model/Model.h` - Main model class
- **State Management**: `src/impl/model/Model/State.h` - Model state
- **Configuration**: `src/impl/model/Model/config.h` - Model configuration
- **CLI Tool**: `tools/cli/Model/Model.cpp` - Command-line interface
- **Test Suite**: `test/src/` - Unit tests

---

## Additional Resources

- **Conan Documentation**: https://docs.conan.io/
- **CMake Documentation**: https://cmake.org/documentation/
- **C++20 Reference**: https://en.cppreference.com/w/cpp/20
