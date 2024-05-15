# SDCSIM - Strand Displacement Cascades simulator
#### Simulator for SIMD||DNA model of computation, as defined in https://doi.org/10.1007/978-3-030-26807-7_12
##### website: https://lukas0025.github.io/sdc_sim/

![](https://lukas0025.github.io/sdc_sim/assets/img/header.gif)

### Building

to set up the application you will at least need make, cmake and g++ with openMP support. If you want to set up a GPU accelerated version, you need NVHPC in addition.

To build an application without GPU acceleration support, just invoke the make command.

```sh
make
```

To build a simulator with support for GPU acceleration, use the gpu argument.

```sh
make gpu
```

> The repository also contains shell files for loading dependencies in the module environment. These scripts are `modules.sh` and `modulesGpu.sh`. If your environment supports modules, you can use these scripts.

### Using

SDCSIM is a simple console application that is controlled using arguments and parameters. To successfully run the simulation, it is necessary to pass one argument to the simulator, which represents the name of the file describing the simulated system. These files are SDCASM files that use the special language described below.

```sh
./build/sdcsim <SDCASM>
```


Using other parameters, it is possible to adjust the behavior of the simulation. These parameters are as follows:

```sh
Usage: sdcsim [--help] [--version] [--format VAR] [--silent] [--decode] [--all] [--nucleotides] [--gpu] [--color VAR] [--break VAR] [--time VAR] [--strands VAR] [--temperature VAR] [--spaceing VAR] sdcasm

Positional arguments:
  sdcasm             File with SDC Assembly code of SDC program. [required]

Optional arguments:
  -h, --help         shows help message and exits 
  -v, --version      prints version information and exits 
  -f, --format       Define output fromat of simulation. [nargs=0..1] [default: "ascii"] {ascii, raw, svg, dummy, assembly}
  --silent           Print only necessary outputs 
  -d, --decode       Print decoded version of registers (by macros) 
  -a, --all          Show all instruction interactions 
  -n, --nucleotides  Run sumulation on nucleotides level 
  -g, --gpu          Run sumulation on nucleotides level on GPU 
  -c, --color        Select color scheme for svg output [nargs=0..1] [default: "domain"] {domain, chain, black}
  -b, --break        Break before instruction [nargs=0..1] [default: 2147483647]
  -t, --time         registr expose time to instruction [nargs=0..1] [default: -1]
  -S, --strands      Number of inserted instruction strand copies [nargs=0..1] [default: 100]
  -T, --temperature  temperature of solution [nargs=0..1] [default: 25]
  -s, --spaceing     Space char/s between domains in ascii formats [nargs=0..1] [default: ""]
```

### Language

The language used in SDCASM files is designed to easily write system simd||dna. Each file consists of at least two sections, which are the data section and the instruction section. Every file must have them to be usable. The data section describes the source registers and the instruction section describes the set of instructions.

```
data:
  <code>

instructions:
  <code>
```

Molecules representing individual registers and instructions are written in these sections on single lines. ASCII characters are used to write DNA molecules, which represent domains, and brackets `[]` representing dsDNA, brackets `{}` representing lower ssDNA and brackets `<>` representing upper ssDNA. Each ascii character can be followed by a `*` character representing that this domain is a complementary domain. The last special character is the `.`, which represents the junction of the top ssDNA with the dsDNA. The `#` character is used for single-line comments.

```
data:
  {A}[BN*].<C*>
  #
  #        /
  #    - / 
  #    | | 
  #  - - - 
  #  A B N*
  #

  {ABC}
  #
  #  - - - 
  #  A B C
  #

instructions:
  {ABC} {CDC} # the space serves to separate the individual DNA molecules in the instruction
```


There are two more sections in the language which are optional. One of the sections is `define`, it is used to define macros, which are used for simple replacement in the rest of the code. Macros always consist of two parts separated by a space. In all other sections, the occurrence of the left part will be replaced by the right part. The second section is the `domains` section, which is used to define the nucleotide representation of domains. The domain representation definition is again composed of two parts separated by a space. The right part represents the nucleotide representation and the left the domain.

```
#
# RULE 110 cellular automaton implementation in DNA|SIMD
# @autor Lukáš Plevač <xpleva07@vutbr.cz>
# @date 11.21.2023
#

domains:
    A CACATAC
    B CTTTACA
    C TCTCCT
    D ACACACA
    E CAAAACT
    F TCAAATC
    G TCCACT

define:
    0 [ABC][DE]
    1 {A}[BCDE]

data: # O(6)
    1001111010

instructions:
    {D*E*A*F*}      # mark 01
    {D*E*A*B*C*G*}  # mark 11
    {DEABCG}        # remove mark 11
    {A*B*C*} {D*E*} # write 0
    {DEAF}          # remove mark 01
    {B*C*D*E*}      # write 1
```

### Support scripts

The repository also contains scripts to facilitate work with the simulation.

#### stats.py

This script is used to collect statistical data from the nucleotide simulation for various simulation parameters. This script can be used, for example, to determine the temperature resistance of the system or its time requirements.

```
usage: stats.py [-h] [-T TEMPERARURE] [-t TIME] [-i INDIVIDUALS] [-r TEMPERARURE_RANGE] [--time_range TIME_RANGE]
                [--time_step TIME_STEP] [--strands STRANDS] [--strands_range STRANDS_RANGE]
                [--strands_step STRANDS_STEP] [-c CSV] [-p PLOT] [-s SAMPLES] [--target TARGET]
                sdcsim sdcasm

Simple statistical evaulator for sdcsim

positional arguments:
  sdcsim                Path to sdcsim simulator
  sdcasm                Path to simulated adcsm file

options:
  -h, --help            show this help message and exit
  -T TEMPERARURE, --temperarure TEMPERARURE
                        Center of tested temperatures
  -t TIME, --time TIME  Center of tested times
  -i INDIVIDUALS, --individuals INDIVIDUALS
                        Number of individuals
  -r TEMPERARURE_RANGE, --temperarure_range TEMPERARURE_RANGE
                        Range for testet temperatures (Temp +- (range/2))
  --time_range TIME_RANGE
                        Range for testet times (Time +- (range/2))
  --time_step TIME_STEP
                        Step between times
  --strands STRANDS     Number of strands
  --strands_range STRANDS_RANGE
                        Range of strands number
  --strands_step STRANDS_STEP
                        step in strends number
  -c CSV, --csv CSV     Store test results in csv file
  -p PLOT, --plot PLOT  Store ratio plot to file
  -s SAMPLES, --samples SAMPLES
                        How many samples per 1deg
  --target TARGET       Expected target output of simulation

Created by Lukas Plevac <xpleva07> BUT FIT
```

> This script requires the following dependencies
>
> * `python3`
> * `matplotlib`
> * `numpy`
> * `packaging`
> * `pillow`
> * `six`
> * `joblib`
#### evo.py


This script is used for the evolutionary design of SIMD||DNA programs using the specification of inital and final configurations. This script is more of a proof of concept and does not allow you to design really complex programs because it is not debugged for them.

```
usage: evo.py [-h] [--pop POP] [--pmut PMUT] [--tour TOUR] [--mgenes MGENES] [--gen GEN] [--drand DRAND]
              [--profile]

Simple evolutionary SIMD||DNA designer

options:
  -h, --help       show this help message and exit
  --pop POP        Size of population
  --pmut PMUT      Probability of mutation for one gene
  --tour TOUR      Base of turnament select GA
  --mgenes MGENES  Number of mutated gene values
  --gen GEN        Number of generations
  --drand DRAND    Randomize data every N generations if set
  --profile        Do profile with output csv

Created by Lukas Plevac <xpleva07> BUT FIT
```

> This script requires the following dependencies
>
> * `python3`
> * `numpy`
### Example SIMD||DNA programs

Inside a folder `sdcasm` in this repository you will find sample sdcasm files that can be used for inspiration or to understand the concepts of SIMD||DNA systems.
