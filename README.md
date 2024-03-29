# resoncalc
Calculate bound states and resonances for potential.

## User manual

### Installation 

```
pip install resoncalc
```

Requirements

```
matplotlib>=3.7.2
numpy>=1.25.1
scipy>=1.11.1
sympy>=1.12
```

### Command line interface

```
usage: resoncalc [-h] [-o OUTPUT] [-v] [-s] [-g] 
                 [-t TITLE] input

Calculate bound states and resonances for potential

positional arguments:
  input                 input file with computation settings

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output directory
  -v, --verbose         verbose mode
  -s, --silent          silent mode
  -g, --generate        generate graphs from data
  -t TITLE, --title TITLE
                        output title, used for generate
```

Typical usage for states calculation and graphs generation
```
resoncalc input/settings.json -o output
resoncalc input/data.csv -o output -g -t test
```

### Sample tests

See folder _samples/settings_ for more samples.

#### With mandatory parameters
```
{
  "potential" : "gaussian",
  "params" : [
    {"start": -0.62, "end": -0.56, "cnt": 13},
    {"start": 0.1,   "end": 0.2,   "cnt": 5}
  ],
  "intervals" : [
    {"start": 0.0,   "end": 9.0,   "elems": 15, 
     "type": "equidistant"},
    {"start": 9.0,   "end": 100.0, "elems": 15, 
     "type": "progressive", "len": 0.6},
    {"start": 100.0, "end": 10000.0, "elems": 15, 
     "type": "progressive", "len": 6.0}
  ]
} 
```

#### With optional parameters
```
{
  "title" : "run1",
  "potential" : "parabolic_gaussian",
  "nquad" : 15,
  "r0" : 0.0,
  "phases" : [40.0, 30.0],
  "prec" : 1e-8,
  "emax" : 1.0,
  "mu" : 1.0,
  "l" : 1,
  "params" : [
    {"list": [0.028]},
    {"start": 0.028, "end": 0.029, "cnt": 10}
  ],
  "params2" : [1.0],
  "intervals" : [
    {"start": 0.0,   "end": 10.0,  "elems": 20, 
     "type": "equidistant"},
    {"start": 10.0,  "end": 150.0, "elems": 15, 
     "type": "progressive", "len": 0.5}
   ],
  "outstates" : ["bound", "resonance"],
  "outfiles" : ["states", "eigenvalues", "potential_grid", 
                "spectrum", "log", "settings"]
} 
```

#### Parameters
- _title_: output directory name, by default according to potential
- _potential_: potential from list
- _nquad_: order of quadrature polynomials
- _r0_: center for ECS method, default 0
- _phases_: 2 phases for ECS method, default 40, 30
- _prec_: states detection precision, default 1e-8
- _emax_: maximum detected energy in atomic units
- _mu_: reduced mass in atomic units
- _l_: secondary quantum number
- _params_: potential parameters definition in interval [start,end] and given count, alternatively list of values       
- _params2_: other potential parameters not changed during calculation
- _intervals_: element definition for FEM-DVR method in interval [start,end] and given count of elements and division
- _outstates_: types of generated states, default all
- _outfiles_: types of generated output files, default all

#### Potentials
See file _potential.py_ for definition. 
- _gaussian_
- _exponential_
- _morse_
- _parabolic_gaussian_
- _parabolic_gaussian2_

You can also add your potential, just create new function and add it to mapping.
Sample settings files for all potential are available in folder _samples/settings_.

### Output

Sample output for states calculation is available in folder _samples/output_.

- _states.csv_: detected bound states and resonances
- _spectrum_par1_par2_l.png_: complex spectrum with highlighted bound states and resonances and eigenvalues for given parameters and both ECS phases
- _potential_par1_par2_l.png_: potential with highlighted bound states and resonances for given parameters
- _eigenvalues_par1_par_l.phase.csv_: eigenvalues for given parameters and phase
- _potential_grid_par1_par2_l_phase.csv_: potential values on grid for given parameters and phase
- _gaussian.json_: copy of input settings
- _log.log_: application log, more detailed in verbose mode

Sample output for graphs generation is available in folde _samples/generate_.

- _gaussian.csv_ input file with resonances for gaussian potential
- _resonances_complex_gaussian.png_: complex spectrum
- _resonances_energy_gaussian.png_: resonance energy
- _resonances_width_gaussian.png_: resonance width