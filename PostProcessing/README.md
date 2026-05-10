# PostProcessing

This folder now contains two separated plotting modules:

- config.py
	- Central default config for paths, positions, variables and state file.

- plot_time_evolution_at_positions.py
	- Plot full-time evolution at user-specified spatial positions.
- plot_ex_time_fft_at_position.py
	- Plot Ex time evolution and FFT spectra at user-specified spatial positions.
- plot_ex_spacetime.py
	- Plot Ex(z,t) spacetime map.
- plot_spatial_profile_for_state.py
	- Plot full-position profile for one specified state file.

A unified entry script is also provided:

- plot_output.py
	- mode=time: call time-evolution module
	- mode=space: call spatial-profile module

## Default paths

- Input: ../output
- Output figures: ../results

You can directly edit config.py to change default values:

- OUTPUT_DIR, RESULTS_DIR
- TIME_POSITIONS_UM
- FFT_POSITIONS_UM
- TIME_VARIABLES
- SPACE_STATE_FILE
- DEFAULT_MODE
- NORMALIZE_OUTPUT
- INITIAL_DENSITY_M3
- VELOCITY_DENSITY_MIN_RATIO

## Normalization mode

Both plotting modules support normalized variables using references:

- velocity: c
- density: n0
- electric field: me*c*omega_pe/e
- magnetic field: me*omega_pe/e

where omega_pe is computed from n0:

- omega_pe = sqrt(n0*e^2/(me*eps0))

Enable normalization either by config or CLI:

- config.py: set NORMALIZE_OUTPUT = True
- CLI: add --normalize
- Optional override: --n0 <value_in_m^-3>

Velocity plots hide vx/vz points where ne/n0 is below VELOCITY_DENSITY_MIN_RATIO.
Use NaN masking so those low-density velocity points are not drawn.

## Usage examples

1) Directly use time-evolution module:

```bash
python PostProcessing/plot_time_evolution_at_positions.py --positions-um 0.5,2.0,4.0
python PostProcessing/plot_time_evolution_at_positions.py --positions-um 0.5,2.0,4.0 --normalize --n0 3e26
python PostProcessing/plot_ex_time_fft_at_position.py --positions-um 2.0,8.0,23.0
python PostProcessing/plot_ex_spacetime.py
python PostProcessing/plot_ex_spacetime.py --z-min-um 0 --z-max-um 30 --time-min 70 --time-max 120
```

2) Directly use spatial-profile module (latest file):

```bash
python PostProcessing/plot_spatial_profile_for_state.py --state-file latest
python PostProcessing/plot_spatial_profile_for_state.py --state-file latest --normalize
```

3) Use unified entry script:

```bash
python PostProcessing/plot_output.py time --positions-um 0.5,2.0,4.0
python PostProcessing/plot_output.py space --state-file state_0018000.csv
python PostProcessing/plot_output.py time --positions-um 0.5,2.0,4.0 --normalize --n0 3e26
```
