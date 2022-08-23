The code provided here, fit_cilium_beat_pattern_to_data.m, is a MATLAB[^matlab]/GNU Octave[^octave] script intended to help fit a truncated Fourier series representation to experimental trace data of cilia beating.

## Running the code

To run the code, the source file itself should first be edited, changing the lines
```
file_name = 'example_cilium_trace';
file_extension = '.csv';
```
to provide the name and extension of the relevant data file[^file_path]. After this, the script is run ***interactively*** in the command-line interface.

Specifically, the user will be required to interact with the script at two stages:
1. The user will be prompted to provide an estimate for the number of separate time points describing a single period of the data

   ```
   Input an (integer) estimate for the number of data per period:
   ```
   
   before being asked to confirm the choice or choose another estimate
   
   ```
   Is this a good choice for the number of data per period? Press enter if it is, or enter a new estimate for the number of data per period if not: 
   ```
   
   Graphics will be displayed at each stage to inform the user's decision.
   
   The intent here is to give the user more control over the process, particularly if they will have direct knowledge of the trace frequency. If this is not the case, and the data set is too large or 'noisy' for the user to choose for themselves based on the graphical output, this step could be automated by calculating the autocorrelation.
   
2. The user will be similarly prompted to specify the number of Fourier modes to truncate at. The intent here is to allow the user to select the minimum number of modes required to accurately capture the qualitative nature of the beat, but this step could be automated to choose the minimal number of modes required to achieve some desired tolerance.

Once the number of Fourier modes to retain has been confirmed, the script will run to completion. The resulting coefficients will be saved in ASCII .dat files and the user will be presented with a graphical depiction of the reconstructed beat described by these coefficients.

[^file_path]: By including the path in the variable `file_name` -- e.g. `file_name = 'path/to/my/data/files/example_cilium_trace'` -- the script will run in-place on data in a different directory, if required.

[^matlab]: https://mathworks.com/products/matlab.html

[^octave]: https://octave.org/
