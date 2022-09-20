The main code provided here, fit_cilium_beat_pattern_to_data.m, is a MATLAB[^matlab]/GNU Octave[^octave] script intended to help fit a truncated Fourier series representation to experimental trace data of cilia beating.

## Running the code

To run the code, the script source file fit_cilium_beat_pattern_to_data.m itself should first be edited, changing the lines
```
file_name = 'example_cilium_trace';
file_extension = '.csv';
```
to provide the name and extension of the relevant data file[^file_path]. After this, the script is run ***interactively*** in the command-line interface (CLI).

Specifically, the user will be required to interact with the script at two stages:
1. The user will be prompted to provide an estimate for the number of separate time points describing a single period of the data

   ```
   Input an (integer) estimate for the number of data per period (data autocorrelation suggests XXX):
   ```
   
   before being asked to confirm the choice or choose another estimate
   
   ```
   Is this a good choice for the number of data per period? Press enter if it is, or enter a new estimate for the number of data per period if not: 
   ```
   
   Graphics will be displayed at each stage to inform the user's decision.
   
   The intent here is to give the user more control over the process, particularly if they will have direct knowledge of the trace frequency. If this is not the case, and the data set is too large or 'noisy' for the user to easily choose for themselves based on the graphical output, an automatic estimate is offered based on calculating the autocorrelation of the data.
   
2. The user will be similarly prompted to specify the number of Fourier modes to truncate at. The intent here is to allow the user to select the minimum number of modes required to accurately capture the qualitative nature of the beat, but this step could be automated to choose the minimal number of modes required to achieve some desired tolerance.

Once the number of Fourier modes to retain has been confirmed, the script will run to completion. The resulting coefficients will be saved in ASCII .dat files and the user will be presented with a graphical depiction of the reconstructed beat described by these coefficients.

## Examples

Included are three example traces stored in .csv files. These were obtained by tracing, by hand, over videos of beating cilia from _Acropora millepora_ larvae -- see Poon _et al._ (2022)[^paper_link] and its supplementary material for more information.

Running the script with
```
file_name = 'trace_realunits_1';
file_extension = '.csv';
```
and selecting 13 data per period and 4 Fourier modes (i.e. the CLI should appear as follows) will produce the coefficients, and hence the reconstructed beat pattern, given in Poon _et al._ (2022).
```
>> fit_cilium_beat_pattern_to_data
Input an (integer) estimate for the number of data per period (data autocorrelation suggests 13): 13
Is this a good choice for the number of data per period? Press enter if it is, or enter a new estimate for the number of data per period if not: 
Enter a number of modes (an integer between 1 and 7) to retain: 4
Does this number of modes offer a good fit? Press enter if it does, or enter a new choice for the number of modes otherwise: 
```

[^file_path]: By including the path in the variable `file_name` -- e.g. `file_name = 'path/to/my/data/files/example_cilium_trace'` -- the script will run in-place on data in a different directory, if required.

[^matlab]: https://mathworks.com/products/matlab.html

[^octave]: https://octave.org/

[^paper_link]: https://doi.org/10.1101/2022.09.19.508546
