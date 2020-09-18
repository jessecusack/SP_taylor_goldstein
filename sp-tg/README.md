Code instructions
-------------------------

`tools`    - Matlab codes to solve the T-G equation by Bill Smyth.

`data`     - time-mean profiles or synoptic profiles

`codes`    - main programs:          TG_SI_mean_ex?.m FGM.m (seek for the FGMs using time-mean mooring profiles)
	        main subroutine:        DP.m, DP_0.m (basic data processing)   DP_1.m, DP_2.m, DP_3.m (advanced data processing)
	        auxiliary subroutines:  DP_pic.m (plot profiles)  SI_pic_1.m, SI_pic_2.m (plot FGM-relavent stuff)

`results`  - directory store results

`analysis` - directory to make analysis

T-G analysis seeking for FGMs (fastest growing modes - unstable)
-------------------------
Install MATLAB or Octave first (python codes are on the way), change `mdirec` in line 47 to the directory in your pc, then run the example script `TG_SI_mean_ex1.m` in the `./codes/` directory. 

T-G analysis seeking for LWMs (long wave modes - idealy stable)
-------------------------
on the way...
