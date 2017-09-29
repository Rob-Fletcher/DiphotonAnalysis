# GP Software Notes
Notes on what will be needed from code that will do the full gaussian process
background estimation

---

# It needs to do:
* Convert histograms to numpy arrays.
    - Data points
    - Array of points on the x-axis (bin centers)
    - error bars (possibly up and down)

* Convert numpy arrays to histograms
    - Error bars
    - x axis range

* Run Gaussian process
    - flags to turn on and off the optimization of hyper-parameters
    - inputs for the hyper-parameters when optimizer turned off.
    - Return plots (or arrays) of the log likelihood landscape
    - Maybe and option to easily change kernels

* Configuration
    - Probably JSON file to set everything up so we dont have a million CL options

* Setup
    - document setup with all of the virtualenv stuff
    - maybe setup a pip package for easy install
