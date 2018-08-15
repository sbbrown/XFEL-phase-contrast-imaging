# XFEL-phase-contrast-imaging
This script creates a forward simulation of a coherent beam as it propagates through a uniform sample with sharp discontinuities. These discontinuities are represented by user-specified phase shift in the target material. The final plot shows the beam profile a specified distance downstream after propagating through the target. Multiple phase shifts at different locations in the target perpendicular to the propagation vector can be introduced and their resulting interference observed in the final output. 

## Motivation
This simulation creates a forward simulation of the expected beam profile at a downstream detector of a coherent light pulse that has passed through a sample with zero, one, or multiple discontinuities (phases shifts) introduced in the plane perpendicular to the x-ray propagation vector. This simulation does not analyze/reconstruct phase contrast images.

The advent of highly-coherent, ultrabright x-ray sources enables imaging through dense solids. When a material undergoes rapid compression, tension, or a change of phase, the localized index of refraction changes; this introduces a phase shift into the portion of an imaging x-ray pulse that traverses the local part of the sample. The magnitude of phase shift introduced in a material can be specified a priori or by measuring the visibility of a recorded data image taken with a coherent beam imaging through a material with sharp discontinuities.

For a full explanation of the setup, mathematical formulation, and analysis of phase contrast imaging, please see 
* **Schropp, Andreas, et al.** "Imaging shock waves in diamond with both high temporal and spatial resolution at an XFEL." Scientific reports 5 (2015): 11089.
* **Nagler, Bob, et al.** "The phase-contrast imaging instrument at the matter in extreme conditions endstation at LCLS." Review of Scientific Instruments 87.10 (2016): 103701.

## Code style
Google python style guide https://github.com/google/styleguide/blob/gh-pages/pyguide.md

## Built with
•	Python 2.7 or newer
•	matplotlib
•	numpy

## Code example
Users can modify the original script to change detector specifications, phase shifts introduced (magnitudes and locations within FOV), and/or plotting options. Alternatively, a user can write her/his own script, import this forward propagation script, and utilize the detectorObject and sampleObject classes themselves. 
```
# Imports
import numpy as np
from matplotlib import pyplot

from forward_pci_simulation import detectorObject
from forward_pci_simulation import sampleObject

# Define variables
FOV        = 200.0          # FOV of propagated beam (micron)                    
R          = 400.0          # distance to detector (cm)               
dpix       = 0.05e-6        # detector pixel size (meters)                
wavelength = 1.5120e-10     # of propagated radiation (meters)                 
dz         = 10e-3          # (meters)                 
    
locations  = [132,114]   # dist. from FOV origin of phase shifts (micron)
shifts     = [1.34,1.11] # values of phase shifts at each location (micron)

# Set up detectors and samples by calling classes
detector1   = detectorObject(R, dpix)
sample1     = sampleObject(FOV, detector1)
    
# Run forward simulation
sample1.intro_phase_shifts(locations, shifts)
sample1.propagate_beam_thru(wavelength, dz, detector1)

# Plot results
# First subplot: uniform, Gaussian beam profile before propagation
pyplot.subplot(2,2,1); pyplot.imshow(abs(sample1.E_in)**2, cmap='gray');
    
# Second subplot: beam profile with phase shifts added, before propagation
pyplot.subplot(2,2,2); pyplot.imshow(np.angle(sample1.E_in)); 
    
# Third subplot: beam profile with phase shifts added, after propagation
pyplot.subplot(2,2,3); pyplot.imshow(abs(sample1.E_out)**2,cmap='gray');  
    
# Fourth subplot: lineout of central region of FOV after propogation
lineout_window = [165/sample1.dist_to_matrix,85/sample1.dist_to_matrix]
pyplot.subplot(2,2,4); 
pyplot.plot((abs(sample1.E_out[lineout_window[1]:lineout_window[0],
                                   1024*0.5])**2))
pyplot.show()    
'''
# Save data to csv
np.savetxt("CentralLineout",(abs(sample1.E_out
                                     [lineout_window[1]:lineout_window[0],
                                      1024*0.5])**2),delimiter="")
'''
```

## Installation
This is a library, so download the script (1 total). To run this python script, navigate in a terminal (mac) or command line (windows) to the folder containing the script and execute the following:
```
python forward_pci_simulation.py
```
## Credits
This code is based upon initial methodology + scripts provided by Bob Nagler. Any and all bugs in this class-based version are the responsibility of Shaughnessy Brennan Brown : )

## License
MIT © Shaughnessy Brennan Brown
