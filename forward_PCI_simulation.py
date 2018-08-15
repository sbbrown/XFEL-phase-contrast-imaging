'''
Name:          forward_PCI_simulation.py
Version:       1.2
Last modified: June 16, 2018 (SBB)
Authors:       Shaughnessy Brown (sbbrown@slac.stanford.edu)
               Bob Nagler (bnagler@slac.stanford.edu)
Description:   Creates detector objects that store detector characteristics. 
               Creates sample objects that consist of an initial electric 
               field. Introduces phase shifts at specified locations in the 
               sample objects. Propagates a beam of specified wavelength 
               "through" the sample object to the detector object located a 
               specified distance downstream. Plots the initial electric field, 
               the electric field after the phase shifts are introduced, the 
               electric field after propagation to the downstream detector 
               object, and a lineout along the center of the propagated 
               electric field across a specified spatial window. Uncomment 
               bottom lines to save lineout data to .csv file.

               This particular example introduces two phase shifts into a 
               single sample object to observe how the shifs interact when 
               propagated a set distance downstream. In a 200 micron FOV, the 
               shifts [1.34,1.11] are introduced at locations [132,114].     
Variables:
    FOV        = 200.0       # FOV of propagated beam (micron)
    R          = 400.0       # distance to detector (cm)
    dpix       = 0.05e-6     # detector pixel size (meters)
    wavelength = 1.5120e-10  # of propagated radiation (meters)
    dz         = 10e-3       # (meters)
    locations  = [132,114]   # dist. from FOV origin of phase shifts (micron)
    shifts     = [1.34,1.11] # values of phase shifts at each location (micron)
'''

'''Imports -----------------------------------------------------------------'''
import numpy as np
from matplotlib import pyplot


'''Class Definitions -------------------------------------------------------'''        

class detectorObject(object):
    """Class of detector instances."""
    
    def __init__(self, R, dpix):
        self.radius = R                               
        self.dpix   = dpix                            
        
        
class sampleObject(object):
    """Class of beam instances."""
    
    def __init__(self, FOV, detector):
        self.FOV    = FOV
        
        # Create initial beam intensity map
        xd = 1024;  yd = 1024                         # Matrix sizes
        xc = 512.0; yc = 512.0                        # Matrix sizes
        
        self.I0 = np.zeros([xd,yd]) 
        for i in range(xd):
            for j in range(yd):
                self.I0[i,j] = np.exp(-((((i-xc)**2+(j-yc)**2)/
                               (0.95*detector.radius)**2))**10) 
        
        # Create inital phase field
        self.W0           = np.ones((xd,yd))+0.0*1J   
        
        # Maps real-life spatial scale of FOV to matrix units for simulation
        self.dist_to_matrix = self.FOV/float((np.shape(self.W0)[0]))  
         
    def intro_phase_shifts(self, locations, shifts):
        """Inserts phase shifts into the sample object at specified locations
        NOTE: Need to insert beam shift in order from bottom to top, since each 
        sequential matrix change overwrites the preceeding
        """
        self.locations  = locations
        self.shifts     = shifts
        
        # Specify phase shift created by each change in density
        for i in range(len(self.locations)):
            current_phase_shift     = self.shifts[i]
            current_phase_location  = self.locations[i]
            current_matrix_location = int(current_phase_location/
                                          self.dist_to_matrix)
            self.W0[0:current_matrix_location,:] = np.exp(1J*
            current_phase_shift*np.pi*np.ones((current_matrix_location,1024)))
        
            # Calculate initial electric field
        self.E_in  = np.sqrt(self.I0)*self.W0
          
    def propagate_beam_thru(self, wavelength, dz, detector):
        """Propagates a beam through the sample object a distance dz in vacuum
        """
        k1           = np.linspace(-0.5*2*np.pi/detector.dpix,0.5*2*np.pi/
                                   detector.dpix,np.shape(self.E_in)[0])
        Pfx          = np.exp(-1J*wavelength*dz*k1*k1/(4*np.pi))# 1D phase term
        Pf           = np.outer(Pfx,Pfx)                        # 2D phase term
        E_ft         = np.fft.fftshift(np.fft.fft2(self.E_in))
        self.E_out   = np.exp(1J*2*np.pi/wavelength*dz)*                      \
                       np.fft.ifft2(np.fft.ifftshift(E_ft*Pf))
             
'''Main --------------------------------------------------------------------'''        
if __name__ == "__main__":

    # Define variables
    FOV        = 200.0          # FOV of propagated beam (micron)                    
    R          = 400.0          # distance to detector (cm)               
    dpix       = 0.05e-6        # detector pixel size (meters)                
    wavelength = 1.5120e-10     # of propagated radiation (meters)                 
    dz         = 10e-3          # (meters)                 
    
    locations  = [132,114]   # dist. from FOV origin of phase shifts (micron)
    shifts     = [1.34,1.11] # values of phase shifts at each location (micron)

    # Set up detectors and samples
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