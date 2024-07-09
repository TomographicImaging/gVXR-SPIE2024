#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright 2024 United Kingdom Research and Innovation
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#   Authored by:    Franck Vidal (UKRI-STFC)


# # CT scan acquisition simulation using [gVXR](https://gvirtualxray.sourceforge.io/) and CT reconstruction with [CIL](https://ccpi.ac.uk/cil/)
# 
# This example makes use of two open source libraries fro X-ray imaging. The first one is [gVXR](https://gvirtualxray.sourceforge.io/). It is used to simulate realistic radiographic images from a CAD model. The second one is [CIL](https://ccpi.ac.uk/cil/). It implements many CT reconstruction algorithms, including the well-known FDK. The details of the CT scan acquisition are given in the table below.
# 
# | Parameter | Value |
# |-----------|-------|
# | source-to-object distance (SOD) | 150 cm |
# | object-to-detector distance (ODD) | 40 cm |
# | source-to-detector distance (SDD) | 190 cm |
# | detector resolution | 512 &times; 512 pixels |
# | pixel pitch | 500 &times; 500 &mu;m |
# | scintillator | 500 &mu;m of CsI|
# | energy response of the detector | ![Plot of the energy response of the detector](../results/dragon-detector-energy_response.png) |
# | detector impulse response | ![Plot of the detector impulse response](../results/dragon-detector-LSF.png) |
# | tube anode material | tungsten |
# | tube anode angle | 10&deg; |
# | tube voltage | 160 kV |
# | beam filtration | 1mm of Sn |
# | tube exposure | 0.5 mAs |
# | beam spectrum | ![Plot of the beam spectrum](../results/dragon-beam-spectrum.png) |
# | sample geometry | ![](../results/dragon-wireframe.png) |
# | sample material composition | Ti90Al6V4 |
# | sample material density | 4.43 g/cm<sup>3</sup>|
# | number of projection | 200 |
# | first angle | 0&deg; |
# | last angle | 360&deg; |
# | number of flat images | 60 |

# Import packages
import os
import numpy as np

# Increase the font size in plots
import matplotlib
font = {'weight' : 'bold',
        'size'   : 25}

matplotlib.rc('font', **font)

import matplotlib.pyplot as plt # Plotting

from gvxrPython3 import gvxr # Simulate X-ray images
from gvxrPython3.utils import loadSpectrumSpekpy
from gvxrPython3 import gvxr2json # Simulate X-ray images

# CT reconstruction using CIL
from cil.io import TIFFStackReader, TIFFWriter
from cil.utilities.display import show2D, show_geometry
from cil.processors import TransmissionAbsorptionConverter
from cil.framework import AcquisitionGeometry, AcquisitionData
from cil.recon import FDK
from cil.optimisation.algorithms import SIRT
from cil.optimisation.functions import IndicatorBox
from cil.plugins.astra.operators import ProjectionOperator


## Set the simulation parameters

# Create an OpenGL context
print("Create an OpenGL context")
gvxr.createOpenGLContext();


# Set up the detector
print("Set up the detector");
gvxr.setDetectorPosition(0.0, 40.0, 0.0, "cm");
gvxr.setDetectorUpVector(0, 0, -1);
gvxr.setDetectorNumberOfPixels(512, 512);
gvxr.setDetectorPixelSize(500, 500, "um");

# Set the impulse response of the detector, a convolution kernel
gvxr.setLSF([0.00110698, 0.00122599, 0.00136522, 0.00152954, 0.00172533, 0.00196116, 0.0022487, 0.00260419, 0.00305074, 0.00362216, 0.00436939, 0.00537209, 0.00676012, 0.0087564, 0.01176824, 0.01659933, 0.02499446, 0.04120158, 0.0767488, 0.15911699, 0.24774516, 0.15911699, 0.0767488, 0.04120158, 0.02499446, 0.01659933, 0.01176824, 0.0087564, 0.00676012, 0.00537209, 0.00436939, 0.00362216, 0.00305074, 0.00260419, 0.0022487, 0.00196116, 0.00172533, 0.00152954, 0.00136522, 0.00122599, 0.00110698]);

# Set the scintillator
gvxr.setScintillator("CsI", 500, "um");


# Plot the energy response of the detector
detector_response = np.array(gvxr.getEnergyResponse("keV"));
plt.figure(figsize= (20,10));
plt.title("Energy response of the detector");
plt.plot(detector_response[:,0], detector_response[:,1]);
plt.xlabel('Incident energy: E (in keV)');
plt.ylabel('Detector energy response: $\\delta$(E) (in keV)');
plt.tight_layout();
plt.savefig("../results/dragon-detector-energy_response.png", dpi=20);
plt.savefig("../results/dragon-detector-energy_response.pdf", dpi=600);


# Plot the energy response of the detector
lsf = np.array(gvxr.getLSF());
half_size = len(lsf) // 2;
x = np.arange(-half_size, half_size + 1);

plt.figure(figsize= (20,10));
plt.title("One dimensional line spread function (LSF)");
plt.plot(x, lsf);
plt.xlabel('Pixels');
plt.ylabel('Intensity');
plt.tight_layout();
plt.savefig("../results/dragon-detector-LSF.png", dpi=20);
plt.savefig("../results/dragon-detector-LSF.pdf", dpi=600);


# Create a source
print("Set up the beam");
gvxr.setSourcePosition(0.0,  -150.0, 0.0, "cm");
gvxr.usePointSource();
#  For a parallel source, use gvxr.useParallelBeam();

# Set its spectrum, here a monochromatic beam
# 1000 photons of 80 keV (i.e. 0.08 MeV) per ray
# gvxr.setMonoChromatic(0.08, "MeV", 1000);
# The following is equivalent: gvxr.setMonoChromatic(80, "keV", 1000);

# Or use a polychromatic beam
# The tube voltage is 160 keV
# The filtration is 1mm of tin (Sn)
# The anode angle is 12 degrees
# mAs is 0,5
# The source to detector distance in 50 cm
loadSpectrumSpekpy(160, filters=[["Sn", 1.0]], th_in_deg=12, mAs=0.5, z=150 - -40);

# Poisson noise will be enable
gvxr.enablePoissonNoise(); # Not needed as mAs was used in the function call above


# Plot the beam spectrum
energy_bins = gvxr.getEnergyBins("keV");
photon_counts = gvxr.getPhotonCountEnergyBins();
plt.figure(figsize=(20,10));
plt.bar(energy_bins, photon_counts);
plt.xlabel('Energy in keV');
plt.ylabel('Photons per pixel');
plt.title('Photon energy distribution');
plt.xlim([0,200]);
plt.tight_layout();

plt.savefig("../results/dragon-beam-spectrum.png", dpi=20);
plt.savefig("../results/dragon-beam-spectrum.pdf", dpi=600);


# Locate the sample STL file from the package directory
path = os.path.dirname(gvxr.__file__);
fname = os.path.join(path, "welsh-dragon-small.stl");

gvxr.loadMeshFile("Dragon", fname, "mm");
gvxr.moveToCentre("Dragon");

# Change the sample's colour
# By default the object is white, which is not always pretty. Let's change it to gold.
red = 255 / 255
green = 215 / 255
blue = 0 / 255
gvxr.setColour("Dragon", red, green, blue, 1.0)

# Material properties

# Iron (Z number: 26, symbol: Fe)
# gvxr.setElement("Dragon", 26);
# gvxr.setElement("Dragon", "Fe");

# Liquid water
# gvxr.setCompound("Dragon", "H2O");
# gvxr.setDensity("Dragon", 1.0, "g/cm3");
# gvxr.setDensity("Dragon", 1.0, "g.cm-3");

# Titanium Aluminum Vanadium Alloy
# gvxr.setMixture("Dragon", "Ti90Al6V4");
gvxr.setMixture("Dragon", [22, 13, 23], [0.9, 0.06, 0.04]);
# gvxr.setMixture("Dragon", ["Ti", "Al", "V"], [0.9, 0.06, 0.04]); # Not yet implemented
# gvxr.setDensity("Dragon", 4.43, "g/cm3");
gvxr.setDensity("Dragon", 4.43, "g.cm-3");


# Compute an X-ray image
x_ray_image = np.array(gvxr.computeXRayImage()).astype(np.single) / gvxr.getWhiteImage();


# Display the corresponding X-ray image
plt.imshow(x_ray_image, cmap="gray");
plt.colorbar();
plt.axis('off');


# Change to a white background (it could be useful for putting screenshots in papers)
gvxr.setWindowBackGroundColour(1.0, 1.0, 1.0, -1);


# Interactive visualisation
# The user can rotate the 3D scene and zoom-in and -out in the visualisation window.
# It can be useful to rotate the visualisation of the 3D environment and zoom in/out
# to take the best posible screenshots

# - Keys are:
#     - Q/Escape: to quit the event loop (does not close the window)
#     - B: display/hide the X-ray beam
#     - W: display the polygon meshes in solid or wireframe
#     - N: display the X-ray image in negative or positive
#     - H: display/hide the X-ray detector
# - Mouse interactions:
#     - Zoom in/out: mouse wheel
#     - Rotation: Right mouse button down + move cursor```
# gvxr.renderLoop()


# Take and display a screenshot
gvxr.setZoom(2500);
gvxr.displayScene();
screenshot = gvxr.takeScreenshot();
plt.imshow(screenshot);
plt.axis('off');


## Simulate the CT acquisition and save the projections

# Simulate a CT scan acquisition
gvxr.computeCTAcquisition("../results/dragon-projs", # Where to save the projections
                          "screenshots", # Where to save the screenshots
                          200, # Total number of projections
                          0, # First angle
                          True, # Include the last angle
                          360, # Last angle
                          60, # Number of flat images
                          0, 0, 0, "mm", # Centre of rotation
                          *gvxr.getDetectorUpVector()); # Rotation axis


# Save a JSON file
gvxr2json.saveJSON("dragon.json");


## Set the CT reconstruction parameters

# Create the TIFF reader by passing the directory containing the files
reader = TIFFStackReader(file_name="../results/dragon-projs", dtype=np.float32);

# Read in file, and return a numpy array containing the data
data_original = reader.read();

# The data is stored as a stack of detector images, we use the CILlabels for the axes
axis_labels = ['angle','vertical','horizontal'];


# Use show2D to visualise a projection
show2D(data_original, origin="upper-left");


# Normalisation
# Not strictly needed as the data was already corrected
data_normalised = data_original / data_original.max();

# Prevent log of a negative value
data_normalised[data_normalised<1e-9] = 1e-9;

# Linearisation
data_absorption = -np.log(data_normalised);


# Create the CIL geoemtry
    

geometry = AcquisitionGeometry.create_Cone3D(source_position=gvxr.getSourcePosition("cm"),
                                             detector_position=gvxr.getDetectorPosition("cm"),
                                             detector_direction_x=gvxr.getDetectorRightVector(),
                                             detector_direction_y=gvxr.getDetectorUpVector(),
                                             rotation_axis_position=gvxr.getCentreOfRotationPositionCT("cm"),
                                             rotation_axis_direction=gvxr.getRotationAxisCT());
                                         
# Set the angles, remembering to specify the units
geometry.set_angles(np.array(gvxr.getAngleSetCT()), angle_unit='degree');

# Set the detector shape and size
geometry.set_panel(gvxr.getDetectorNumberOfPixels(), gvxr.getDetectorPixelSpacing("cm"));

# Set the order of the data
geometry.set_labels(axis_labels);

# Set the angles, remembering to specify the units
geometry.set_angles(np.array(gvxr.getAngleSetCT()), angle_unit='degree');

# Set the detector shape and size
geometry.set_panel(gvxr.getDetectorNumberOfPixels(), gvxr.getDetectorPixelSpacing("cm"));


# Shutdown the simulation engine
gvxr.terminate();


# Display and save the geometry, does it look like a feasible CT scan set up?
show_geometry(geometry).save("geometry.png");


# Print details of the scanning geometry
print(geometry);


# Prepare the data for the reconstruction
acquisition_data = AcquisitionData(data_absorption, deep_copy=False, geometry=geometry);
acquisition_data.reorder(order='tigre');
ig = acquisition_data.geometry.get_ImageGeometry();


# Print details of the reconstructed volume
print(ig);


## Perform the CT reconstruction using the FDK algorithm and save the reconstructed volume

# Perform the FDK reconstruction
fdk =  FDK(acquisition_data, ig);
recon = fdk.run();


# Update the font size
font = {'weight' : 'bold',
        'size'   : 12};

matplotlib.rc('font', **font);


# Show a 2D slice
show2D(recon);


# Save the CT volume as a TIFF stack
TIFFWriter(data=recon, file_name=os.path.join("../results/dragon-recons-FDK", "out")).write();


## Perform the CT reconstruction using the SIRT algorithm and save the reconstructed volume

# Create projection operator using Astra-Toolbox.
acquisition_data.reorder(order='astra');
A = ProjectionOperator(ig, geometry, "gpu");


# Create the initial guess
x0 = ig.allocate();

# non-zero constraint
constraint = IndicatorBox(lower=0);

# Instantiate the reconstruction algorithm
sirt = SIRT(initial=x0, operator=A, data=acquisition_data, constraint=constraint, max_iteration=500);


# Perform 200 iterations
sirt.update_objective_interval = 50;
sirt.run(500);

recon_sirt_noisy = sirt.solution;


# Show a 2D slice
show2D(recon_sirt_noisy);


# Show a comparison
show2D([recon, recon_sirt_noisy, (recon-recon_sirt_noisy).abs()], \
       ['FBP', 'SIRT', 'difference'], \
       cmap="gray", num_cols=3, size=(15,15), origin='bottom-left', fix_range=True);


# Save the CT volume as a TIFF stack
TIFFWriter(data=recon_sirt_noisy, file_name=os.path.join("../results/dragon-recons-SIRT", "out")).write();

# Show the plots
plt.show();
