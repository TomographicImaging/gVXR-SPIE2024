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


# CT scan acquisition simulation using [gVXR](https://gvirtualxray.sourceforge.io/) and CT reconstruction with [CIL](https://ccpi.ac.uk/cil/)
# 
# This example makes use of two open source libraries fro X-ray imaging. The first one is [gVXR](https://gvirtualxray.sourceforge.io/). It is used to simulate realistic radiographic images from a CAD model. This time **a JSON file is used to set the simulation parameters.** The JSON file format is relatively user friendly. It is much easier to describe the simulation parameters as no code is needed. The second one is [CIL](https://ccpi.ac.uk/cil/). It implements many CT reconstruction algorithms, including the well-known FDK. The details of the CT scan acquisition are given in the table below. Both FDK and SIRT (iterative method) reconstructions were performed. 
# 
# ![Image comparison of the FDK and SIRT reconstructions](../results/dragon-comparison.png)
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
from gvxrPython3 import json2gvxr # Simulate X-ray images

# CT reconstruction using CIL
from gvxrPython3.JSON2gVXRDataReader import *

from cil.io import TIFFWriter
from cil.utilities.display import show2D, show_geometry
from cil.processors import TransmissionAbsorptionConverter
from cil.framework import AcquisitionGeometry, AcquisitionData
from cil.recon import FDK
from cil.optimisation.algorithms import SIRT
from cil.optimisation.functions import IndicatorBox
from cil.plugins.astra.operators import ProjectionOperator


## Set the simulation parameters

# Initialise gVXR using our JSON file
json_fname = "../results/dragon.json"

# MS Windows
if os.name == "nt":
    json2gvxr.initGVXR(json_fname, renderer="EGL")
# MacOS
elif str(os.uname()).find("Darwin") >= 0:
    json2gvxr.initGVXR(json_fname, renderer="OPENGL")
# GNU/Linux
else:
    json2gvxr.initGVXR(json_fname, renderer="EGL")


# Set up the detector
json2gvxr.initDetector(json_fname, verbose=0)


# Plot the energy response of the detector
detector_response = np.array(gvxr.getEnergyResponse("keV"))
plt.figure(figsize= (20,10))
plt.title("Energy response of the detector")
plt.plot(detector_response[:,0], detector_response[:,1])
plt.xlabel('Incident energy: E (in keV)')
plt.ylabel('Detector energy response: $\\delta$(E) (in keV)')
plt.tight_layout()
plt.savefig("../results/dragon-detector-energy_response.png", dpi=20)
plt.savefig("../results/dragon-detector-energy_response.pdf", dpi=600)


# Plot the energy response of the detector
lsf = np.array(gvxr.getLSF())
half_size = len(lsf) // 2
x = np.arange(-half_size, half_size + 1)

plt.figure(figsize= (20,10))
plt.title("One dimensional line spread function (LSF)")
plt.plot(x, lsf)
plt.xlabel('Pixels')
plt.ylabel('Intensity')
plt.tight_layout()
plt.savefig("../results/dragon-detector-LSF.png", dpi=20)
plt.savefig("../results/dragon-detector-LSF.pdf", dpi=600)


# Create a source
json2gvxr.initSourceGeometry(verbose=0)
json2gvxr.initSpectrum(verbose=0);


# Plot the beam spectrum
energy_bins = gvxr.getEnergyBins("keV")
photon_counts = gvxr.getPhotonCountEnergyBins()
plt.figure(figsize=(20,10))
plt.bar(energy_bins, photon_counts)
plt.xlabel('Energy in keV')
plt.ylabel('Photons per pixel')
plt.title('Photon energy distribution')
plt.xlim([0,200])
plt.tight_layout()

plt.savefig("../results/dragon-beam-spectrum.png", dpi=20)
plt.savefig("../results/dragon-beam-spectrum.pdf", dpi=600)


# Load the sample
json2gvxr.initSamples(verbose=0)


# Compute an X-ray image
x_ray_image = np.array(gvxr.computeXRayImage()).astype(np.single) / gvxr.getWhiteImage()


# Display the corresponding X-ray image
plt.imshow(x_ray_image, cmap="gray");
plt.colorbar()
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
json2gvxr.initScan();
angles = json2gvxr.doCTScan()


## Set the CT reconstruction parameters

# Create the JSON2gVXR reader by passing the filename
reader = JSON2gVXRDataReader(file_name=json_fname)

# Read in file, and return a numpy array containing the data
data_original = reader.read()


# Update the font size
font = {'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)


# Use show2D to visualise a projection
show2D(data_original, origin="upper-left")


# Normalisation and linearisation
data_absorption = TransmissionAbsorptionConverter(white_level=data_original.max(), min_intensity=1e-9)(data_original)


# Create the CIL geoemtry
geometry = data_absorption.geometry


# Shutdown the simulation engine
gvxr.terminate()


# Display and save the geometry, does it look like a feasible CT scan set up?
fig = show_geometry(geometry);
fig.save("../results/dragon-geometry.png");
fig.save("../results/dragon-geometry.pdf");


# Print details of the scanning geometry
print(geometry)


# Prepare the data for the reconstruction
acquisition_data = AcquisitionData(data_absorption, deep_copy=False, geometry=geometry)
acquisition_data.reorder(order='tigre')
ig = acquisition_data.geometry.get_ImageGeometry()


# Print details of the reconstructed volume
print(ig)


## Perform the CT reconstruction using the FDK algorithm and save the reconstructed volume

# Perform the FDK reconstruction
fdk =  FDK(acquisition_data, ig)
recon = fdk.run()


# Show a 2D slice
show2D(recon)


# Save the CT volume as a TIFF stack
TIFFWriter(data=recon, file_name=os.path.join("../results/dragon-recons-FDK", "out")).write()


## Perform the CT reconstruction using the SIRT algorithm and save the reconstructed volume

# Create projection operator using Astra-Toolbox.
acquisition_data.reorder(order='astra')
A = ProjectionOperator(ig, geometry, "gpu")


# Create the initial guess
x0 = ig.allocate()

# non-zero constraint
constraint = IndicatorBox(lower=0)

# Instantiate the reconstruction algorithm
sirt = SIRT(initial=x0, operator=A, data=acquisition_data, constraint=constraint, max_iteration=500)


# Perform 200 iterations
sirt.update_objective_interval = 50
sirt.run(500)

recon_sirt_noisy = sirt.solution


# Show a 2D slice
show2D(recon_sirt_noisy)


# Show a comparison
fig = show2D([recon, recon_sirt_noisy, (recon-recon_sirt_noisy).abs()], \
       ['FBP', 'SIRT', 'difference'], \
       cmap="gray", num_cols=3, size=(15,15), origin='bottom-left', fix_range=True);

fig.save("../results/dragon-comparison.png", dpi=300)
fig.save("../results/dragon-comparison.pdf", dpi=600)


# Save the CT volume as a TIFF stack
TIFFWriter(data=recon_sirt_noisy, file_name=os.path.join("../results/dragon-recons-SIRT", "out")).write()

# Show the plots
plt.show();
