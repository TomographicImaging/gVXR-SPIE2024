#!/usr/bin/env python3
import os
import numpy as np

from gvxrPython3 import gvxr # Simulate X-ray images
from gvxrPython3 import json2gvxr # Simulate X-ray images

# CT reconstruction using CIL
from gvxrPython3.JSON2gVXRDataReader import *

from cil.io import TIFFWriter
from cil.processors import TransmissionAbsorptionConverter
from cil.framework import AcquisitionGeometry, AcquisitionData
from cil.recon import FDK
from cil.optimisation.algorithms import SIRT
from cil.optimisation.functions import IndicatorBox
from cil.plugins.astra.operators import ProjectionOperator

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

# Create a source
json2gvxr.initSourceGeometry(verbose=0)
json2gvxr.initSpectrum(verbose=0);

# Load the sample
json2gvxr.initSamples(verbose=0)

# Compute an X-ray image
x_ray_image = np.array(gvxr.computeXRayImage()).astype(np.single) / gvxr.getWhiteImage()

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

# Simulate a CT scan acquisition
json2gvxr.initScan();
angles = json2gvxr.doCTScan()

# Create the JSON2gVXR reader by passing the filename
data_original = JSON2gVXRDataReader(file_name=json_fname).read()

# Normalisation and linearisation
data_absorption = TransmissionAbsorptionConverter(white_level=data_original.max(), min_intensity=1e-9)(data_original)

# Prepare the data for the reconstruction
acquisition_data = AcquisitionData(data_absorption, deep_copy=False, geometry=data_absorption.geometry);
acquisition_data.reorder(order='tigre');
ig = acquisition_data.geometry.get_ImageGeometry();

# Perform the FDK reconstruction
fdk =  FDK(acquisition_data, ig);
recon = fdk.run();

# Save the CT volume as a TIFF stack
TIFFWriter(data=recon, file_name=os.path.join("../results/dragon-recons-FDK", "out")).write();

# Perform the CT reconstruction using the SIRT algorithm and save the reconstructed volume
# Create projection operator using Astra-Toolbox.
acquisition_data.reorder(order='astra');
A = ProjectionOperator(ig, data_absorption.geometry, "gpu");

# Create the initial guess
x0 = ig.allocate();

# non-zero constraint
constraint = IndicatorBox(lower=0);

# Instantiate the reconstruction algorithm
sirt = SIRT(initial=x0, operator=A, data=acquisition_data, constraint=constraint, max_iteration=500);

# Perform 500 iterations
sirt.update_objective_interval = 50;
sirt.run(500);

recon_sirt_noisy = sirt.solution;

# Save the CT volume as a TIFF stack
TIFFWriter(data=recon_sirt_noisy, file_name=os.path.join("../results/dragon-recons-SIRT", "out")).write();
