#!/usr/bin/env python3
import os
import numpy as np

from gvxrPython3 import gvxr # Simulate X-ray images
from gvxrPython3.utils import loadSpectrumSpekpy
from gvxrPython3 import gvxr2json # Simulate X-ray images

# CT reconstruction using CIL
from cil.io import TIFFStackReader, TIFFWriter
from cil.processors import TransmissionAbsorptionConverter
from cil.framework import AcquisitionGeometry, AcquisitionData
from cil.recon import FDK
from cil.optimisation.algorithms import SIRT
from cil.optimisation.functions import IndicatorBox
from cil.plugins.astra.operators import ProjectionOperator

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

# Create a source
print("Set up the beam");
gvxr.setSourcePosition(0.0,  -150.0, 0.0, "cm");
gvxr.usePointSource();
#  For a parallel source, use gvxr.useParallelBeam();

# Set its spectrum, here a monochromatic beam
# 1000 photons of 80 keV (i.e. 0.08 MeV) per ray
# gvxr.setMonoChromatic(80, "keV", 1000);

# Or use a polychromatic beam
# The tube voltage is 160 keV
# The filtration is 1mm of tin (Sn)
# The anode angle is 12 degrees
# mAs is 0,5
# The source to detector distance in 50 cm
loadSpectrumSpekpy(160, filters=[["Sn", 1.0]], th_in_deg=12, mAs=0.5, z=150 - -40);

# Poisson noise will be enable
gvxr.enablePoissonNoise(); # Not needed as mAs was used in the function call above

# Locate the sample STL file from the package directory
path = os.path.dirname(gvxr.__file__);
fname = os.path.join(path, "welsh-dragon-small.stl");

gvxr.loadMeshFile("Dragon", fname, "mm");
gvxr.moveToCentre("Dragon");

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

# Simulate the CT acquisition and save the projections
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
gvxr2json.saveJSON("../results/dragon.json");

# Set the CT reconstruction parameters
# Create the TIFF reader by passing the directory containing the files
reader = TIFFStackReader(file_name="../results/dragon-projs", dtype=np.float32);

# Read in file, and return a numpy array containing the data
data_original = reader.read();

# Normalisation
# Not strictly needed as the data was already corrected
data_normalised = data_original / data_original.max();

# Prevent log of a negative value
data_normalised[data_normalised<1e-9] = 1e-9;

# Linearisation
data_absorption = -np.log(data_normalised);

# The data is stored as a stack of detector images, we use the CILlabels for the axes
axis_labels = ['angle','vertical','horizontal'];

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

# Prepare the data for the reconstruction
acquisition_data = AcquisitionData(data_absorption, deep_copy=False, geometry=geometry);
acquisition_data.reorder(order='tigre');

# Perform the FDK reconstruction
ig = acquisition_data.geometry.get_ImageGeometry();
fdk =  FDK(acquisition_data, ig);
recon = fdk.run();

# Save the CT volume as a TIFF stack
TIFFWriter(data=recon, file_name=os.path.join("../results/dragon-recons-FDK", "out")).write();

# Perform the CT reconstruction using the SIRT algorithm and save the reconstructed volume
# Create projection operator using Astra-Toolbox.
acquisition_data.reorder(order='astra');
A = ProjectionOperator(ig, geometry, "gpu");

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