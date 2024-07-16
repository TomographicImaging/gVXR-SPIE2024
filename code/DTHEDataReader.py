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
# Authors:
# CIL Developers, listed at: https://github.com/TomographicImaging/CIL/blob/master/NOTICE.txt
# Franck P. Vidal (Science and Technology Facilities Council)


from cil.framework import AcquisitionGeometry #, AcquisitionData, ImageData, ImageGeometry, DataOrder
from cil.io.TIFF import TIFFStackReader

import numpy as np
import pandas as pd
import os
from pathlib import Path
from xml.etree import ElementTree
from tifffile import imread

class DTHEDataReader(object):

    '''
    Create a reader for DTHE files
    
    Parameters
    ----------
    directory_path: str
        file name to read

    normalise: bool, default=True
        normalises loaded projections by detector white level (I_0)

    mode: str: {'bin', 'slice'}, default='bin'
        In bin mode, 'step' number of pixels is binned together,
        values of resulting binned pixels are calculated as average. 
        In 'slice' mode 'step' defines standard numpy slicing.
        Note: in general output array size in bin mode != output array size in slice mode

    fliplr: bool, default = False,
        flip projections in the left-right direction (about vertical axis)

    Notes
    -----
    `roi` behaviour:
        Files are stacked along axis_0. axis_1 and axis_2 correspond
        to row and column dimensions, respectively.
        
        Files are stacked in alphabetic order. 
        
        To skip projections or to change number of projections to load, 
        adjust 'angle'. For instance, 'angle': (100, 300)
        will skip first 100 projections and will load 200 projections.
        
        ``'angle': -1`` is a shortcut to load all elements along axis.
            
        ``start`` and ``end`` can be specified as ``None`` which is equivalent
        to ``start = 0`` and ``end = load everything to the end``, respectively.
        Start and end also can be negative.
    '''
    
    def __init__(self,
                 directory_path: str=None,
                 normalise: bool=True,
                 mode: str="bin",
                 fliplr: bool=False,
                 tiff_directory_path: str=None,
                 binning: int=1
                ):

        # Make sure both the files as follwos exist: os.path.join(data_path, "unireconstruction.xml") and os.path.join(data_path, "geom.csv");
        unireconstruction_fname = os.path.abspath(os.path.join(directory_path, "unireconstruction.xml"));
        if not os.path.exists(unireconstruction_fname):
            raise IOError(unireconstruction_fname + " does not exist");

        geomdata_fname = os.path.abspath(os.path.join(directory_path, "geom.csv"));
        if not os.path.exists(geomdata_fname):
            raise IOError(geomdata_fname + " does not exist");

        # Initialise class attributes to None
        self.directory_path = None
        self.unireconstruction_fname = None
        self.geomdata_fname = None
        self.projection_fname_set = None
        
        self.normalise = normalise
        self.mode = mode
        self.fliplr = fliplr
        self.binning = binning
        self._ag = None # The acquisition geometry object
        
        if tiff_directory_path:
            self.tiff_directory_path = Path(tiff_directory_path);
        else:
            self.tiff_directory_path = None

        # The file name is set
        if directory_path is not None:

            # Initialise the instance
            self.set_up(directory_path=directory_path,
                normalise=normalise,
                fliplr=fliplr)


    def set_up(self,
               directory_path: str=None,
               normalise: bool=True,
               mode: str="bin",
               fliplr: bool=False,
               binning: int=1):

        '''Set up the reader
        
        Parameters
        ----------
        data_path: str
            path of the directory where both unireconstruction.xml and geom.csv are

        normalise: bool, default=True
            normalises loaded projections by detector white level (I_0)

        mode: str: {'bin', 'slice'}, default='bin'
            In bin mode, 'step' number of pixels is binned together,
            values of resulting binned pixels are calculated as average. 
            In 'slice' mode 'step' defines standard numpy slicing.
            Note: in general output array size in bin mode != output array size in slice mode

        fliplr: bool, default = False,
            flip projections in the left-right direction (about vertical axis)
        '''

        # Save the attributes
        self.directory_path = os.path.abspath(directory_path)
        # self.roi = roi
        self.normalise = normalise
        self.mode = mode
        self.fliplr = fliplr
        self.binning = binning

        # Check a file name was provided
        if self.directory_path is None:
            raise ValueError('Directory path to unireconstruction.xml and geom.csv files is required.')

        # Make sure both the files as follwos exist: os.path.join(data_path, "unireconstruction.xml") and os.path.join(data_path, "geom.csv");
        self.unireconstruction_fname = os.path.join(self.directory_path, "unireconstruction.xml");
        if not os.path.exists(self.unireconstruction_fname):
            raise IOError(self.unireconstruction_fname + " does not exist");

        self.geomdata_fname = os.path.join(self.directory_path, "geom.csv");
        if not os.path.exists(self.geomdata_fname):
            raise IOError(self.geomdata_fname + " does not exist");

        # Look for projections
        if not self.tiff_directory_path:
            self.tiff_directory_path = os.path.join(self.directory_path, "Proj");

        if not os.path.isdir(self.tiff_directory_path):
            raise ValueError(f"The projection directory '{self.tiff_directory_path}' does not exist")

            
        # Load the CSV file
        geom_data = self.__load_geom_data()
        
        # Convert the geom_data into SOUV notation
        SOUV_data = self.__geom_data_to_SOUV(geom_data);
        number_of_projections = SOUV_data.shape[0];
        
        # Extract the arrays
        source_position_set = SOUV_data[:, 0, :];
        detector_position_set = SOUV_data[:, 1, :];
        detector_direction_x_set = SOUV_data[:, 2, :];
        detector_direction_y_set = SOUV_data[:, 3, :];

        # Extract the pixel size
        pixel_size_in_mm = self.__SOUV_to_pixel_size(SOUV_data);

        self._ag = AcquisitionGeometry.create_Cone3D_SOUV(
            source_position_set=source_position_set,
            detector_position_set=detector_position_set,
            detector_direction_x_set=detector_direction_x_set,
            detector_direction_y_set=detector_direction_y_set,
            volume_centre_position=None,
            units='mm');

        # Read the first projection to extract its size in nmber of pixels
        first_projection_data = imread(self.projection_fname_set[0])
        projections_shape = (number_of_projections, *first_projection_data.shape)

        # Set the detector shape and size
        # Panel is width x height
        self._ag.set_panel(first_projection_data.shape[::-1], pixel_size_in_mm, origin='top-left')

        # Set the order of the data
        self._ag.set_labels(['angle','vertical','horizontal']);            
            

    def __load_geom_data(self) -> pd.DataFrame:

        # Load the CSV file
        df = pd.read_csv(self.geomdata_fname,
                         sep=";",
                         skiprows=2,
                         index_col=False,
                         names=["fname",
                                "source position (x)", "source position (y)", "source position (z)",
                                "imager centre (x)", "imager centre (y)", "imager centre (z)",
                                "imager v vector (x)", "imager v vector (y)", "imager v vector (z)",
                                "imager u vector (x)", "imager u vector (y)", "imager u vector (z)"]);

        # Make sure every column is processed as floating-point data
        for col in df.columns:

            # Except for the file names
            if col != "fname":
                df[col] = df[col].astype(np.single);

        # Load the projections
        self.projection_fname_set = [];

        for i, tiff_fname in enumerate(df["fname"]):
            split_fname = tiff_fname.split("\\");

            self.projection_fname_set.append(os.path.join(self.tiff_directory_path, split_fname[-1]));
            
        # Read the first projection to extract its size in nmber of pixels
        first_projection_data = imread(self.projection_fname_set[0]);
        projection_shape = first_projection_data.shape;
            
        df["imager u vector (x)"] -= df["imager centre (x)"];
        df["imager u vector (y)"] -= df["imager centre (y)"];
        df["imager u vector (z)"] -= df["imager centre (z)"];

        df["imager v vector (x)"] -= df["imager centre (x)"];
        df["imager v vector (y)"] -= df["imager centre (y)"];
        df["imager v vector (z)"] -= df["imager centre (z)"];

        df["imager v vector (x)"] = -df["imager v vector (x)"];
        df["imager v vector (y)"] = -df["imager v vector (y)"];
        df["imager v vector (z)"] = -df["imager v vector (z)"];

        df["imager u vector (x)"] /= projection_shape[1] / 2.0;
        df["imager u vector (y)"] /= projection_shape[1] / 2.0;
        df["imager u vector (z)"] /= projection_shape[1] / 2.0;

        df["imager v vector (x)"] /= projection_shape[0] / 2.0;
        df["imager v vector (y)"] /= projection_shape[0] / 2.0;
        df["imager v vector (z)"] /= projection_shape[0] / 2.0;

        def swapAxis(df, a, b):
            temp = df[a];
            df[a] = df[b];
            df[b] = temp;

        swapAxis(df, "source position (x)", "source position (z)");
        swapAxis(df, "imager centre (x)", "imager centre (z)");
        swapAxis(df, "imager u vector (x)", "imager u vector (z)");
        swapAxis(df, "imager v vector (x)", "imager v vector (z)");

        swapAxis(df, "source position (y)", "source position (z)");
        swapAxis(df, "imager centre (y)", "imager centre (z)");
        swapAxis(df, "imager u vector (y)", "imager u vector (z)");
        swapAxis(df, "imager v vector (y)", "imager v vector (z)");

        return df;

            
    def __geom_data_to_SOUV(self, geom_data):

        SOUV = [];

        for i in range(geom_data.shape[0]):

            source_position = np.array([
                geom_data["source position (x)"][i],
                geom_data["source position (y)"][i],
                geom_data["source position (z)"][i]
            ], dtype=np.single);

            imager_centre = np.array([
                geom_data["imager centre (x)"][i], 
                geom_data["imager centre (y)"][i], 
                geom_data["imager centre (z)"][i]
            ], dtype=np.single);

            detector_direction_x = np.array([
                geom_data["imager u vector (x)"][i],
                geom_data["imager u vector (y)"][i],
                geom_data["imager u vector (z)"][i]
            ], dtype=np.single);

            detector_direction_y = np.array([
                geom_data["imager v vector (x)"][i],
                geom_data["imager v vector (y)"][i],
                geom_data["imager v vector (z)"][i]
            ], dtype=np.single);

            row = [
                    source_position,
                    imager_centre,
                    detector_direction_x,
                    detector_direction_y
                ];

            SOUV.append(row);

        return np.array(SOUV, dtype=np.single);

    def __SOUV_to_pixel_size(self, SOUV):
        initial_detector_direction_x = np.copy(SOUV[0][2]);
        initial_detector_direction_y = np.copy(SOUV[0][3]);

        initial_pixel_width = np.linalg.norm(initial_detector_direction_x);
        initial_pixel_height = np.linalg.norm(initial_detector_direction_y);

        return [initial_pixel_width, initial_pixel_height];

    def read(self):
        
        '''
        Reads projections and returns AcquisitionData with corresponding geometry,
        arranged as ['angle', horizontal'] if a single slice is loaded
        and ['vertical, 'angle', horizontal'] if more than 1 slice is loaded.
        '''

        # Check a file name was provided
        if self.tiff_directory_path is None:
            raise ValueError('The reader was not set properly.')

        # Create the TIFF reader
        reader = TIFFStackReader()

        reader.set_up(file_name=self.projection_fname_set,
                    #   roi=roi,
                      mode=self.mode)

        ad = reader.read_as_AcquisitionData(self._ag)
              
        if (self.normalise):
            white_level = np.max(ad.array)
            ad.array[ad.array < 1] = 1

            # cast the data read to float32
            ad = ad / np.float32(white_level)
            
        
        if self.fliplr:
            dim = ad.get_dimension_axis('horizontal')
            ad.array = np.flip(ad.array, dim)
        
        return ad

    def load_projections(self):
        '''alias of read for backward compatibility'''
        return self.read()


    def get_geometry(self):
        
        '''
        Return AcquisitionGeometry object
        '''
        
        return self._ag

    def get_geometry(self):
        '''
        Return the acquisition geometry object
        '''
        return self._ag
