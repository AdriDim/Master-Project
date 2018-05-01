# Master Project | Adrien Dimech | Hydrogeophysics

Matlab codes to prepare - process - visualize and interpret 3D time-lapse geolelectrical monitoring of a waste rock pile.

Feel free to visit : https://www.researchgate.net/profile/Adrien_Dimech for more information about my research or contact me for more information and data files : adrien.dimech@gmail.com

## Content of this repository :

#### I. Preliminary approach
- 01- 3D modeling of the waste rock pile (476 lines)
- 02- 3D modeling of the pile for COMSOL (1172 lines)

#### II. Static measurements (2016)
- 03- 3D geoelectrical database to RES2D\3DINV (638 lines)
- 04- 3D visualization of RES3DINV inversion results (293 lines)

#### III. Time lapse measurements (2017)
- 05- Optimized protocols for 3D geoelectrical monitoring (1278 lines)
- 06- Optimized protocols for ABEM Terrameter LS (2158 lines)

#### IV. Data processing and inversions
- 07- Hydrogeological database processing (988 lines)
- 08- 3D modeling of the pile for E4D (694 lines)
- 09- 3D time-lapse ERT monitoring (5220 lines)

#### V. Laboratory measurements
- 10- Laboratory column measurements toolbox (4039 lines)


## Description of the Matlab codes :

### I. Preliminary approach

#### 01- 3D modeling of the waste rock pile (476 lines)
This Matlab code was designed to use surveying of the pile to create a complex 3D model of the pile with both external topography and internal structure and instrumentation positions. The 3D model is used to provide better inversion results for the 3D time-lapse geoelectrical monitoring of the experimental waste rock pile.

#### 02- 3D modeling of the pile for COMSOL (1172 lines)
This Matlab code was designed to create a 3D model of the waste rock pile for the COMSOL Multiphysics modelisation software. External topography and internal structure was used to calculate geometric factors for both standard and optimized protocols. The approach developped to model the pile can be applied to any complex structure to generate a complex 3D high-resolution model for COMSOL Multiphysics with LiveLink for Matlab.

### II. Static measurements (2016)

#### 03- 3D geoelectrical database to RES2D\3DINV (638 lines)
This Matlab code was designed to process 3D geoelectrical database of the experimental waste rock pile to RES2DINV and RES3DINV inversion software. Both surface and borehole electrodes are used for meausurements with standard and optimized protocols. 2D and 3D txt files are created for 2D and 3D inversion with user-selected protocols.

#### 04- 3D visualization of RES3DINV inversion results (293 lines)
This Matlab code was used to load inversion results from RES3DINV inversion software and provide 3D visualizations of resistivity and sensitivity distribution. Several tools are provided to the user to identify the best 3D visualizations of the data. 

### III. Time lapse measurements (2017)

#### 05- Optimized protocols for 3D geoelectrical monitoring (1278 lines)
This Matlab code was used to select optimized protocols for the geoelectrical monitoring of the experimental waste rock pile. The methodology follows the example of the 'COMPARE R' method developed by the British Geological Survey (Stummer et al., 2004; Wilkinson et al., 2006b, 2012; Loke et al.,2014a, 2014b, 2015). Optimized protocols with 1000 configurations provide better coverage and senstivities than standard protocols of 4000 configurations. This Matlab code can be used for any 3D electrode distribution and the number of configurations can be selected by the user. 3D visualizations of the optimized protocol sensitivies are used to assess the quality of the selected protocols.

#### 06- Optimized protocols for ABEM Terrameter LS (2158 lines)
This Matlab code was used to select optimized protocols for any Electrical Resistivity Tomography measurements with surface electrodes in 1D or 2D grids (regular or not). The optimized protocols are identified following the methodology of the 'COMPARE R' method developed by the  British Geological Survey (Stummer et al., 2004; Wilkinson et al., 2006b,  2012; Loke et al.,2014a, 2014b, 2015). Optimized protocols provide better resolution and coverage compared to standard protocols and they need less acquisition time. The Matlab code automatically generates .xml files to upload the optimized protocols on the Terramter LS (ABEM) which can be completed by one or more ES1064C electrode selecter (ABEM).

### IV. Data processing and inversions

#### 07- Hydrogeological database processing (988 lines)
This Matlab code was used to load and process the hydrogeological database measured in the experimental waste rock pile. Both GS3 and MPS sensors are used (Decagon Devices) to monitor water content, temperature, resistivity and suction in the waste rock. Precipitation data are also processed. This code provide useful tools to handle complexe database and to visualize different datasets. Hydrogeological data are then used as a validation of geoelectrical results.

#### 08- 3D modeling of the pile for E4D (694 lines)
This Matlab code was used to load surveying data to create a complex 3D model for the inversion software E4D (Johnson, 2010). The 3D tetrahedron mesh is generated with Tetgen. This code provides 3D visualization of the different geometries and structure you wish to model with Tetgen and helps verifying the consistency of the model to find meshes errors for instance. Some useful tools are also provided to complement and interpolate surveying datasets. 

#### 09- 3D time-lapse ERT monitoring (5220 lines)
This Matlab code was designed to process inversion results (E4D) of the 3D time-lapse geoelectrical monitoring of an experimental waste rock pile. Part I loads geoelectrical database from Excel files and pre-processes the measures. Files of data are created to be inversed with E4D (Johnson, 2010). Part II loads inversion files and post-processes the geophysical images of the waste rock pile overt time in 3D. 1D, 2D, 3D and 4D images are generated from the inversion and geoelectrical data are compared to hydrogeophysical measurements.

### V. Laboratory measurements

#### 10- Laboratory column measurements toolbox (4039 lines)
This Matlab code was used to carry laboratory column measurements with rock samples from the experimental waste rock pile (sand, ilmenite and anorthosite waste rocks). ERT monitoring of the column was used to recover the relationship between bulk electrical resistivity, water electrical resistivity and moisture content. These empirical relationships are then used to calculate the moisture content distribution in the waste rock pile from bulk electrical resistivity monitoring.






