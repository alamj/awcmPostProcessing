
# Scale Adaptive Large Eddy Simulation (SALES)
The SALES framework has been implemented over libMesh+PETSc and OpenFOAM. The libMesh-version is based on wavelet-based higher order discretization. The OpenFOAM-version is based on second-order finite volume discretization, which peserves the skew-symmetry of Navier-Stokes equation. 

SALES adopts three principles: i) implicit filtering of Navier-Stokes equation should preserve the skew-symmetry of nonlinear differential operator ii) adopt an explicit filter to construct the structural for of subgrid-scale stress, and iii) find an optimal eddy viscosity that follows Kolmogorov refined similarity hypothesis. 

The following document aims to help post-processing the results of OpenFOAM-based SALES methodology.

## AWCMViewer - Probes -> HDF5 Quickstart

Probes are fixed sensors placed in the computational domain. Each sensor collects time series of select fields. Due to Langrangian view of fluid flows and Taylor's frozen turbulence hypothesis, each sensor is affected by eddies passing through multiple sensor. In datascience view, it is crucial to project the sensed flow onto a low-dimensional vector space. SALES framework has developed a SURE-WT algorithm for extracting coherent structures from the probe data. Some relevant tools are indicated below, which is a work in-progress. 

### An example: how to read probe data and convert to HDF5 format.

Raw data is located in case directory under postProcessing/probes/time/fields, where fields maybe U, UPrime2Mean, UMean, etc.
We want to save time, probe coordinates, and time series of fields within the case directory as statistics/xdmf/fowf15mwR0.h5. In the following command, post_probes.inp maybe an empty file or the sample file (post_process.inp) available in the case directory. 

/project/def-alamj/shared/bin/v2306/awcmviewer -INP post_probes.inp -analysis probes -fields 'U UPrime2Mean UMean' -time 0 -out fowf15mwR0 -pwd

/project/def-alamj/shared/bin/v2306/awcmviewer -INP post_probes.inp -analysis probes -fields 'U UPrime2Mean UMean' -time 0 -out fowf15mwR0 -pwd -pod 1

The last option [-pod value] indicates to apply a POD filter. If value = 0, only the strongest mode is used. If value = 1, no POD filtering is used. Default value is 1. 

The HDF5 data can be proced with phython (h5py), Matlab, or C++ xtensor. Some helper python code is given in awcmviewerutils.py. The idea is to provide the desried probe coordinate and extract the time series for plotting purposes. 

Below is an example of how to read timeseries of streamwise velocity (t, Ux) at probe (1680, 1320, 150):

path = "/scratch/alamj/WindFarms/fowf15mwR0/statistics/xdmf/fowf15mwR0.h5"

t, u = av.get_time_series(out, [1680, 1320, 150], 3, 0, "U") # get Ux

t, u = av.get_time_series(out, [1680, 1320, 150], 6, 0, "UPrime2Mean") # get xx component of Reynolds stress 

t, u = av.get_time_series(out, [1680, 1320, 150], 6, 1, "UPrime2Mean") # xy

t, u = av.get_time_series(out, [1680, 1320, 150], 6, 2, "UPrime2Mean") # xz

t, u = av.get_time_series(out, [1680, 1320, 150], 6, 3, "UPrime2Mean") # yy

----------------------------------------------------------------------------------------------------------------------------------------------------------

### Parallel Reconstruction and Sampling

awcmviewer can read the entire flow decomposed over many processors without needing the costly reconstruction step. Currently, we select a few fields over a single time step and from which, we can extract fields values on a line connecting pointA to pointB. This can be done from vdi node as:

mpirun -np 16 /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis line -fields 'U UPrime2Mean UMean' -time 3600 -out fowf15mwR0 -pwd -pointA '10 1320 0' -pointB '10 1320 960'

It would be more efficient to submit a batch job for extracting data over multiple lines as:

srun /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis line -fields 'U UPrime2Mean UMean' -time 3600 -out fowf15mwR0 -pwd -pointA '10 1320 0' -pointB '10 1320 960' > line1.out

srun /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis line -fields 'U UPrime2Mean UMean' -time 3600 -out fowf15mwR0 -pwd -pointA '90 1320 0' -pointB '90 1320 960' > line2.out

add more lines as needed

----------------------------------------------------------------------------------------------------------------------------------------------------------

### Sampling on a plane

awcmviewer can extract flow field on a 2D plane and save as xdmf + h5. Note, XDMF/HDF5 is a highly efficient data structure to deal with big-data problems using MPI tools. When the simulation is done over a several millions grid points, it better to submit batch jobs and then plot the HDF5 data, which may be more efficient to loading the entire data into paraview. The batch job can be submitted from compute nodes, and plotted from vdi-nodes. 

srun /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis line -fields 'U UPrime2Mean UMean Q Vort' -time 3600 -out fowf15mwR0 -pwd -pointA '90 1320 150' -pointB '0 0 1' > slice.out

here, pointA is a reference point and pointB is a unit vector, which defines a plane over which the flow is extracted. 

### Reconstruct 3D flow field

Use the following command to save as xdmf + h5 fomat:

srun /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis xdmf -fields 'U UPrime2Mean UMean Q Vort' -time 3600 -out fowf15mwR0 -pwd  > xdmf.out

XDMF/HDF5 will create a text file, containing META data and .h5 file containing the actual data. This these files are large, further post-processing can be done using e.g. c++ xtensor tool with MPI framework. 









