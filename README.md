
# Scale Adaptive Large Eddy Simulation (SALES)
The SALES framework has been implemented over libMesh+PETSc and OpenFOAM. The libMesh-version is based on wavelet-based higher order discretization. The OpenFOAM-version is based on second-order finite volume discretization, which peserves the skew-symmetry of Navier-Stokes equation. 

SALES adopts three principles: i) implicit filtering of Navier-Stokes equation should preserve the skew-symmetry of nonlinear differential operator ii) adopt an explicit filter to construct the structural for of subgrid-scale stress, and iii) find an optimal eddy viscosity that follows Kolmogorov refined similarity hypothesis. 

The following document aims to help post-processing the results of OpenFOAM-based SALES methodology.

## AWCMViewer - Probes -> HDF5 Quickstart


### An example: how to read probe data and convert to HDF5 format.

Raw data is located in case directory under postProcessing/probes/time/fields, where fields maybe U, UPrime2Mean, UMean, etc.
We want to save time, probe coordinates, and time series of fields within the case directory as statistics/xdmf/fowf15mwR0.h5

/project/def-alamj/shared/bin/v2306/awcmviewer -INP post_probes.inp -analysis probes -fields 'U UPrime2Mean UMean' -time 0 -out fowf15mwR0 -pwd

/project/def-alamj/shared/bin/v2306/awcmviewer -INP post_probes.inp -analysis probes -fields 'U UPrime2Mean UMean' -time 0 -out fowf15mwR0 -pwd -pod 1

The last option [-pod value] indicates to apply a POD filter. If value = 0, only the strongest mode is used. If value = 1, no POD filtering is used. 

Read timeseries of streamwise velocity (t, Ux) at probe (1680, 1320, 150):

path = "/scratch/alamj/WindFarms/fowf15mwR0/statistics/xdmf/fowf15mwR0.h5"
t, u = av.get_time_series(out, [1680, 1320, 150], 3, 0, "U") // Ux
t, u = av.get_time_series(out, [1680, 1320, 150], 6, 0, "UPrime2Mean") // xx 
t, u = av.get_time_series(out, [1680, 1320, 150], 6, 1, "UPrime2Mean") // xy
t, u = av.get_time_series(out, [1680, 1320, 150], 6, 2, "UPrime2Mean") // xz
t, u = av.get_time_series(out, [1680, 1320, 150], 6, 3, "UPrime2Mean") // yy

----------------------------------------------------------------------------------------------------------------------------------------------------------

awcmviewer can read the entire flow over a single time step and extract fields values on a line connecting pointA to pointB. This can be done from vdi node as:

mpirun -np 16 /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis line -fields 'U UPrime2Mean UMean' -time 3600 -out fowf15mwR0 -pwd -pointA '10 1320 0' -pointB '10 1320 960'

It would be more efficient to submit a batch job for extracting data over multiple lines as:

srun /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis line -fields 'U UPrime2Mean UMean' -time 3600 -out fowf15mwR0 -pwd -pointA '10 1320 0' -pointB '10 1320 960' > line1.out

srun /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis line -fields 'U UPrime2Mean UMean' -time 3600 -out fowf15mwR0 -pwd -pointA '90 1320 0' -pointB '90 1320 960' > line2.out

add more lines as needed

----------------------------------------------------------------------------------------------------------------------------------------------------------

extract flow field on a 2D plane and save as xdmf + h5 

srun /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis line -fields 'U UPrime2Mean UMean Q Vort' -time 3600 -out fowf15mwR0 -pwd -pointA '90 1320 150' -pointB '0 0 1' > slice.out

here, pointA is a reference point and pointB is a unit vector, which defines a plane over which the flow is extracted. 

Reconstruct 3D flow field and save as xdmf + h5 fomat:

srun /project/def-alamj/shared/bin/v2306/awcmviewer -INP post_process.inp -analysis xdmf -fields 'U UPrime2Mean UMean Q Vort' -time 3600 -out fowf15mwR0 -pwd  > xdmf.out









