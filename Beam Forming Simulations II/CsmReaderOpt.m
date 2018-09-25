% Compatible with HDF5 CSM definitions Rev 2.2

h5in_file='Benchmark0_optional.h5';

conventionalSolution = h5read(h5in_file,'/GridSolution/conventionalSolution');
gridPointCoordinatesM = h5read(h5in_file,'/GridSolution/gridPointCoordinatesM');
steeringArrayImaginary = h5read(h5in_file,'/GridSolution/steeringArrayImaginary');
steeringArrayReal = h5read(h5in_file,'/GridSolution/steeringArrayReal');