
 MPP Coursework: Reverse edge-detection algorithm
 ------------------------------------------------------------------------------

 1. Application building commands (mpiversion by defaut) :

 $ make serialversion   : Compiles and links the serial implementation
                          of the application in binserial directory.

 
 $ make mpiversion      : Compiles and links the MPI implementation
                          of the application in bin directory.

 Note: the makefile is configured to compile the code using the PGI
 compiler on Morar, with the default MPI library (mpif90).


 2. Execute the appliaction
 
 For the MPI version:
 $ ./bin/invertedges  INPUT_FILE NUM_ITERATIONS IT_BTW_RED [OUTPUT_FILE]
 For the serial version:
 $ ./binserial/invertedges  INPUT_FILE NUM_ITERATIONS IT_BTW_RED [OUTPUT_FILE]
 
 where:
   - INPUT_FILE: Is the input pgm image
   - NUM_ITERATIONS: Number of iterations the algorithm will be executed, if
                     it is set to 0, it will use the stopping criterion with
                     a 0.1 threshold of maximum pixel change.
   - IT_BTW_RED: Number of iterations performed between the reduction operations
                 needed to print out the pixel average or calculate the stopping
                 criterion value.
   - OUTPUT_FILE: Optional argument to specify the output image path, if it is
                  not specified, the default value is 'output.pgm'


 3. To execute the correctness test:
 $ make test
 or
 $ ./tests/correctness/ctest.sh BINARY_FILE INPUT_FOLDER PROPER_RESULTS

 4. To execute the performance tests:
 $ ./tests/performance/reductions.sh BINARY_FILE INPUT_FILE NUM_THREADS
 $ ./tests/performance/image_size.sh BINARY_FILE NUM_THREADS INPUT_FOLDER
 $ ./tests/performance/scalability.sh BINARY_FILE INPUT_FILE

 Note: more information about the arguments available in the scripts.
 Note: Some examples in mpibatch.sge
