#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void relax(double*, double*, int, int);
int getPrecision(double*, double*, int, double);
double * init_array(double*, int);
void print_arrays(double*, double*, int);
void print_matrices(double*, double*, int);
double * read_array(int n);

double *a1;
double *a2;

int main(int argc, char **argv) {

  //Init vars
  int rc, rank, nprocs, n;
  double p, t_start, t_end;
  n = atoi(argv[1]); //Matrix dimensions
  p = atof(argv[2]); //Precision
  rc = MPI_Init(&argc, &argv); //Init MPI system

  //Error check
  if (rc != MPI_SUCCESS) {//If MPI INIT unsuccessful
    printf ("Error sending MPI program\n"); //Print to console
    MPI_Abort(MPI_COMM_WORLD, rc); //Kill program w/ error code
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);//Get rank of current processor
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);//Get number of processors

  MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */

  t_start = MPI_Wtime();

  //Init matrix
//  a1 = init_array(a1, n); //Create populated matrix
//  a2 = init_array(a2, n); //Create empty matrix with borders

  a1 = read_array(n);
  a2 = read_array(n);

  //Split work
  //Based on feedback from Shared Implementation:
  // -> Assign entire rows to communicators to reduce complexity
  // -> Compared to calculating exact number of cells and assigning data sets which span between rows

  int nrows = (n-2); //Number of computable rows (Excludes borders)
  int quo = nrows/nprocs; //Number of rows per thread
  int rem = nrows % nprocs;


  //Main thread
  if (rank==0) {
    //4 Processors
    printf("MPI Communicator initialised with %d processors", nprocs);
    printf("\n");

    //6 total rows, 4 computational rows => 1 row per processor
    printf("Starting computation on %d total rows, with %d rows per thread and a remainder of %d", nrows, quo, rem);
    printf("\n");

  }

  //Allocate memory for arrays to hold sending and receiving work indices
  int *shift_send = malloc(nprocs*sizeof(int));
  int *shift_recv = malloc(nprocs*sizeof(int));
  int *ncells_send = malloc(nprocs*sizeof(int));
  int *ncells_recv = malloc(nprocs*sizeof(int));

  //Loop through processors in communicator:
    for (int i = 0; i < nprocs; i++){
      //Get index for each row:
      //-> Multiply loop index by quotient & cells per row to get index
      //E.g. 0*2*10 = 0, 1*2*10 = 20, 2*2*10 = 40 => 20 cells per processor
      int rowIndices = i * quo * n;
      shift_send[i] =  rowIndices; //Assign send index for each row to array -> Scatter
      shift_recv[i] = n + rowIndices; //Assign retrieval index for each row to array -> Gather

      //Special case for last processor (n-1):
      if (i == nprocs-1) {
        //Get send index of unallocated cells:
        // -> Remainder given by number of rows MOD number of processor
        // -> Add cells per row (Includes borders) to remainder, and multiply by length
        ncells_send[i] = (quo + rem + 2) * n; //Includes borders
        ncells_recv[i] = (quo + rem) * n; //Excludes border
      } else {
        //Determine number of cells to distribute to each processor
  			ncells_send[i] = (quo + 2) * n; //Includes border
  			ncells_recv[i] = quo * n; //Excludes border
      }
    }//End Loop


  //Allocate memory for recv and send buffers
  //Once again: final processor will have a different number of cells
  //This is calculated using rem like before
  int nrecv, nsend;

  //Calculate size of buffer for each processor in communicator
	if (rank == (nprocs-1)) {//Special case
    //Final process
    // -> Remaining: unallocated cells
    nrecv = (quo + rem + 2) * n; //Includes borders
    nsend = (quo + rem) * n; //Excludes border
  } else {
    //All processes besides final
    nrecv = (quo + 2) * n; //Includes border
    nsend = quo * n; //Excludes border
  }



  //Now we can allocate buffer array using sizes calculated above
  double *buffer_recv = malloc(nrecv * sizeof(double));
	double *buffer_send = malloc(nsend * sizeof(double));

  int precise = 0;
  int iterations = 0;

  //Main program Loop
  //Repeates until precision threshold is met
  if (rank == 0) {
    //printf("\n!!!!!");
  //  print_arrays(a1, a2, n);
    //printf("!!!!!\n");
  }

  while (1) { //Loop until told otherwise
    //scatter data
    /*
    MPI_Scatterv:
    Scatters data of varying sizes to processes in Communicator
    Different to MPI_Scatter, which can only send data of equal sizes

    @param shift_send: Displacement applied to message

    */
    MPI_Scatterv(a1, ncells_send, shift_send, MPI_DOUBLE, buffer_recv,
            nrecv, MPI_DOUBLE, 0, MPI_COMM_WORLD);



    //Execute relax function with scattered data
    relax(buffer_recv, buffer_send, n, nrecv);



    //Retrieve scattered & relaxed data
    /*
    MPI_Gatherv:
    Scatters data of varying sizes to processes in Communicator
    Different to MPI_Gatherv, which can only receive data of equal sizes

    @param nsend: Number of items to gather
    */
    MPI_Gatherv(buffer_send, nsend, MPI_DOUBLE, a2, ncells_recv,
          shift_recv, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // if (rank == 0) {
    //   printf("\n!!!!!\n");
    //   print_matrices(a1, a2, n);
    //   printf("!!!!!\n");
    //   break;
    // }

    //Root process checks precision
    //-> Root process has rank 0
    if (rank == 0) {
      //printf("checking");
      precise = getPrecision(a1, a2, n, p);

      //printf("%i",precise);
      //Precise:
      //-> 1: Matrices are within precision
      //-> 0: Matrices not within precision yet
    }

    // Broadcast single value: precision- to entire communicator
    MPI_Bcast(&precise, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //End of iteration
    //Swap matrices
    double *temp = a1; //Temp acts as a buffer
    a1 = a2; //Copy a2 to a1
    a2 = temp; //Copy temp(a1) to a2

    iterations++;




    //Determine break condition:
    //If precision reached
		if (precise){
      //printf("processor with rank: %d reported precision", rank);
			break;//Stop loop
		}

  }

  //End of main program loop
  if (rank == 0) {
  //print_arrays(a1, a2, n);
  }

  MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
  t_end = MPI_Wtime();

  //Root process summary
  if (rank == 0) {
    printf("Completed in %d iterations", iterations);
    printf("\n");

    //printf("Runtime = %f\n", t_end-t_start);
    double t = t_end - t_start;
    printf("Duration:");
    printf("%f",t);
    printf("\n");
  }

  //Destroy MPI Communicator:
  MPI_Finalize();

  //Deallocate Memory
	free(a1);
	free(a2);
  //free(temp);
	free(shift_send);
	free(shift_recv);
	free(buffer_recv);
	free(buffer_send);
  free(ncells_send);
  free(ncells_recv);
}

//RELAXATION FUNCTION:

//RELAX FUNCTION
void relax(double* input, double* output, int nrow, int nrecv) {

  int j = nrow; //Do not relax first row
  int i = nrecv - nrow - 1; //Do not relax final row
  int c = 0; //counter
  float cross;
  //printf("\nThread: %d:Starting at [%d,%d] for %d \n", id, i, j, c);
  for (int k = j; k <= i; k++) {//Loop through provided section
    //Only operate on cells inside the borders
    //Border detection:
    // -> Calculate remainder after division of consecutive indices by row length
    // -> If neither divide equally, index k is not an edge case
    if ((k % nrow != 0) && ((k+1) % nrow != 0)) {
      //Sum: Add vertical and horizontal neighbours together
      cross = input[k-1] + input[k+1] + input[k+nrow] + input[k-nrow];
      //Average: Divide result by 4
      output[c] = (cross/4);
    } else { //Edge cases
      output[c] = input[k]; //Do nothing: blindly copy the input data
    }
    c++; //Increment counter
  }
}

//HELPER METHODS-------------------------------------------

//Returns whether new matrix is within precision of previous matrix
int getPrecision(double *a1, double *a2, int n, double precision){

  //Loop through matrix
    for (int i=0; i<n*n; i++) {
      //Calculate numerical difference and compare with given precision
  		if (fabs(a1[i]-a2[i]) > precision){
  			return 0;//If outside precision, return 0/False
  		}
  	}
  //If none were outside precision
	return 1;//Return 1/True
}



// //Function to initialise arrays with starting values
// double * init_array(double *arr, int n){
//
//
//   //Allocate memory for array
//   arr = malloc(n*n*sizeof(double));
//
//   //Loop through indices
// 	for (int i = 0; i < (n*n); i++){
//     //Check if index lies on borders
//     if (i < n || ((i) % n == 0) || ((i+1) % n == 0) || (i>= (n*n)-n)) {
//       arr[i] = 1;//Assign 1 if so
//     } else {
//       arr[i] = 0;//Otherwise assign 0
//     }
// 	}
//   return arr;
// }



double * read_array(int n) {
  FILE *datafile;
  datafile = fopen("textfile.txt", "r");

    //read file into array
    double *data_array = malloc(n*n*sizeof(double));

    int i;

    for (i = 0; i < n*n; i++)
    {
      fscanf(datafile, "%lf\n", &data_array[i]);
    }


    return data_array;

}





//PRETTY PRINT MATRIX
void print_arrays(double * a1, double * a2, int n) {
  printf("------------------------------------MATRIX 1------------------------------------");
  printf("\n");
  //Print with nice lines and shapes :)
  for (int i = 0; i < n*n; i++) {
    if (i % n == 0){
			printf("\n");
		}
		printf("%f ",a1[i]);
	}
  printf("\n");
  printf("\n");


printf("----------------------------------------MATRIX 2---------------------------------------");

    for (int i = 0; i < n*n; i++) {
      if (i % n == 0){
        printf("\n");
      }
      printf("%f ",a2[i]);
    }

    printf("\n");
  }
