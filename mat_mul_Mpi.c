/*
    Plan:
        Take input values from the user for both matrices.
        Scatter the first matrix amoung the processes and broad cast
        the second matrix to the all of the processes. Each process 
        calculates its own local submatrix and send it to process 0.
        Process 0 collects the local submatrices and and combine them
        in a matrix and displays it.
*/


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>

void get_size(int my_rank, int* rows, int* columns){
    if(my_rank == 0){
        
        printf("Rows : ");
        fflush(stdout);
        scanf("%d",rows);

        printf("Columns : ");
        fflush(stdout);
        scanf("%d",columns);

    }
    MPI_Bcast(rows,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(columns,1,MPI_INT,0,MPI_COMM_WORLD);
}

int check_dim(int my_rank,int cols, int rows){
    if(my_rank == 0){
        if(cols != rows){
            printf("\nDimension Error!!!!\n");
            printf("Can't multiply these two matrices.\n");
            return 1;
        }
        return 0;
    }
    return 0;
}

void get_elements(int my_rank, int r1, int c1, double mat_a[], int local_row, int counts[], int displs[],int r2, int c2, double mat_b[]){
    double* matrix1 = NULL;
    if(my_rank == 0){

        matrix1 = malloc(r1 * c1 * sizeof(double));

        printf("Enter the elements of First Matrix: \n");
        fflush(stdout);

        for(int i=0 ;i< r1;i++){
            for(int j=0;j<c1;j++){
                scanf("%lf",&matrix1[i*c1 + j]);
            }
        }

        printf("Enter the elements of Second Matrix: \n");
        fflush(stdout);
        for(int i=0 ;i< r2;i++){
            for(int j=0;j<c2;j++){
                scanf("%lf",&mat_b[i*c2 + j]);
            }
        }
    }

    MPI_Scatterv(matrix1,counts,displs,MPI_DOUBLE,mat_a,local_row,MPI_DOUBLE,0,MPI_COMM_WORLD); // Scatter the First matrix amoung the processes.
    MPI_Bcast(mat_b,r2*c2,MPI_DOUBLE,0,MPI_COMM_WORLD); // Broad cast the second matrix to all the processes in the communicator.

    if (my_rank == 0) free(matrix1); // Free the memory of the matrix that was allocated before.

}

int main(){

    int my_rank, comm_size;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int row1, col1, row2, col2;

    if(my_rank == 0){
        printf("Enter the values for First Matrix: \n");
    }
    get_size(my_rank,&row1,&col1);  // Gets the number of rows and columns for the first matirx
    
    if(my_rank == 0){
        printf("Enter the values for Second Matrix: \n");
    }
    get_size(my_rank,&row2,&col2);  // Gets the number of rows and columns for the second matirx

    // // This code is for debugging
    // printf("Process %d, r1 : %d, c1 %d, r2: %d, c2: %d\n",my_rank,row1,col1,row2,col2);
    // fflush(stdout);
  
    if(check_dim(my_rank,col1,row2) == 1){              // Check the dimensions for multiplication              
        MPI_Finalize();                                 // If not valid then program exits with a warning!!!!
        return 0;
    }

    int base_rows = row1 / comm_size;
    int remainder = row1 % comm_size;
    int local_rows = base_rows;

    if(my_rank == comm_size - 1){
        local_rows += remainder;
    }

    double matrix_A[local_rows * col1], matrix_B[row2 * col2];
    double local_r[local_rows * col2], global_r[row1 * col2];

    int counts[comm_size]; // number of elements each process will send
    int displs[comm_size]; // Displacements for each process
    int current_displ = 0; 

    for(int i = 0; i < comm_size; i++){
        int current_rows = base_rows; // Number of rows each process will get
        if(i == comm_size - 1){
            // If number rows is not completely divisible by the number of processes
            // then the extra rows will be assigned to the last process of the 
            // comunicator
            current_rows += remainder;
        }

        counts[i] = current_rows * col1; // Element count to be sent by each process
        displs[i] = current_displ; // Displacement for each element sent by the process
        current_displ += counts[i]; // Displcement for the next element
    }
    
    // Each process will get its local values for the matrix and vector
    get_elements(my_rank, row1, col1, matrix_A, counts[my_rank], counts, displs, row2, col2, matrix_B);


    // For testing. Printing the matrices from each process
    // if(my_rank == comm_size -1){
    //     printf("First Local Matrix from process %d:\n",my_rank);
    //     fflush(stdout);
    //     for(int i=0;i<local_rows;i++){
    //         for(int j=0;j<col1;j++){
    //             printf("%f ",matrix_A[i*col1 + j]);
    //             fflush(stdout);
    //         }
    //         printf("\n");
    //         fflush(stdout);
    //     }

    //     printf("\nSecond Matrix from process %d:\n",my_rank);
    //     fflush(stdout);
    //     for(int i=0;i<row2;i++){
    //         for(int j=0;j<col2;j++){
    //             printf("%f ",matrix_B[i*col2 + j]);
    //             fflush(stdout);
    //         }
    //         printf("\n");
    //         fflush(stdout);
    //     }
    // }


    // Below is the matrix multiplication logic
    for(int i=0;i<local_rows;i++){
        for(int j=0;j<col2;j++){
            local_r[i * col2 + j] = 0.0;
            for(int k=0;k<col1;k++){
                local_r[i * col2 + j] += matrix_A[i * col1 + k] * matrix_B[k * col2 + j];
            }
        }
    }


    // Below code was used for testing the correctness multiplication logic
    // if (my_rank == comm_size-1) {
    //     printf("\nFinal Result Matrix\n");
    //     fflush(stdout);
    //     for (int i = 0; i < row1; i++) {
    //         for (int j = 0; j < col2; j++) {
    //             printf("%8.2f ", local_r[i * col2 + j]);
    //         }
    //         printf("\n");
    //     }
    //     fflush(stdout);
    // }


    // Now we need to send the local resultant matrix to process 0
    int gather_counts[comm_size];
    int gather_displs[comm_size];
    int current_gather_displ = 0;

    for(int i = 0; i < comm_size; i++){
        int current_rows = base_rows; 
        if(i == comm_size - 1){
            current_rows += remainder;
        }
        gather_counts[i] = current_rows * col2; 
        gather_displs[i] = current_gather_displ; 
        current_gather_displ += gather_counts[i]; 
    }

    MPI_Gatherv(local_r, gather_counts[my_rank], MPI_DOUBLE, 
                global_r, gather_counts, gather_displs, MPI_DOUBLE, 
                0, MPI_COMM_WORLD);


    if (my_rank == 0) {
        printf("\n--- Final Result Matrix ---\n");
        fflush(stdout);
        for (int i = 0; i < row1; i++) {
            for (int j = 0; j < col2; j++) {
                printf("%8.2f ", global_r[i * col2 + j]);
            }
            printf("\n");
        }
        printf("---------------------------\n");
        fflush(stdout);
    }

    MPI_Finalize();
    return 0;
}
