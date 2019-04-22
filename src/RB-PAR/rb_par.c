#include <mpi.h>
#include "utils.h"

void updateValues(float plate[][N*N], int last, int offset, int rows, int rank, int no_procs, MPI_Status status){
    if(!rank){
        MPI_Recv( &(plate[last][rows*N]), N, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &status );
        MPI_Send( &(plate[last][(rows-1)*N]), N, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD);
    }else{
        MPI_Send( &(plate[last][offset*N]), N, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD);

        if(rank!=no_procs-1)
            MPI_Recv( &(plate[last][offset*N + rows*N]), N, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &status );

        if(rank!=no_procs-1)
            MPI_Send( &(plate[last][offset*N + (rows-1)*N]), N, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD);

        MPI_Recv( &(plate[last][(offset-1)*N]), N, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &status );
    }
}

float phase(float plate[][N*N], int last, int offset, int rows, int rank, int no_procs, MPI_Status status){
    float dif, temp;
    dif = 0.0f;

    //Black
    for(int i = (!rank) ? 1 : offset; i < offset + ((rank==no_procs-1) ? rows-1 : rows); i ++)
        for(int j = 2 - i%2; j < N-1; j +=2){
            plate[!last][i*N+j] = (plate[last][(i-1)*N+j] + plate[last][i*N+j-1] + plate[last][i*N+j+1] + plate[last][(i+1)*N+j]) / 4.0f;
            temp = fabs(plate[!last][i*N+j] - plate[last][i*N+j]);
            if(temp > dif) dif = temp;
        }

    updateValues(plate, !last, offset, rows, rank, no_procs, status);

    //Red
    for(int i = (!rank) ? 1 : offset; i < offset + ((rank==no_procs-1) ? rows-1 : rows); i ++)
        for(int j = 1 + i%2; j < N-1; j +=2){
            plate[!last][i*N+j] = (plate[!last][(i-1)*N+j] + plate[!last][i*N+j-1] + plate[!last][i*N+j+1] + plate[!last][(i+1)*N+j]) / 4.0f;
            temp = fabs(plate[!last][i*N+j] - plate[last][i*N+j]);
            if(temp > dif) dif = temp;
        }

    updateValues(plate, !last, offset, rows, rank, no_procs, status);

    if(rank){
        MPI_Send( &dif, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        MPI_Recv( &dif, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
    }else{
        for(int i=1; i<no_procs; i++){
            MPI_Recv( &temp, 1, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &status);
            if(temp > dif) dif = temp;
        }
        for(int i=1; i<no_procs; i++){
            MPI_Send( &dif, 1, MPI_FLOAT, i, 1, MPI_COMM_WORLD);
        }
    }
    return dif;
}

/* 
 * Returns the number of iterations performed
 */
int poissongs(float plate[][N*N], float tol, int rank, int no_procs, MPI_Status status){
    int it = 0, last=0, rows_per_proc, remaining_rows, work_load = N-2, begin, offset;
    float dif = tol + 1.0f; 

    rows_per_proc = work_load/(no_procs-1);      //truncates result
    remaining_rows = work_load % (no_procs-1);   //remaining rows

    offset = rank * rows_per_proc;

    while(dif > tol){
        if(rank==no_procs-1)
            dif = phase(plate, last, offset, remaining_rows, rank, no_procs, status);
        else
            dif = phase(plate, last, offset, rows_per_proc, rank, no_procs, status);
        
        it ++;
        last = !last; //update matrix
    }

    //collect data
    if(!rank){
        begin = rows_per_proc;

        for(int i=1; i<no_procs-1; i++, begin+=rows_per_proc)
            MPI_Recv(&(plate[last][begin]), rows_per_proc*N, MPI_FLOAT, i, 3, MPI_COMM_WORLD, &status);

        MPI_Recv(&(plate[last][begin]), remaining_rows*N, MPI_FLOAT, no_procs-1, 3, MPI_COMM_WORLD, &status);
    }else if(rank == no_procs-1) {
        MPI_Send( &(plate[last][offset]), remaining_rows*N, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
    }else{
        MPI_Send( &(plate[last][offset]), rows_per_proc*N, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
    }

    return it;
}


int main(int argc, char *argv[]){
    float plate[2][N*N];
    float tol = 1.0f / (N*N);
    int it, rank, no_procs;

    initPlateForMPI(plate);

    //Move MPI process creation to poissongs
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &no_procs);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    it = poissongs(plate, tol, rank, no_procs, status);

    MPI_Finalize();

    if(rank == 0)
        printf("Sequential Poisson GS with Red-Back strategy\n Iteration Count: %d\n", it);

    return 0;
}
