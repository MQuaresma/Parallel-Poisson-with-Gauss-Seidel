#include <mpi.h>
#include "utils.h"

void updateValues(int rows_per_proc, float tempProc[2][(rows_per_proc+2)*N], int last, int rank, int no_procs, MPI_Status status){
    if(!rank){
        //comunica com o ultimo
        MPI_Recv( tempProc[last], 1*N, MPI_FLOAT, no_procs-1, 2, MPI_COMM_WORLD, &status );
        MPI_Send( &(tempProc[last][1*N]), 1*N, MPI_FLOAT, no_procs-1, 2, MPI_COMM_WORLD);
    }else{
        if(rank!=1)
            MPI_Send( &(tempProc[last][1*N]), 1*N, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD);

        if(rank!=no_procs-1)
            MPI_Recv( &(tempProc[last][(rows_per_proc+1)*N]), 1*N, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &status );

        if(rank!=no_procs-1)
            MPI_Send( &(tempProc[last][rows_per_proc*N]), 1*N, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD);

        if(rank!=1)
            MPI_Recv( tempProc[last], 1*N, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &status );

        if(rank==no_procs-1){
            //comunica com o 0
            MPI_Send( &(tempProc[last][1*N]), 1*N, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
            MPI_Recv( tempProc[last], 1*N, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status );
        }
    }

}


float phase(int rows_per_proc, float tempProc[2][(rows_per_proc+2)*N], int rank, int last, int no_procs, MPI_Status status){
    float dif, temp;
    dif = 0.0f;

    //Black
    for(int i = 1; i < rows_per_proc-1; i ++)
        for(int j = 2 - i%2; j < N-1; j +=2){
            tempProc[!last][i*N+j] = (tempProc[last][(i-1)*N+j] + tempProc[last][i*N+j-1] + tempProc[last][i*N+j+1] + tempProc[last][(i+1)*N+j]) / 4.0f;
            temp = fabs(tempProc[!last][i*N+j] - tempProc[last][i*N+j]);
            if(temp > dif) dif = temp;
        }

    updateValues(rows_per_proc, tempProc, !last, rank, no_procs, status);

    //Red
    for(int i = 1; i < rows_per_proc-1; i ++)
        for(int j = 1 + i%2; j < N-1; j +=2){
            tempProc[!last][i*N+j] = (tempProc[!last][(i-1)*N+j] + tempProc[!last][i*N+j-1] + tempProc[!last][i*N+j+1] + tempProc[!last][(i+1)*N+j]) / 4.0f;
            temp = fabs(tempProc[!last][i*N+j] - tempProc[last][i*N+j]);
            if(temp > dif) dif = temp;
        }

    updateValues(rows_per_proc, tempProc, !last, rank, no_procs, status);

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
    int it = 0, last=0, rows_per_proc, remaining_rows, work_load = N-1, begin;
    float dif = tol + 1.0f; 

    rows_per_proc = work_load/(no_procs-1);      //truncates result
    remaining_rows = work_load % (no_procs-1);   //remaining rows

    float tempProc[2][( ((rank==0) ? remaining_rows : rows_per_proc) +2)*N];   

    //send data
    if(rank == 0){
        begin = 0;
        for(int i=1; i<no_procs; i++, begin+=rows_per_proc) {
            MPI_Send(plate[last]+begin*N*sizeof(float), (rows_per_proc+2)*N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
        }
    }else{
        MPI_Recv( tempProc[last], (rows_per_proc + 2)*N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status );

        for(int i=0; i<rows_per_proc+2; i++)
            for(int j=0; j<N; j++)
                tempProc[!last][i*N+j] = tempProc[last][i*N+j];
    }

    while(dif > tol){
        if(!rank)
            dif = phase(remaining_rows,tempProc,rank,last, no_procs, status);
        else
            dif = phase(rows_per_proc,tempProc,rank,last, no_procs, status);
        
        it ++;
        last = !last; //update matrix
    }

    //receive data
    if( rank == 0 ){
        begin = 1;
        for(int i=1; i<no_procs; i++, begin+=rows_per_proc){
            MPI_Recv( plate[begin], (rows_per_proc)*N, MPI_FLOAT, i, 3, MPI_COMM_WORLD, &status);
        }
    }else{
        MPI_Send( &(tempProc[last][1*N]), rows_per_proc*N, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
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
