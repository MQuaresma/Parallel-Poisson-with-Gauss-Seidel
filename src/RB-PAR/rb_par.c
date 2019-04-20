#include "utils.h"
#include <mpi.h>

void phase(float **tempProc, int rank, float temp, int last, int ph){
    if (rank == 0) {

        //TODO: Processar as remaining lines
         
        //comunica com o ultimo
        MPI_Recv( tempProc[last], 1*N, MPI_FLOAT, no_procs-1, 2, MPI_COMM_WORLD, &status );
        MPI_Send( &(tempProc[last][1*N]), 1*N, MPI_FLOAT, no_procs-1, 2, MPI_COMM_WORLD);
    }else{
        MPI_Recv( tempProc[last], (rows_per_proc + 2)*N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status );
        copy(tempProc[last],0,tempProc[!last],0,rows_per_proc+2);

        if(ph==0){
            //Black
            for(int i = 1; i < N-1; i ++)
                for(int j = 2 - i%2; j < N-1; j +=2)
                    tempProc[!last][i][j] = (tempProc[last][i-1][j] + tempProc[last][i][j-1] + tempProc[last][i][j+1] + tempProc[last][i+1][j]) / 4.0f;
        else{
            //Red
            for(int i = 1; i < N-1; i ++)
                for(int j = 1 + i%2; j < N-1; j +=2)
                    plate[!last][i][j] = (plate[!last][i-1][j] + plate[!last][i][j-1] + plate[!last][i][j+1] + plate[!last][i+1][j]) / 4.0f;
        }

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

/* 
 * Returns the number of iterations performed
 */
int poissongs(float plate[][N][N], float tol, int rank, int no_procs){
    int it = 0, last=0, rows_per_proc, remaining_rows, work_load = N-1, begin;
    float dif = tol + 1.0f; 

    rows_per_proc = work_load/(no_procs-1);      //truncates result -> by defect
    remaining_rows = work_load % (no_procs-1);   //remaining rows

    float tempProc[2][(rows_per_proc+2)*N];

    //send data
    if(rank == 0){
        begin = 0;
        for(int i=1; i<no_procs; i++, begin+=rows_per_proc) {
            MPI_Send(plate[begin], (rows_per_proc+2)*N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
        }
    }

    while(dif > tol){
        dif = 0.0f;             //differences at the borders are 0.0f

        //Black
        phase(tempProc,rank,temp,last,0);
        
        //TODO: wait
        
        //Red
        phase(tempProc,rank,temp,last,1);

        //TODO: process 0 calculate dif, and send to others processes(maybe)

        it ++;
        last = !last; //update matrix
    }

    //receive data
    if( rank == 0 ){
        begin = 1;
        for(int i=1; i<no_procs; i++, begin+=rows_per_proc){
            MPI_Recv( plate[begin], (rows_per_proc)*N, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &status);
        }
    }else{
        MPI_Send( &(tempProc[last][1*N]), rows_per_proc*N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
    }

    return it;
}


int main(int argc, char *argv[]){
    float plate[2][N][N];
    float tol = 1.0f / (N*N);
    int it, rank, no_procs,

    initPlate(plate);

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &no_procs);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    it = poissongs(plate, tol, rank, no_procs);

    MPI_Finalize();

    printf("Sequential Poisson GS with Red-Back strategy\n Iteration Count: %d\n", it);

    return 0;
}
