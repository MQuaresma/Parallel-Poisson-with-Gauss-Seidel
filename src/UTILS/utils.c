#include "utils.h"

/**
 * Initialize heat plate with temperature values
 */
void initPlate(float plate[][N][N]){
    for(int i = 0; i < N; i ++){
        plate[0][0][i] = 100.0f;
        plate[0][i][0] = 100.0f;
        plate[0][i][N-1] = 100.0f;
        plate[0][N-1][i] = 0.0f;

        plate[1][0][i] = 100.0f;
        plate[1][i][0] = 100.0f;
        plate[1][i][N-1] = 100.0f;
        plate[1][N-1][i] = 0.0f;
    }

    for(int i = 1; i < N-1; i ++)
        for(int j = 1; j < N-1; j ++)
            plate[0][i][j] = 50.0f;

    plate[0][N-1][0] = 0.0f;
    plate[1][N-1][0] = 0.0f;
}

void initPlateForMPI(float plate[][N*N]){
    for(int i = 0; i < N; i ++){
        plate[0][0*N+i] = 100.0f;
        plate[0][i*N+0] = 100.0f;
        plate[0][i*N+N-1] = 100.0f;
        plate[0][(N-1)*N+i] = 0.0f;

        plate[1][0*N+i] = 100.0f;
        plate[1][i*N+0] = 100.0f;
        plate[1][i*N+N-1] = 100.0f;
        plate[1][(N-1)*N+i] = 0.0f;
    }

    for(int i = 1; i < N-1; i ++)
        for(int j = 1; j < N-1; j ++)
            plate[0][i*N+j] = 50.0f;

    plate[0][(N-1)*N+0] = 0.0f;
    plate[1][(N-1)*N+0] = 0.0f;
}
