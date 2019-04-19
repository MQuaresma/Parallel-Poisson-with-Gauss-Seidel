#include "../UTILS/utils.h"

/* 
 * Returns the number of iterations performed
 */
int poissongs(float plate[][N][N], float tol){
    int it = 0, last=0;
    float dif = tol + 1.0f, temp;

    while(dif > tol){
        dif = 0.0f;             //differences at the borders are 0.0f
        for(int i = 1; i < N-1; i ++)
            for(int j = 1; j < N-1; j ++){
                plate[!last][i][j] = (plate[!last][i-1][j] + plate[!last][i][j-1] + plate[last][i][j+1] + plate[last][i+1][j]) / 4.0f;
                temp = fabs(plate[!last][i][j] - plate[last][i][j]);
                if(temp > dif) dif = temp;
            }
        it ++;
        last = !last; //update matrix
        printf("%f\n", dif);
        printf("%d\n", last);
    }
    return it;
}


int main(int argc, char *argv[]){
    float plate[2][N][N];
    float tol = 1.0f / (N*N);
    int it;

    initPlate(plate);

    it = poissongs(plate, tol);
    printf("Sequential Poisson GS\n Iteration Count: %d\n", it);

    return 0;
}
