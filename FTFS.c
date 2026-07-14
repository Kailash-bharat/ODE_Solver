#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.141592653589793

double initial_condition(double x){
    return sin(2 * PI * x);  
}

double exact_solution(double x, double t, double c){
    return sin(2 * PI * (x - t));
}

void write_to_file(char *filename, double x[], double u[], int nx){
    FILE *fp = fopen(filename, "w");
    for (int i = 0; i < nx; i++){
        fprintf(fp, "%lf %lf\n", x[i], u[i]);
    }
    fclose(fp);
}

int main(){

    double c = 1.0;
    double L = 1.0;
    double dx = 0.01;
    double v = 0.72; //Use 0.7 for sake of Problem 3 later
    double dt = v * dx / c;
    double T = 1.0;
    int NX = 100;
    int NT = (int)(T/dt);

    double x[NX];
    double u_ftfs[NX], u_ftbs[NX];
    double u_new[NX];

    for (int i = 0; i < NX; i++){
        x[i] = i * dx;
        u_ftfs[i] = initial_condition(x[i]);
        u_ftbs[i] = initial_condition(x[i]);
    }

    write_to_file("ftfs_0.dat", x, u_ftfs, NX);
    write_to_file("ftbs_0.dat", x, u_ftbs, NX);

    for (int n = 0; n < NT; n++){

        // FTFS
        for (int i = 0; i < NX - 1; i++){
            u_new[i] = u_ftfs[i] - v * (u_ftfs[i + 1] - u_ftfs[i]);
        }
        u_new[NX - 1] = u_new[0];

        for (int i = 0; i < NX; i++)
            u_ftfs[i] = u_new[i];

        // FTBS
        for (int i = 1; i < NX; i++){
            u_new[i] = u_ftbs[i] - v * (u_ftbs[i] - u_ftbs[i - 1]);
        }
        u_new[0] = u_new[NX - 1];

        for (int i = 0; i < NX; i++)
            u_ftbs[i] = u_new[i];

        if(n*dt > 0.25 - dt/2 && n*dt <= 0.25 + dt/2){
            write_to_file("ftfs_0.25.dat", x, u_ftfs, NX);
            write_to_file("ftbs_0.25.dat", x, u_ftbs, NX);
        }
        else if(n*dt > 0.5 - dt/2 && n*dt <= 0.5 + dt/2){
            write_to_file("ftfs_0.5.dat", x, u_ftfs, NX);
            write_to_file("ftbs_0.5.dat", x, u_ftbs, NX);
        }
        else if(n*dt > 0.75 - dt/2 && n*dt <= 0.75 + dt/2){
            write_to_file("ftfs_0.75.dat", x, u_ftfs, NX);
            write_to_file("ftbs_0.75.dat", x, u_ftbs, NX);
        }
        else if(n == NT-1){
            write_to_file("ftfs_1.0.dat", x, u_ftfs, NX);
            write_to_file("ftbs_1.0.dat", x, u_ftbs, NX);
        }
    }

    double error_ftfs = 0.0, error_ftbs = 0.0;

    for (int i = 0; i < NX; i++){

        double u_exact = exact_solution(x[i], T, c);

        error_ftfs += pow(u_ftfs[i] - u_exact, 2);
        error_ftbs += pow(u_ftbs[i] - u_exact, 2);
    }

    error_ftfs = sqrt(error_ftfs * dx);
    error_ftbs = sqrt(error_ftbs * dx);

    printf("L2 Error (FTFS) at T = %lf\n", error_ftfs);
    printf("L2 Error (FTBS) at T = %lf\n", error_ftbs);

    return 0;
}
