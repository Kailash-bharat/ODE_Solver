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
    double v = 1.58;
    double dt = v * dx / c;
    double T = 1.0;

    int NX = 100;
    int NT = (int)(T / dt);

    double x[NX];
    double u[NX], u_new[NX];

    for (int i = 0; i < NX; i++){
        x[i] = i * dx;
        u[i] = initial_condition(x[i]);
    }

    char filename[50];
    sprintf(filename, "lw_0_v_%.2f.dat", v);
    write_to_file(filename, x, u, NX);

    for (int n = 0; n < NT; n++){

        for (int i = 1; i < NX - 1; i++){
            u_new[i] = u[i]
                     - 0.5 * v * (u[i+1] - u[i-1])
                     + 0.5 * v * v * (u[i+1] - 2*u[i] + u[i-1]);
        }

        u_new[0] = u[0]
                 - 0.5 * v * (u[1] - u[NX-1])
                 + 0.5 * v * v * (u[1] - 2*u[0] + u[NX-1]);

        u_new[NX-1] = u[NX-1]
                    - 0.5 * v * (u[0] - u[NX-2])
                    + 0.5 * v * v * (u[0] - 2*u[NX-1] + u[NX-2]);

        for (int i = 0; i < NX; i++)
            u[i] = u_new[i];

        if(n*dt > 0.25 - dt/2 && n*dt <= 0.25 + dt/2){
            sprintf(filename, "lw_0.25_v_%.2f.dat", v);
            write_to_file(filename, x, u, NX);
        }
        else if(n*dt > 0.5 - dt/2 && n*dt <= 0.5 + dt/2){
            sprintf(filename, "lw_0.5_v_%.2f.dat", v);
            write_to_file(filename, x, u, NX);
        }
        else if(n*dt > 0.75 - dt/2 && n*dt <= 0.75 + dt/2){
            sprintf(filename, "lw_0.75_v_%.2f.dat", v);
            write_to_file(filename, x, u, NX);
        }
        else if(n == NT-1){
            sprintf(filename, "lw_1.0_v_%.2f.dat", v);
            write_to_file(filename, x, u, NX);
        }
    }

    double error_lw = 0.0;

    for (int i = 0; i < NX; i++){
        double u_exact = exact_solution(x[i], T, c);
        error_lw += pow(u[i] - u_exact, 2);
    }

    error_lw = sqrt(error_lw * dx);

    printf("L2 Error (Lax-Wendroff) at T = %lf\n", error_lw);

    return 0;
}
