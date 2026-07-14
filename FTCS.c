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

int main()
{
    double c = 1.0;
    double L = 1.0;
    double dx = 0.01;
    double v = 1.15;
    double dt = v * dx / c;
    double T = 1.0;
    int NX = 100;
    int NT = (int)(T/dt);

    double x[NX];
    double u_ftcs[NX];
    double u_new[NX];
    double max_ftcs[NT];

    for (int i = 0; i < NX; i++){
        x[i] = i * dx;
        u_ftcs[i] = initial_condition(x[i]);
    }

    char filename[50];
    sprintf(filename, "ftcs_0_v_%.2f.dat", v);
    write_to_file(filename, x, u_ftcs, NX);

    for (int n = 0; n < NT; n++)
    {
        for (int i = 1; i < NX - 1; i++){
            u_new[i] = u_ftcs[i] - 0.5 * v * (u_ftcs[i+1] - u_ftcs[i-1]);
        }

        u_new[0] = u_ftcs[0] - 0.5 * v * (u_ftcs[1] - u_ftcs[NX-1]);
        u_new[NX-1] = u_ftcs[NX-1] - 0.5 * v * (u_ftcs[0] - u_ftcs[NX-2]);

        double max_val = fabs(u_new[0]);

        for (int i = 1; i < NX; i++){
            if (fabs(u_new[i]) > max_val){
                max_val = fabs(u_new[i]);
            }
        }

        max_ftcs[n] = max_val;

        for (int i = 0; i < NX; i++)
            u_ftcs[i] = u_new[i];
        
        if(n*dt > 0.25 - dt/2 && n*dt <= 0.25 + dt/2){
            sprintf(filename, "ftcs_0.25_v_%.2f.dat", v);
            write_to_file(filename, x, u_ftcs, NX);
        }
        else if(n*dt > 0.5 - dt/2 && n*dt <= 0.5 + dt/2){
            sprintf(filename, "ftcs_0.5_v_%.2f.dat", v);
            write_to_file(filename, x, u_ftcs, NX);
        }
        else if(n*dt > 0.75 - dt/2 && n*dt <= 0.75 + dt/2){
            sprintf(filename, "ftcs_0.75_v_%.2f.dat", v);
            write_to_file(filename, x, u_ftcs, NX);
        }
        else if(n == NT-1){
            sprintf(filename, "ftcs_1.0_v_%.2f.dat", v);
            write_to_file(filename, x, u_ftcs, NX);
        }
    }

    sprintf(filename, "ftcs_absmax_v_%.2f.dat", v);
    FILE *fp = fopen(filename, "w");

    for (int n = 0; n < NT; n++){
        fprintf(fp, "%lf %lf\n", n * dt, max_ftcs[n]);
    }

    fclose(fp);

    return 0;
}
