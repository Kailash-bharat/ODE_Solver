#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.141592653589793

int NX = 100;
double v = 2.0;
double dx, dt, c = 1.0, T = 1.0;

double initial_condition(double x){
    return sin(2 * PI * x);
}

double exact_solution(double x, double t){
    return sin(2 * PI * (x - t));
}

void write_to_file(char *filename, double *x, double *u, int Nx){
    FILE *fp = fopen(filename, "w");

    for(int i=0;i<Nx;i++){
        fprintf(fp, "%lf %lf\n", x[i], u[i]);
    }

    fclose(fp);
}

void MacCormack(double *u, double *u_new){
    double u_dash[NX];

    for(int i=0;i<NX;i++){
        int ip = (i+1)%NX;
        u_dash[i] = u[i] - v * (u[ip] - u[i]);
    }

    for(int i=0;i<NX;i++){
        int im = (i-1+NX)%NX;
        u_new[i] = 0.5 * (u[i] + u_dash[i] 
                    - v * (u_dash[i] - u_dash[im]));
    }
}

void CrankNicolson(double *u, double *u_new){
    double a = -v/4.0, b = 1.0, c = v/4.0;
    double rhs[NX];

    for(int i=0;i<NX;i++){
        int ip = (i+1)%NX;
        int im = (i-1+NX)%NX;
        rhs[i] = u[i] - (v/4.0)*(u[ip] - u[im]);
    }

    for(int i=0;i<NX;i++) u_new[i] = u[i];

    for(int iter=0; iter<500; iter++){
        for(int i=0;i<NX;i++){
            int ip = (i+1)%NX;
            int im = (i-1+NX)%NX;

            u_new[i] = (rhs[i] - a*u_new[im] - c*u_new[ip]) / b;
        }
    }
}

double L2_error(double *u, double *x, double t){
    double e = 0.0;

    for(int i=0;i<NX;i++){
        double u_exact = exact_solution(x[i], T);
        e += pow(u[i] - u_exact, 2);
    }

    return sqrt(dx * e);
}

int main(){
    dx = 0.01;
    dt = v * dx / c;

    int NT = (int)(T / dt);

    double x[NX];
    double u_mac[NX], u_cn[NX], u_new[NX];

    for(int i=0;i<NX;i++){
        x[i] = i * dx;
    }

    for(int i=0;i<NX;i++){
        u_mac[i] = initial_condition(x[i]);
        u_cn[i]  = initial_condition(x[i]);
    }

    char filename[50];

    sprintf(filename, "mac_0_v_%.2f.dat", v);
    write_to_file(filename, x, u_mac, NX);

    sprintf(filename, "cn_0_v_%.2f.dat", v);
    write_to_file(filename, x, u_cn, NX);

    for(int n=1; n<=NT; n++){

        MacCormack(u_mac, u_new);
        for(int i=0;i<NX;i++) u_mac[i] = u_new[i];

        CrankNicolson(u_cn, u_new);
        for(int i=0;i<NX;i++) u_cn[i] = u_new[i];


        if(n*dt > 0.25 - dt/2 && n*dt <= 0.25 + dt/2){
            sprintf(filename, "mac_0.25_v_%.2f.dat", v);
            write_to_file(filename, x, u_mac, NX);
        }
        if(n*dt > 0.5 - dt/2 && n*dt <= 0.5 + dt/2){
            sprintf(filename, "mac_0.5_v_%.2f.dat", v);
            write_to_file(filename, x, u_mac, NX);
        }
        if(n*dt > 0.75 - dt/2 && n*dt <= 0.75 + dt/2){
            sprintf(filename, "mac_0.75_v_%.2f.dat", v);
            write_to_file(filename, x, u_mac, NX);
        }
        if(n == NT - 1){
            sprintf(filename, "mac_1.0_v_%.2f.dat", v);
            write_to_file(filename, x, u_mac, NX);
        }

        if(n*dt >= 0.25 - dt/2 && n*dt <= 0.25 + dt/2){
            sprintf(filename, "cn_0.25_v_%.2f.dat", v);
            write_to_file(filename, x, u_cn, NX);
        }
        if(n*dt >= 0.5 - dt/2 && n*dt <= 0.5 + dt/2){
            sprintf(filename, "cn_0.5_v_%.2f.dat", v);
            write_to_file(filename, x, u_cn, NX);
        }
        if(n*dt >= 0.75 - dt/2 && n*dt <= 0.75 + dt/2){
            sprintf(filename, "cn_0.75_v_%.2f.dat", v);
            write_to_file(filename, x, u_cn, NX);
        }
        if(n == NT - 1){
            sprintf(filename, "cn_1.0_v_%.2f.dat", v);
            write_to_file(filename, x, u_cn, NX);
        }
    }

    double error_mac = L2_error(u_mac, x, T);
    double error_cn  = L2_error(u_cn, x, T);

    printf("L2 Error (MacCormack) at T =     %lf\n", error_mac);
    printf("L2 Error (Crank-Nicolson) at T = %lf\n", error_cn);
    if(isnan(error_mac) || isinf(error_mac)){
        printf("MacCormack unstable for v = %.2f\n", v);
    }

    return 0;
}
