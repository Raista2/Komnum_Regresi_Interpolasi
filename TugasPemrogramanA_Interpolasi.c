#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void decompose(double e[], double f[], double g[], int n);
void substitute(double e[], double f[], double g[], double r[], int n, double d2x[]);
void tridiagonalMatrix(double x[], double y[], int n, double e[], double f[], double g[], double r[]);
void interpolation(double x[], double y[], int n, double d2x[], double xu, double *yu, double *dy, double *d2y);
void spline(double x[], double y[], int n, double xu, double *yu, double *dy, double *d2y);
int readData(char *filename, double x[], double y[], int *n);

void main() {
    int n;
    double x[100], y[100], yint[100], dy, d2y;
    double years_to_predict[] = {2005.0, 2006.0, 2015.0, 2016.0};
    int num_years = 4;

    FILE *output_file = fopen("TugasPemrogramanA_nomor1.csv", "w");
    if (output_file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    // Prediksi populasi
    n = readData("DataNo1_Populasi.csv", x, y, &n);
    if (n <= 0) {
        printf("Failed to read population data!\n");
        return;
    }
    
    fprintf(output_file, "Cubic Spline Interpolation for Population Missing Years\n\n");
    fprintf(output_file, "Year,Predicted Population,Growth Rate,Acceleration\n");
    
    printf("Population predictions using Cubic Spline Interpolation:\n");
    printf("------------------------------------------------------\n");
    
    for (int year_idx = 0; year_idx < num_years; year_idx++) {
        double xi = years_to_predict[year_idx];
        dy = 0.0; d2y = 0.0;
        
        spline(x, y, n, xi, &yint[0], &dy, &d2y);
        
        fprintf(output_file, "%.0f,%.0f,%.2f,%.2f\n", 
                xi, yint[0], dy, d2y);
        
        printf("Year %.0f: %.0f people (growth rate: %.0f people/year)\n", xi, yint[0], dy);
    }
    
    // Prediksi persentase internet
    fprintf(output_file, "\n\nCubic Spline Interpolation for Internet Usage Percentage\n\n");
    fprintf(output_file, "Year,Predicted Percentage,Rate of Change,Acceleration\n");
    
    n = readData("DataNo1_Persentasi.csv", x, y, &n);
    if (n <= 0) {
        printf("Failed to read internet percentage data!\n");
        return;
    }

    printf("\nInternet usage percentage predictions using Cubic Spline Interpolation:\n");
    printf("--------------------------------------------------------------\n");
    
    for (int year_idx = 0; year_idx < num_years; year_idx++) {
        double xi = years_to_predict[year_idx];
        dy = 0.0; d2y = 0.0;
        
        spline(x, y, n, xi, &yint[0], &dy, &d2y);
        
        if (yint[0] < 0) {yint[0] = 0;}
        if (yint[0] > 100) {yint[0] = 100;}

        fprintf(output_file, "%.0f,%.2f,%.2f,%.2f\n", xi, yint[0], dy, d2y);
        
        printf("Year %.0f: %.2f%% (change rate: %.2f%%/year)\n", xi, yint[0], dy);
    }
    
    fclose(output_file);
    printf("\nInterpolation complete. Results saved to TugasPemrograman4_nomor2.csv\n");
}

void decompose(double e[], double f[], double g[], int n) {
    int i;
    for (i = 2; i <= n; i++) {
        e[i] = e[i] / f[i - 1];
        f[i] = f[i] - e[i] * g[i - 1];
    }
}

void substitute(double e[], double f[], double g[], double r[], int n, double d2x[]) {
    int i;
    r[1] = r[1] / f[1];
    for (i = 2; i <= n; i++) {
        r[i] = (r[i] - e[i] * r[i - 1]) / f[i];
    }
    d2x[n] = r[n];
    for (i = n - 1; i >= 1; i--) {
        d2x[i] = r[i] - g[i] * d2x[i + 1];
    }
}

void tridiagonalMatrix(double x[], double y[], int n, double e[], double f[], double g[], double r[]) {
    int i;
    f[1] = 2 * (x[2] - x[0]);
    g[1] = x[2] - x[1];
    r[1] = 6 / (x[2] - x[1]) * (y[2] - y[1]);
    r[1] += 6 / (x[1] - x[0]) * (y[0] - y[1]);

    for (i = 2; i < n - 1; i++) {
        e[i] = x[i] - x[i - 1];
        f[i] = 2 * (x[i + 1] - x[i - 1]);
        g[i] = x[i + 1] - x[i];
        r[i] = 6 / (x[i + 1] - x[i]) * (y[i + 1] - y[i]);
        r[i] += 6 / (x[i] - x[i - 1]) * (y[i-1] - y[i]);
    }

    e[n - 1] = x[n - 1] - x[n - 2];
    f[n - 1] = 2 * (x[n] - x[n - 2]);
    r[n - 1] = 6 / (x[n] - x[n - 1]) * (y[n] - y[n - 1]);
    r[n - 1] += 6 / (x[n - 1] - x[n - 2]) * (y[n - 2] - y[n - 1]);
}

void interpolation(double x[], double y[], int n, double d2x[], double xu, double *yu, double *dy, double *d2y) {
    int flag = 0;
    int i;
    i = 1;

    do {
        if (xu >= x[i - 1] && xu <= x[i]) {
            double c1 = d2x[i - 1] / 6.0 * (x[i] - x[i - 1]);
            double c2 = d2x[i] / 6.0 * (x[i] - x[i - 1]);
            double c3 = y[i - 1] / (x[i] - x[i - 1]) - d2x[i - 1] * (x[i] - x[i - 1]) / 6.0;
            double c4 = y[i] / (x[i] - x[i - 1]) - d2x[i] * (x[i] - x[i - 1]) / 6.0;

            double t1 = c1 * pow(x[i] - xu, 3);
            double t2 = c2 * pow(xu - x[i - 1], 3);
            double t3 = c3 * (x[i] - xu);
            double t4 = c4 * (xu - x[i - 1]);
            *yu = t1 + t2 + t3 + t4;

            t1 = -3 * c1 * pow(x[i] - xu, 2);
            t2 = 3 * c2 * pow(xu - x[i - 1], 2);
            t3 = -c3;
            t4 = c4;
            *dy = t1 + t2 + t3 + t4;

            t1 = 6 * c1 * (x[i] - xu);
            t2 = 6 * c2 * (xu - x[i - 1]);
            *d2y = t1 + t2;

            flag = 1;
        } else {
            i = i + 1;
        }
        if (i == n + 1 || flag == 1) break;
    } while (flag != 1);

    if (flag == 0) {
        printf("Error: xu is out of range.\n");
        *yu = 0;
        *dy = 0;
        *d2y = 0;
    }
}

void spline(double x[], double y[], int n, double xu, double *yu, double *dy, double *d2y) {
    double e[100], f[100], g[100], r[100], d2x[100];
    tridiagonalMatrix(x, y, n, e, f, g, r);
    decompose(e, f, g, n - 1);
    substitute(e, f, g, r, n - 1, d2x);
    interpolation(x, y, n, d2x, xu, yu, dy, d2y);
}

int readData(char *filename, double x[], double y[], int *n) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file %s!\n", filename);
        return 0;
    }

    char line[256];
    int count = 0;

    fgets(line, sizeof(line), file);
    
    while (fgets(line, sizeof(line), file) != NULL && count < 100) {
        char *token = strtok(line, ",");
        if (token != NULL) {
            x[count] = atof(token);
            token = strtok(NULL, ",");
            if (token != NULL) {
                y[count] = atof(token);
                count++;
            }
        }
    }

    fclose(file);
    *n = count;
    return count;
}