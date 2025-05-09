#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void decompose(double e[], double f[], double g[], int n);
void substitute(double e[], double f[], double g[], double r[], int n, double d2x[]);
void tridiagonalMatrix(double x[], double y[], int n, double e[], double f[], double g[], double r[]);
void interpolation(double x[], double y[], int n, double d2x[], double xu, double *yu, double dy, double d2y);

void newtonInterpolation(double x[], double y[], int n, double xi, double yint[], double ea[], int iteration) {
    double fdd[n][n];
    int i, j, order;
    for (i = 0; i < n; i++) {
        fdd[i][0] = y[i];
    }
    for (j = 1; j < n; j++) {
        for (i = 0; i < n - j; i++) {
            fdd[i][j] = (fdd[i + 1][j - 1] - fdd[i][j - 1]) / (x[i + j] - x[i]);
        }
    }
    double xterm = 1.0;
    yint[0] = fdd[0][0];
    double yint2;
    for (order = 1; order < n; order++) {
        xterm *= (xi - x[order - 1]);
        yint2 = yint[order - 1] + fdd[0][order] * xterm;
        ea[order] = fabs((yint2 - yint[order - 1]) / yint2) * 100.0;
        yint[order] = yint2;
    }
}

void lagrangeInterpolation(double x[], double y[], int n, double xx, double *yint) {
    int i, j;
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        double product = y[i];
        for (j = 0; j < n; j++) {
            if (i != j) {
                product *= (xx - x[j]) / (x[i] - x[j]);
            }
        }
        sum += product;
    }
    *yint = sum;
}

void spline(double x[], double y[], int n, double xu, double *yu, double *dy, double *d2y) {
    double e[100], f[100], g[100], r[100], d2x[100];
    tridiagonalMatrix(x, y, n, e, f, g, r);
    decompose(e, f, g, n - 1);
    substitute(e, f, g, r, n - 1, d2x);
    interpolation(x, y, n, d2x, xu, yu, *dy, *d2y);
}

int readData(char *filename, double x[], double y[], int *n) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file!\n");
        return 0;
    }

    char line[256];
    int count = 0;

    while (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] < '0' || line[0] > '9') {
        } else {
            char *token = strtok(line, ",");
            if (token != NULL) {
                x[count] = atof(token);
                token = strtok(NULL, ",");
                if (token != NULL) {
                    y[count] = atof(token);
                    count++;
                } else {
                    printf("Warning: Missing y value in line %d\n", count + 1);
                }
            } else {
                printf("Warning: Missing x value in line %d\n", count + 1);
            }
        }
    }

    fclose(file);
    return count;
}

void main() {
    int i, n;
    double x[100], y[100], yint[100], ea[100], xi, dy, d2y;
    FILE *output_file = fopen("TugasPemrograman4_nomor2.csv", "w");
    if (output_file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    // Membaca data dari file CSV
    n = readData("data_interpolasi.csv", x, y, &n);
    if (n <= 0) {
        printf("Failed to read interpolation data or file is empty!\n");
        return;
    }

    xi = 2.0;
    double y_true = 0.69314718056; // True value for comparison
    double error_true;

    newtonInterpolation(x, y, n, xi, yint, ea, n);
    fprintf(output_file, "Newton Interpolation\n");
    fprintf(output_file, "Iteration,xi,yint,ea,et\n");
    for (i = 0; i < n; i++) {
        error_true = fabs((yint[i] - y_true) / y_true) * 100.0; // Calculate true error
        fprintf(output_file, "%d,%.10f,%.10f,%.10f,%.10f\n", i + 1, xi, yint[i], ea[i], error_true);
    }
    fprintf(output_file, "\n");

    yint[0] = 0.0; // Reset yint for Lagrange interpolation
    lagrangeInterpolation(x, y, n, xi, &yint[0]);
    error_true = fabs((yint[0] - y_true) / y_true) * 100.0; // Calculate true error
    fprintf(output_file, "Lagrange Interpolation\n");
    fprintf(output_file, "xi,yint,et\n");
    fprintf(output_file, "%.10f,%.10f,%.10f\n", xi, yint[0], error_true);
    fprintf(output_file, "\n");

    yint[0] = 0.0; // Reset yint for Spline interpolation
    dy = 0.0; d2y = 0.0; // Reset dy and d2y for Spline interpolation
    spline(x, y, n, xi, &yint[0], &ea[0], &ea[1]);
    error_true = fabs((yint[0] - y_true) / y_true) * 100.0; // Calculate true error
    fprintf(output_file, "Spline Interpolation\n");
    fprintf(output_file, "xi,yint,dy,d2y,et\n");
    fprintf(output_file, "%.10f,%.10f,%.10f,%.10f,%.10f\n", xi, yint[0], ea[0], ea[1], error_true);

    fclose(output_file);
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

void interpolation(double x[], double y[], int n, double d2x[], double xu, double *yu, double dy, double d2y) {
    int flag = 0;
    int i, j;
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
            dy = t1 + t2 + t3 + t4;

            t1 = 6 * c1 * (x[i] - xu);
            t2 = 6 * c2 * (xu - x[i - 1]);
            d2y = t1 + t2;

            flag = 1;
        } else {
            i = i + 1;
        }
        if (i == n + 1 || flag == 1) break;
    } while (flag != 1);

    if (flag == 0) {
        printf("Error: xu is out of range.\n");
    }
}