#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX 100
#define MAXORDER 10

void gaussElimination(int order, double a[order + 1][order + 2], double coeffs[]);

void linRegress (double x[], double y[], int n, double *a1, double *a0, double *syx, double *r2, FILE *output_file) {
    fprintf(output_file, "Iteration, x, y, y_pred, MSE\n");
    double sumx, sumy, sumxy, sumx2, st, sr;
    sumx = 0;
    sumy = 0;
    sumxy = 0;
    sumx2 = 0;
    st = 0;
    sr = 0;
    int i;
    for (i = 0; i < n; i++) {
        sumx += x[i];
        sumy += y[i];
        sumxy += x[i] * y[i];
        sumx2 += x[i] * x[i];
    }
    double xmean = sumx / n;
    double ymean = sumy / n;
    *a1 = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);
    *a0 = ymean - *a1 * xmean;

    for (i = 0; i < n; i++) {
        double y_pred = *a1 * x[i] + *a0;
        double mse = (y[i] - y_pred) * (y[i] - y_pred);
        fprintf(output_file, "%d,%.2f,%.2f,%.2f,%.2f\n", i + 1, x[i], y[i], y_pred, mse);
        st += (y[i] - ymean) * (y[i] - ymean);
        sr += (y[i] - y_pred) * (y[i] - y_pred);
    }

    *syx = sqrt(sr / (n - 2));
    *r2 = (st - sr) / st;
}

void polyRegress(double x[], double y[], int order, int n, double a[order + 1][order + 2], double *syx, double *r2, FILE *output_file) {
    fprintf(output_file, "Iteration, x, y, y_pred, MSE\n");

    double coeffs[order + 1], st, sr, sum_y;
    int i, j, k;
    st = 0;
    sr = 0;
    // Mengisi matriks normal equation
    for (i = 0; i <= order; i++) {
        for (j = 0; j <= order; j++) {
            double sum = 0;
            for (k = 0; k < n; k++) {
                sum += pow(x[k], i + j);
            }
            a[i][j] = sum;
        }

        sum_y = 0;
        for (k = 0; k < n; k++) {
            sum_y += y[k] * pow(x[k], i);
        }
        a[i][order + 1] = sum_y;
    }
    // Menyelesaikan persamaan normal menggunakan Eliminasi Gauss
    gaussElimination(order, a, coeffs);

    // Menghitung y_pred dan MSE
    double ymean = sum_y / n;
    for (i = 0; i < n; i++) {
        double y_pred = 0;
        for (j = 0; j <= order; j++) {
            y_pred += coeffs[j] * pow(x[i], j);
        }
        double mse = (y[i] - y_pred) * (y[i] - y_pred);
        fprintf(output_file, "%d,%.2f,%.2f,%.2f,%.2f\n", i + 1, x[i], y[i], y_pred, mse);
        st += (y[i] - ymean) * (y[i] - ymean);
        sr += (y[i] - y_pred) * (y[i] - y_pred);
    }
    
    *syx = sqrt(sr / (n - order - 1));
    *r2 = (st - sr) / st;
}

void multiRegress (double x[][MAXORDER], double y[], int order, double a[][order+2], int n, double *syx, double *r2, FILE *output_file) {
	int i, j, k;
    fprintf(output_file, "Iteration,");
    for (i = 1; i <= order; i++) {
        fprintf(output_file, "x[%d],", i);
    }
    fprintf(output_file, "y,y_pred,MSE\n");
    double coeffs[order + 1], st, sr, sum_y;
    st = 0;
    sr = 0;

    for (i = 0; i <= order; i++) {
        for (j = 0; j <= i; j++) {
            double sum = 0;
            for (k = 0; k < n; k++) {
                sum += x[k][i] * x[k][j];
            }
            a[i][j] = sum;
            a[j][i] = sum;
        }
        sum_y = 0;
        for (k = 0; k < n; k++) {
            sum_y += y[k] * x[k][i];
        }
        a[i][order + 1] = sum_y;
    }

    gaussElimination(order, a, coeffs);
    // Menghitung y_pred dan MSE
    double ymean = sum_y / n;
    for (i = 0; i < n; i++) {
        double y_pred = 0;
        fprintf(output_file, "%d,", i+1);
        for (j = 1; j <= order; j++) {
            fprintf(output_file, "%.2f,", x[i][j]);
        }
        for (j = 0; j <= order; j++) {
            y_pred += coeffs[j] * x[i][j];
        }
        double mse = (y[i] - y_pred) * (y[i] - y_pred);
        fprintf(output_file, "%.2f,%.2f,%.2f\n", y[i], y_pred, mse);
        st += (y[i] - ymean) * (y[i] - ymean);
        sr += (y[i] - y_pred) * (y[i] - y_pred);
    }
    *syx = sqrt(sr / (n - order - 1));
    *r2 = (st - sr) / st;
}

int readDataRegresi(char *filename, double x[], double y[], int *n) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file!\n");
        return 0;
    }

    char line[256];
    int count = 0;

    if (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] < '0' || line[0] > '9') {
        } else {
            char* token = strtok(line, ",");
            if (token != NULL) {
                x[count] = atof(token);
                token = strtok(NULL, ",");
                if (token != NULL) {
                    y[count] = atof(token);
                    count++;
                }
            }
        }
    }

    while (fgets(line, sizeof(line), file) != NULL && count < MAX) {
        char* token = strtok(line, ",");
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
    return count;
}

int readDataMultiple(char *filename, double x[][MAXORDER], double y[], int max_size, int order) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file! %s\n", filename); 
        return 0;
    }

    char line[256];
    int count = 0;

    if (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] < '0' || line[0] > '9') {
        } else {
            x[count][0] = 1.0;
            
            char* token = strtok(line, ",");
            int col = 1;
            
            while (token != NULL && col <= order) {
                x[count][col] = atof(token);
                token = strtok(NULL, ",");
                col++;
            }
            
            if (token != NULL) {
                y[count] = atof(token);
                count++;
            } else {
                printf("Warning: Missing y value in line %d\n", count+1);
            }
        }
    }

    while (fgets(line, sizeof(line), file) != NULL && count < max_size) {
        x[count][0] = 1.0; 
        
        char* token = strtok(line, ",");
        int col = 1;
        
        while (token != NULL && col <= order) {
            x[count][col] = atof(token);
            token = strtok(NULL, ",");
            col++;
        }
        
        if (token != NULL) {
            y[count] = atof(token);
            count++;
        } else {
            printf("Warning: Missing y value in line %d\n", count+1);
        }
    }

    fclose(file);
    return count;
}

int main() {
    int i, j;
    FILE *output_file;
    output_file = fopen("TugasPemrograman4_nomor1.csv", "w");
    if (output_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }  

    // Membaca data dari file CSV
    double x_linear[MAX], y_linear[MAX];
    double x_polynomial[MAX], y_polynomial[MAX];
    double x_multi[MAX][MAXORDER], y_multi[MAX];
    
    // Declare variables to hold the number of data points
    int n_linear_count;
    int n_poly_count;
    
    // Membaca data untuk regresi linear
    n_linear_count = readDataRegresi("data_linear.csv", x_linear, y_linear, &n_linear_count);
    if (n_linear_count <= 0) {
        printf("Failed to read linear regression data!\n");
        return 1;
    }
    
    // Membaca data untuk regresi polinomial
    n_poly_count = readDataRegresi("data_polynomial.csv", x_polynomial, y_polynomial, &n_poly_count);
    if (n_poly_count <= 0) {
        printf("Failed to read polynomial regression data!\n");
        return 1;
    }
    
    // Membaca data untuk regresi berganda
    int order = 2; // Order sesuai dengan kebutuhan
    int n_multi_count;
    n_multi_count = readDataMultiple("data_multiple.csv", x_multi, y_multi, MAX, order);
    if (n_multi_count <= 0) {
        printf("Failed to read multiple regression data!\n");
        return 1;
    }
    
    // Process the read data using the correct count variables
    double a1, a0, syx, r2;
    double a[MAXORDER][MAXORDER + 2];

    // Linear Regression
    fprintf(output_file, "Linear Regression:\n");
    linRegress(x_linear, y_linear, n_linear_count, &a1, &a0, &syx, &r2, output_file);
    if (a0 < 0) {
    	a0 *= -1;
    	fprintf(output_file, "\nf(x) = %.6fx - %.6f\n syx, %.6f\n r2, %.6f\n\n\n", a1, a0, syx, r2);
	} else {
		fprintf(output_file, "\nf(x) = %.6fx + %.6f\n syx, %.6f\n r2, %.6f\n\n\n", a1, a0, syx, r2);
	}
    
    // Polynomial Regression
    fprintf(output_file, "Polynomial Regression:\n");
    polyRegress(x_polynomial, y_polynomial, order, n_poly_count, a, &syx, &r2, output_file);
    double coeffs[order + 1];  // Array untuk menyimpan koefisien regresi
    gaussElimination(order, a, coeffs);
    fprintf(output_file, "\nOrder, %d\nf(x) = ", order);
    for (i = order; i > 0; i--) {
    	if (i == order) {
    		fprintf(output_file, "%.6fx^%d ", coeffs[i], i);
		} else {
			if (coeffs[i] < 0) {
	    		coeffs[i] *= -1;
	    		fprintf(output_file, "- %.6fx^%d ", coeffs[i], i);
	    	} else {
	    		fprintf(output_file, "+ %.6fx^%d ", coeffs[i], i);
			}
		}
    	
	}
	if (coeffs[0] < 0) {
    	coeffs[0] *= -1;
    	fprintf(output_file, "- %.6f\n syx, %.6f\n r2, %.6f\n\n\n", coeffs[0], syx, r2);
    } else {
    	fprintf(output_file, "+ %.6f\n syx, %.6f\n r2, %.6f\n\n\n", coeffs[0], syx, r2);	
	}

    // Multiple Regression
    fprintf(output_file, "Multiple Regression:\n");
    multiRegress(x_multi, y_multi, order, a, n_multi_count, &syx, &r2, output_file);
    gaussElimination(order, a, coeffs);
    fprintf(output_file, "\nOrder, %d\nf(x) = ", order);
    for (i = order; i > 0; i--) {
    	if (i == order) {
    		fprintf(output_file, "%.6fx%d ", coeffs[i], i);
		} else {
			if (coeffs[i] < 0) {
	    		coeffs[i] *= -1;
	    		fprintf(output_file, "- %.6fx%d ", coeffs[i], i);
	    	} else {
	    		fprintf(output_file, "+ %.6fx%d ", coeffs[i], i);
			}
		}
    	
	}
	if (coeffs[0] < 0) {
    	coeffs[0] *= -1;
    	fprintf(output_file, "- %.6f\n syx, %.6f\n r2, %.6f\n\n\n", coeffs[0], syx, r2);
    } else {
    	fprintf(output_file, "+ %.6f\n syx, %.6f\n r2, %.6f\n\n\n", coeffs[0], syx, r2);	
	}
    fclose(output_file);
    return 0;
}

void gaussElimination(int order, double a[order + 1][order + 2], double coeffs[]) {
    int i, j, k;
    for (i = 0; i <= order; i++) {
        // Pivoting (opsional, untuk stabilitas numerik)
        for (j = i + 1; j <= order; j++) {
            if (fabs(a[j][i]) > fabs(a[i][i])) {
                for (k = 0; k <= order + 1; k++) {
                    double temp = a[i][k];
                    a[i][k] = a[j][k];
                    a[j][k] = temp;
                }
            }
        }

        // Eliminasi Gauss
        for (j = i + 1; j <= order; j++) {
            double factor = a[j][i] / a[i][i];
            for (k = i; k <= order + 1; k++) {
                a[j][k] -= factor * a[i][k];
            }
        }
    }

    // Back-substitution
    for (i = order; i >= 0; i--) {
        coeffs[i] = a[i][order + 1];
        for (j = i + 1; j <= order; j++) {
            coeffs[i] -= a[i][j] * coeffs[j];
        }
        coeffs[i] /= a[i][i];
    }
}
