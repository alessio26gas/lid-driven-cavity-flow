#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define NX 128
#define NY 128


int main (void) {
    double X, Y, dx, dy;
    X = 1.0;
    Y = 1.0;
    dx = X/(NX-1);
    dy = Y/(NY-1);

    double Re, rho, nu;
    Re = 100.0;
    rho = 1.225;
    nu = 1.55E-5;

    double delta = 1.5;
    double CFL = 0.02;

    double u0 = 0.0;
    double v0 = 0.0;
    double p0 = 101325.0;

    int maxiter = 300000;
    double tol = 1.0E-5;

    clock_t start, end;
    double cpu_time;
    int i, j;

    int iter = 0;
    double error = 1.0;
    double maxerror = 0.0;

    double u_lid = Re * nu / X;
    double dt = CFL * dx / u_lid;

    double u[NX][NY+1], un[NX][NY+1], uc[NX][NY];
    double v[NX+1][NY], vn[NX+1][NY], vc[NX][NY];
    double p[NX+1][NY+1], pn[NX+1][NY+1], pc[NX][NY];
    double m[NX+1][NY+1];

    // Initial conditions
    for (i=0; i<=NX-1; i++) {
        for (j=0; j<=NY; j++) {
            u[i][j] = u0;
        }
    }

    for (i=0; i<=NX; i++) {
        for (j=0; j<=NY-1; j++) {
            v[i][j] = v0;
        }
    }

    for (i=0; i<=NX; i++) {
        for (j=0; j<=NY; j++) {
            p[i][j] = p0;
        }
    }

    start = clock();
    while (iter < maxiter && (iter > 0 || error/maxerror > tol)) {
        // Boundary conditions
        for (i=0; i<=(NX-1); i++) {
            un[i][0] = -un[i][1];                   // Bottom boundary u
            un[i][NY] = 2*u_lid - un[i][NY-1];      // Top boundary u
        }
        for (j=1; j<=(NY-1); j++) {
            un[0][j] = 0.0;                         // Left boundary u
            un[NX-1][j] = 0.0;                      // Right boundary u
        }
        for (i=0; i<=(NX); i++) {
            vn[i][0] = 0.0;                         // Bottom boundary v
            vn[i][NY-1] = 0.0;                      // Top boundary v
        }
        for (j=1; j<=(NY-2); j++) {
            vn[0][j] = -vn[1][j];                   // Left boundary v
            vn[NX][j] = -vn[NX-1][j];               // Right boundary v
        }        
        for (i=1; i<=(NX-1); i++) {
            pn[i][0] = pn[i][1];                    // Bottom boundary p
            pn[i][NY] = pn[i][NY-1];                // Top boundary p
        }
        for (j=0; j<=(NY); j++) {
            pn[0][j] = pn[1][j];                    // Left boundary p
            pn[NX][j] = pn[NX-1][j];                // Right boundary p
        }

        // X-momentum equation
        for (i=1; i<=(NX-2); i++) {
            for (j=1; j<=(NY-1); j++) {
                un[i][j] = u[i][j] - dt*((u[i+1][j]*u[i+1][j]-u[i-1][j]*u[i-1][j])/2.0/dx 
                           + 0.25*((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]))/dy)
                           - dt/rho/dx*(p[i+1][j]-p[i][j]) 
                           + dt*nu*((u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx/dx +(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy/dy);
            }
        }

        // Y-momentum equation
        for (i=1; i<=(NX-1); i++) {
            for (j=1; j<=(NY-2); j++) {
                vn[i][j] = v[i][j] - dt*(0.25*((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i-1][j]+u[i-1][j+1])*(v[i][j]+v[i-1][j]))/dx 
                           +(v[i][j+1]*v[i][j+1]-v[i][j-1]*v[i][j-1])/2.0/dy) 
                           - dt/rho/dy*(p[i][j+1]-p[i][j]) 
                           + dt*nu*((v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy/dy);
            }
        }

        // Continuity equation
        for (i=1; i<=(NX-1); i++) {
            for (j=1; j<=(NY-1); j++)
            {
                pn[i][j] = p[i][j]-rho*u_lid*u_lid*dt*delta*((un[i][j]-un[i-1][j])/dx + (vn[i][j]-vn[i][j-1]) /dy);
            }
        }
    
        // Display normalized error
        error = 0.0;
        for (i=1; i<=(NX-1); i++) {
            for (j=1; j<=(NY-1); j++) {
                m[i][j] = ((un[i][j]-un[i-1][j])/dx + (vn[i][j]-vn[i][j-1])/dy);
                error = error + fabs(m[i][j]);
            }
        }

        if (error > maxerror) {
            maxerror = error;
        }

        if (iter % 1000 == 0) {
            printf("Step #%d \t%5.6lf\n", iter, error/maxerror);
        }

        for (i=0; i<=(NX-1); i++) {
            for (j=0; j<=(NY); j++) {
                u[i][j] = un[i][j];
            }
        }

        for (i=0; i<=(NX); i++) {
            for (j=0; j<=(NY-1); j++) {
                v[i][j] = vn[i][j];
            }
        }

        for (i=0; i<=(NX); i++) {
            for (j=0; j<=(NY); j++) {
                p[i][j] = pn[i][j];
            }
        }

        iter++;
    }
    end=clock();
    cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Simulation completed in %lf seconds at step #%d\n", cpu_time, iter);

    for (i=0; i<=(NX-1); i++) {
        for (j=0; j<=(NY-1); j++) {    
            uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
            vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
            pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
        }
    }
    
    // Write output data
    FILE *file;
    file = fopen("output.vtk", "w");
    if (file == NULL) {
        printf("\nERROR when opening files\n");
    } else {
        fprintf(file, "# vtk DataFile Version 4.2\n");
        fprintf(file, "Velocity and Pressure Data\n");
        fprintf(file, "ASCII\n");
        fprintf(file, "DATASET STRUCTURED_GRID\n");
        fprintf(file, "DIMENSIONS %d %d 1\n", NX, NY);
        fprintf(file, "POINTS %d float\n", NX * NY);

        // Writing point coordinates
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                double xpos = i * dx;
                double ypos = j * dy;
                fprintf(file, "%5.8lf %5.8lf 0.0\n", xpos, ypos);
            }
        }
        fprintf(file, "\nPOINT_DATA %d\n", NX * NY);

        // Writing velocity vector
        fprintf(file, "VECTORS V float\n");
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                fprintf(file, "%5.8lf %5.8lf 0.0\n", uc[i][j], vc[i][j]);
            }
        }

        // Writing normalized velocity vector
        fprintf(file, "VECTORS V/V0 float\n");
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                fprintf(file, "%5.8lf %5.8lf 0.0\n", uc[i][j]/u_lid, vc[i][j]/u_lid);
            }
        }

        // Writing pressure
        fprintf(file, "SCALARS p float 1\n");
        fprintf(file, "LOOKUP_TABLE default\n");
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                fprintf(file, "%5.8lf\n", pc[i][j]);
            }
        }

        // Writing normalized pressure
        fprintf(file, "SCALARS p/q float 1\n");
        fprintf(file, "LOOKUP_TABLE default\n");
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                fprintf(file, "%5.8lf\n", pc[i][j]/rho/u_lid/u_lid);
            }
        }

        fclose(file);
        printf("Data successfully written to output.vtk\n");
    }
}