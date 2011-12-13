#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
/* Dependency on GSL: to install on OSX, sudo port install gsl-devel */
#include "gsl/gsl_linalg.h"

/* GENERIC FOURTH-ORDER RUNGE-KUTTA IMPLEMENTATION */

//derivs takes arguments (t,*x,*dxdt)
void rk4(double* x, double* dxdt, int n, double t, double h, double* xout, void (*derivs)(double, double*, double*)) {
	int i;
	double* dxm = malloc(n*sizeof(double));
	double* dxt = malloc(n*sizeof(double));
	double* xt  = malloc(n*sizeof(double));
	double hh = h*0.5;
	double h6 = h/6.0;
	double th = t+hh;
	for(i=0;i<n;i++) xt[i]=x[i]+hh*dxdt[i];
	(*derivs)(th,xt,dxt);
	for(i=0;i<n;i++) xt[i]=x[i]+hh*dxt[i];
	(*derivs)(th,xt,dxm);
	for(i=0;i<n;i++) {
		xt[i]=x[i]+h*dxm[i];
		dxm[i] += dxt[i];
	}
	(*derivs)(t+h,xt,dxt);
	for(i=0;i<n;i++) {
		xout[i] = x[i]+h6*(dxdt[i]+dxt[i]+2.0*dxm[i]);
	}
	free(dxm);
	free(dxt);
	free(xt);
}

//Remember, tout and xout must have steps+1 cells allocated!
void rkdumb(double* xinit, double* tout, double** xout, int n, double t0, double t1, int steps, void (*derivs)(double, double*, double*)) {
	int i,k;
	double t,h;
	double *dv = malloc(n*sizeof(double));
	for(i=0;i<n;i++) {
		xout[0][i]=xinit[i];
	}
	tout[0]=t0;
	t=t0;
	h=(t1-t0)/steps;
	for(k=0;k<steps;k++) {
		(*derivs)(t,xout[k],dv);
		rk4(xout[k],dv,n,t,h,xout[k+1],derivs);
		if((double)(t+h) == t) {
			fprintf(stderr, "Step size too small!\n");
		}
		t+=h;
		tout[k+1]=t;
	}
	free(dv);
}

void alloc_out(int n, int steps, double** tout, double*** xout) {
	int i;
	*tout = malloc((steps+1)*sizeof(double));
	*xout = malloc((steps+1)*sizeof(double*));
	for(i=0;i<=steps;i++) {
		(*xout)[i] = malloc(n*sizeof(double));
	}
}

/* END GENERIC RUNGE-KUTTA IMPLEMENTATION */

int main(int argc, char** argv) {

/* Algorithm and values from (Diependaal 1988) */

  double ell_max=35.0; // (mm) length of basilar membrane
	int n=512; // Number of discrete steps for ell
  double delta = ell_max / (double)n; // Size of discrete steps in ell
	int i; //General purpose loop index, generally ranges from 0 to n-1
  double ell[n+1];
  for(i=0;i<=n;i++) {
    ell[i] = delta*i;
  }

  /* Time-independent parameters */
  double m(double ell) { // mass of membrane unit
    return 0.5; // (mg/mm^2)
  }
  double alpha(double ell) { // alpha = 2 * density * membrane_thickness / mass * channel_area
    return 0.4; // (mm^-2)
  }
  
  /* Behavior of basilar membrane */
  double g(double ell, double x, double v) {
    double k(double ell, double x, double v) {
      return 20000*exp(-0.3*ell); // (mg mm^-2 ms^-2)
    }
    double b(double ell, double x, double v) {
      return 10*exp(-0.15*ell); // (mg mm^-2 ms^-1)
    }
    return k(ell,x,v)*x + b(ell,x,v)*v;
  }

  /* Forcing/driving/input function */
	double input(double t) {
		double scale=1.0;
    return scale*sin(t*8);
	}

  /* Membrane displacement and velocity */
  double xv[2*n];
	// Given a segment i, its current position is stored in xv[2*i],
	// while its current velocity is stored in xv[2*i+1].
  // The initial condition is rest:
  for(i=0;i<n;i++) {
    xv[2*i] = 0.0; // zero displacement
    xv[2*i+1] = 0.0; // zero velocity
  }

  /* Spatial system */
  // First compute the tridiagonal A matrix (not time-dependent)
  gsl_vector* A_diag = gsl_vector_alloc(n);
  gsl_vector* A_upper = gsl_vector_alloc(n-1);
  gsl_vector* A_lower = gsl_vector_alloc(n-1);
  gsl_vector_set(A_diag,0,0.5*alpha(ell[0])*delta+1/delta);
  for(i=1;i<n;i++) {
    gsl_vector_set(A_diag,i,delta*alpha(ell[i])+2/delta);
  }
  gsl_vector_set_all(A_upper,-1/delta);
  gsl_vector_set_all(A_lower,-1/delta);

  //Reusable allocations
  gsl_vector* k = gsl_vector_alloc(n);
  gsl_vector p_vec = {.data = NULL, .size = n, .stride = 1, .block = NULL, .owner = 0};

  void solve_spatial_system(double* p, double* xv, double forcing) {
    // Compute k vector
    gsl_vector_set(k,0,0.5*alpha(ell[0])*delta*g(ell[0],xv[0],xv[1])-forcing);
    for(i=1;i<n;i++) {
      gsl_vector_set(k,i,alpha(ell[i])*delta*g(ell[i],xv[2*i],xv[2*i+1]));
    }

    // Set up p_vec struct
    p_vec.data = p;

    // Solve system
    gsl_linalg_solve_tridiag(A_diag, A_upper, A_lower, k, &p_vec);
  }

  /* Temporal system */
  double p[n];
  //from t and xv (x and \dot{x}), we must calculate dx/dt (\dot{x} and \ddot{x}).
  void f(double t, double* xv, double* dxdt) {
    //First we must solve the spatial system in P (transmembrane pressure).
    solve_spatial_system(p,xv,input(t));
    for(i=0;i<n;i++) {
      double x_i = xv[2*i];
      double v_i = xv[2*i+1];

      double x_i_dot = v_i;
      double v_i_dot = (1 / m(ell[i])) * (p[i] - g(ell[i], x_i, v_i));

      dxdt[2*i] = x_i_dot;
      dxdt[2*i+1] = v_i_dot;
    }
  }
  

  /* Simulation parameters */
  double t_max = 41.0; // ms
  int iterations = (int)(500.0 * t_max);
  double h = t_max / (double)iterations;
  int output_filter = 16; // output the result of only every Nth iteration

  double t = 0.0;
  double dv[2*n];
  double xv_new[2*n];

  //JSON output setup
	FILE* out;
	if(argc==2) {
		out=fopen(argv[1],"w");
	} else {
		out=stdout;
	}
  double clamp(double z) {
    if(z > 1e9) return 1e9;
    else if(z< -1e9) return -1e9;
    else if isnan(z) return 1e9;
    else return z;
  }
  fprintf(out,"[\n");

  //It's go time!
  int j;
  for(j=0;j<iterations;j++) {
    f(t,xv,dv);
    rk4(xv,dv,2*n,t,h,xv_new,f);
    // Output==
    if(j%output_filter==0) {
      fprintf(out,"  {\n");
      fprintf(out,"    \"t\": %lf,\n",t);
      fprintf(out,"    \"x\": [ ");
      for(i=0;i<n;i++) {
        fprintf(out,(i!=n-1)?"%lf, ":"%lf ],\n",clamp(xv[2*i]));
      }
      fprintf(out,"    \"v\": [ ");
      for(i=0;i<n;i++) {
        fprintf(out,(i!=n-1)?"%lf, ":"%lf ],\n",clamp(xv[2*i+1]));
      }
      fprintf(out,"    \"p\": [ ");
      for(i=0;i<n;i++) {
        fprintf(out,(i!=n-1)?"%lf, ":"%lf ]\n",clamp(p[i]));
      }
      fprintf(out,(j+output_filter<iterations)?"  },\n":"  }\n");
    }
    // ========
    t+=h;
    memcpy(xv,xv_new,sizeof(double)*2*n);
  }

  fprintf(out,"]\n");

  return 0;
}
