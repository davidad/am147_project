#include <stdlib.h>
#include <stdio.h>

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

int main(int argc, char** argv) {
	int n=3; // Number of masses in the system
	int i; //General purpose loop index, generally ranges from 0 to n-1
	double m[n]; // Masses (constants)
	double k[n]; // Spring constants
	double b[n]; // Damping constants
	double c[n]; // Coupling spring constants
	
	// Given a mass i, its current position is stored in xv[2*i],
	// while its current velocity is stored in xv[2*i+1].
	// Similarly, the time derivative of position is in dxdt[2*i]
	// while the time derivative of velocity is in dxdt[2*i+1]

	//This function, f, represnts all the equations of the dynamical system.
	void f(double t, double* xv, double* dxdt) {
		for(i=0;i<n;i++) {
			double x_i = xv[2*i];
			double v_i = xv[2*i+1];
			
			double x_im1;
			//xim1 is x_{i-1}. If i is zero, it's a boundary condition.
			if(i==0) {
				//For now the boundary condition is going to be 1.
				x_im1=1.0;
			} else {
			  x_im1 = xv[2*(i-1)];
			}

			double x_i_dot = v_i;
			double v_i_dot = - (1 / m[i]) * ( k[i]*x_i + c[i]*(x_i - x_im1) + b[i]*v_i);

			dxdt[2*i] = x_i_dot;
			dxdt[2*i+1] = v_i_dot;
		}
	}

	//Initialize constants
	for(i=0;i<n;i++) {
		m[i] = 1.0;
		k[i] = 1.0;
		b[i] = 0.0;
		c[i] = 1.0;
	}

	//Initialize initial values
	double x_init[2*n];
	for(i=0;i<n;i++) {
		x_init[2*i]   = 0.0; //Initial position of mass i
		x_init[2*i+1] = 0.0; //Initial velocity of mass i
	}
	double t_init = 0.0;

	int steps = 200; // How many steps (like frames of animation) to simulate
	double t_end = 5.0; // How many seconds those steps should cover
	double *t_out, **x_out; //This is where the results for each step get stored.
	alloc_out(2*n,steps,&t_out,&x_out); // Allocate memory for the above

	rkdumb(x_init,t_out,x_out,2*n,t_init,t_end,steps,f); //Actually run the simulation

	//Output results
	FILE* out;
	if(argc==2) {
		out=fopen(argv[1],"w");
	} else {
		out=stdout;
	}

	printf("[\n");
	int j; //frame number
	for(j=0;j<steps;j++) {
		fprintf(out,"  {\n");
		fprintf(out,"    \"t\": %lf,\n",t_init+j*((t_end-t_init)/steps));
		fprintf(out,"    \"x\": [ ");
		for(i=0;i<n;i++) {
			fprintf(out,(i!=n-1)?"%lf, ":"%lf ],\n",x_out[j][2*i]);
		}
		fprintf(out,"    \"v\": [ ");
		for(i=0;i<n;i++) {
			fprintf(out,(i!=n-1)?"%lf, ":"%lf ]\n",x_out[j][2*i+1]);
		}
		fprintf(out,(j!=steps-1)?"  },\n":"  }\n");
	}
	printf("]\n");
}
