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
	int steps = 240;
	double t0 = 0.0;
	double t1 = 6.0;

	double r = 0.2;
	double xinit[1] = {2.0};

	if(argc>1) {
		sscanf(argv[1],"%lf",&r);
	}
	if(argc>2) {
		sscanf(argv[2],"%lf",xinit);
	}

	double *tout, **xout;
	alloc_out(1,steps,&tout,&xout);

	void f(double t, double* xv, double* dxdt) {
		double x = xv[0];
		dxdt[0] = r*x - x/(1+x*x);
	}
	
	rkdumb(xinit,tout,xout,1,t0,t1,steps,f);

	int i;
	for(i=0;i<=steps;i++) {
		printf("%lf\t%lf\n",tout[i],xout[i][0]);
	}
}
