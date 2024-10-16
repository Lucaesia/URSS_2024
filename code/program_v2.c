#include "grid/multigrid.h"
#include "run.h"
#include "diffusion.h"
#include "view.h"
#include "runge-kutta.h"
#include <math.h>

#define EPS 1e-14
#define ALPHA 1
#define BETA 1
#define Q 1
#define M 1
#define D_a 1
#define c_0 1
#define a_0 1

#define resolution 100

scalar c[], a[];
double s,s_t;


double t_end = 1.01;


// Runtime statistics output stream
FILE * fp_stats;


// function for runge-kutter
static void dc (scalar * ul, double t, scalar * kl)
{
  fprintf(stderr,"%g\n",s);
  scalar u = ul[0], k = kl[0];
  foreach() {
    k[] = u[] + dt*(u[1] + u[-1] -2*u[])/sq(Delta*s) + dt*x*(u[0]-u[-1])*(s_t)/(Delta*s);    
  }
}

static void da (scalar * ul, double t, scalar * kl)
{
  fprintf(stderr,"%g\n",s);
  scalar u = ul[0], k = kl[0];
  foreach() {
    k[] = u[] + dt*D_a*(u[1] + u[-1] -2*u[])/sq(Delta*(1-s)) - dt*x*(u[0]-u[-1])*(s_t)/(Delta*(1-s));
  }
}



int main()
{
  init_grid (64);
  size (1);
  TOLERANCE = 1e-4;

  // Pointer of the file to save stats
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  }

  

  run();

  fclose(fp_stats);
}
event init (i = 0)
{
  
  s = 0.5;
  dt = 0.001;
  foreach() {
    c[] = c_0; 
    a[] = a_0;
  }
	c[left] = neumann(0.);
	a[left] = neumann(0.);

}

event runninglog (i++)
{
  //output_ppm (C1, linear = true, spread = 2, file = "f.mp4", n = 200);
  //fprintf (stderr, "%d %g %g %g\n", i, t, dt,s);
  if (s<=EPS || s>1-EPS){
    return 1;
  }
   
}

event logstats (t += 0.1, t<t_end) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}

event logconcentration (t += 0.002){
    int j;
    double c_conc = 0;
    double a_conc = 0;
    double x_step = 0;
    char str[80];

    sprintf(str, "data_output/time=%.3f.dat", t);
    FILE * fp = fopen (str, "w");

    if (fp == NULL) {
      fprintf(stderr,"Error opening file!\n");
      exit(1);
    }

    for (j=0;  j<resolution; j++) {
      c_conc = interpolate(c, x_step);
      a_conc = interpolate(a, x_step);
      // And write them to file
      fprintf (fp, "%g %g %g %g %g\n", x_step, c_conc, a_conc, s, s_t);
      
      x_step += (1./100.);
      
    } 
    fclose(fp);
}



/*
event movie(i+=100){
  char timestring[100];

  view (fov = 27.0, tx = -0.5, ty = -0.4, width = 2000, height = 2000);
  clear();
  squares ("a", map = cool_warm, min = 0.0, max = 1.0);
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("cleanser concentration.mp4");
}
*/


event final (t = end)
{
  //char name[80];
  printf("#DONE at t=%g",t);
  //output_ppm (C1, file = name, n = 200, linear = true, spread = 2);
}
event integration (i++)
{
  dt = dtnext (0.00001);
  double c_right = 0;
  double a_right = 0;
  foreach_boundary (right){
    c_right = c[];
    a_right = a[];
  }
  s_t = -ALPHA*Q*M*pow(c_right,ALPHA)*pow(a_right,BETA);
  c[right] = c[] + Delta*s*s_t*((-1/M)+c[]);
  
  a[right] = a[] + Delta*(1-s)*s_t*((BETA/(M*ALPHA))+a[])/D_a;
  fprintf(stderr,"%g\n",s);
  runge_kutta ({c}, t, dt, dc, 4);
  runge_kutta ({a}, t, dt, da, 4);
  s = s + dt*s_t;
}