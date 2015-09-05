/* Outputs analytic solutions to the Hard Sphere Fokker Planck equation
 * with power law cooling coefficients.
 * as in Park & Petrosian 1995a eq 71,72, but from Steinacker et al 88
 *
 * Compilation: 
 * gcc make_hard_sphere_solutions.c -std=c99 -I$MYPREFIX/include -L$MYPREFIX/lib -lm -lgsl -lgslcblas  -o steady_state
 * */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_bessel.h>

#define MODE 5 // Input, change 

#if MODE == 0
#define fout "./hard_sphere_0/solution.txt"
#define r 2.0
#define B -1.0
#elif MODE == 1 
#define fout "./hard_sphere_1/solution.txt"
#define r 0.0
#define B 1.0
#elif MODE == 2 
#define fout "./hard_sphere_2/solution.txt"
#define r 2.0
#define B 1.0
#elif MODE == 3 
#define fout "./hard_sphere_3/solution.txt"
#define r 0.0
#define B -1.0
#elif MODE == 4 
#define fout "./hard_sphere_4/solution.txt"
#define q 2.0
#define s 1.0
#elif MODE == 5 
#define fout "./hard_sphere_5/solution_0.3.txt" // change tau and fout here
#define tau 0.3
#define q 3.0
#define s 0.0
#endif

#define nbins 1000
#define xmin 1e-3
#define xmax 1e3

#if MODE <=3
int main(int argc, char *argv[])
{
    const double lambda = 0;  // these are constant 
    const double x0 = 0.1;
    const double D = 1;
    const double A = 1;
    const double T = 1;

    FILE *fp = NULL;

    double *result = malloc(nbins * sizeof(*result));
    if (result == NULL)
        exit(-1);

    double *xarr = malloc(nbins * sizeof(*xarr));
    if (xarr == NULL)
        exit(-1);

    double xstep = log10(xmax/xmin)/(nbins-1);
    for (int i=0; i<nbins; i++)   
        xarr[i] = xmin * pow(10,i * xstep);

    double a = A / D;
    
    double b = B / D;
    
    double alpha = r - 1;
    
    double beta = -b/alpha;
    
    double theta = 1/(D*T);
    
    double lambda0 = 0.25*pow(alpha+3,2) + theta;

    double mu = sqrt(lambda0-lambda);

    double deltp = 0.5 * (a + 1) + mu;
    
    double deltm = 0.5 * (a + 1) - mu;
    
    double delta = 0, ddelta = 0;
    
    if (alpha > 0 ) {
        
        delta = deltp;
        
        ddelta = deltm;
    } else {

        delta = deltm;
        
        ddelta = deltp;
    }

    double gamma = (1 - delta)/alpha;

    double aa = 1 + delta/alpha;
    double bb = 1 + (delta - ddelta)/alpha;

    double y0 = fabs(beta) * pow(x0,alpha);
    
    double xc = pow(fabs(a)/fabs(b),1/alpha);

    printf("alpha=%g beta=%g gamma=%g \nxc = %g \n "
            ,alpha, beta, gamma, xc);

    
    for (int i=0; i<nbins; i++) { // Park & Petrosian 1995  

        double y = fabs(beta) * pow(xarr[i],alpha);
    
        if (beta > 0 && y < y0){    

            result[i] = 1/(D*fabs(b)) * pow(x0,-r) * pow(y/y0,-gamma)
                        * pow(y0,bb) * exp(-y0)
                        * gsl_sf_gamma(aa) / gsl_sf_gamma(bb)
                        * gsl_sf_hyperg_1F1(aa,bb,y)
                        * gsl_sf_hyperg_U(aa,bb,y0);
        } else if (beta > 0 && y > y0) {  

            result[i] = 1/(D*fabs(b)) * pow(x0,-r) * pow(y/y0,-gamma)
                        * pow(y0,bb) * exp(-y0)
                        * gsl_sf_gamma(aa) / gsl_sf_gamma(bb)
                        * gsl_sf_hyperg_1F1(aa,bb,y0)
                        * gsl_sf_hyperg_U(aa,bb,y);
        } else if (beta < 0 && y < y0) {  

            result[i] = 1/(D*fabs(b)) * pow(x0,-r) * pow(y/y0,-gamma)
                        * pow(y0,bb) * exp(-y)
                        * gsl_sf_gamma(bb-aa) / gsl_sf_gamma(bb)
                        * gsl_sf_hyperg_1F1(bb-aa,bb,y)
                        * gsl_sf_hyperg_U(bb-aa,bb,y0);
            
        } else if (beta < 0 && y > y0) {  

            result[i] = 1/(D*fabs(b)) * pow(x0,-r) * pow(y/y0,-gamma)
                        * pow(y0,bb) * exp(-y)
                        * gsl_sf_gamma(bb-aa) / gsl_sf_gamma(bb)
                        * gsl_sf_hyperg_1F1(bb-aa,bb,y0)
                        * gsl_sf_hyperg_U(bb-aa,bb,y);
        } 
    }

    printf("Output to %s \n", fout);

    if (!(fp = fopen(fout, "w")))
        printf("Can't open <%s> for writing ! \n",fout);
    
    for (int i=0; i<nbins; i++)
        fprintf(fp, "%g %g \n",xarr[i], result[i]);
   
    fclose(fp);
    free(xarr); 
    free(result);

    return EXIT_SUCCESS;
} 
#endif

#if MODE == 4
int main(int argc, char *argv[])
{
    const double x0 = 0.1;
    const double D = 1;
    const double A = 1;
    const double T = 1;
    const double theta = 1/(D*T);

    FILE *fp = NULL;

    double *result = malloc(nbins * sizeof(*result));
    if (result == NULL)
        exit(-1);

    double *xarr = malloc(nbins * sizeof(*xarr));
    if (xarr == NULL)
        exit(-1);

    double xstep = log10(xmax/xmin)/(nbins-1);
    for (int i=0; i<nbins; i++)   
        xarr[i] = xmin * pow(10,i * xstep);

    if (q == 2-s) {
        printf("Solution not valid here ! \n");
        exit(-1);
    }

    double a = A / D;
    
    double alpha = 0.5 * (2 -q - s);

    double beta = sqrt(theta)/fabs(alpha);

    double gamma = (q-1-a) / (2*alpha);
    
      printf("alpha=%g beta=%g gamma=%g \nx0 = %g \n "
            ,alpha, beta, gamma, x0);
    
    for (int i=0; i<nbins; i++) { // Park & Petrosian 1995  
        
        double y = beta * pow(xarr[i],alpha);

        double y0 = beta * pow(x0,alpha);

        double nu = fabs((q-1+a)/(2*alpha));

        result[i] = 1.0/D * 1/fabs(alpha) * pow(xarr[i], 0.5*(1-q+a)) 
            * pow(x0, 0.5*(1-q-a));

        if (y < y0) 
            result[i] *= gsl_sf_bessel_Inu (nu, y) 
                       * gsl_sf_bessel_Knu (nu, y0);    
        else 
            result[i] *= gsl_sf_bessel_Inu (nu, y0) 
                       * gsl_sf_bessel_Knu (nu, y);    
    }

    printf("Output to %s \n", fout);

    if (!(fp = fopen(fout, "w")))
        printf("Can't open <%s> for writing ! \n",fout);
    
    for (int i=0; i<nbins; i++)
        fprintf(fp, "%g %g \n",xarr[i], result[i]);
   
    fclose(fp);
    free(xarr); 
    free(result);

    return EXIT_SUCCESS;
} 
#endif

#if MODE == 5
int main(int argc, char *argv[])
{
    const double x0 = 0.1;
    const double D = 1;
    const double A = 1;
    const double T = 1;
    const double theta = 1/(D*T);

    FILE *fp = NULL;

    double *result = malloc(nbins * sizeof(*result));
    if (result == NULL)
        exit(-1);

    double *xarr = malloc(nbins * sizeof(*xarr));
    if (xarr == NULL)
        exit(-1);

    double xstep = log10(xmax/xmin)/(nbins-1);
    for (int i=0; i<nbins; i++)   
        xarr[i] = xmin * pow(10,i * xstep);

    if (q == 2-s) {
        printf("Solution not valid here ! \n");
        exit(-1);
    }

    double a = A / D;
    
    double alpha = 0.5 * (2 -q - s);

    double beta = sqrt(theta)/fabs(alpha);

    double gamma = (q-1-a) / (2*alpha);
    
      printf("alpha=%g beta=%g gamma=%g \nx0 = %g \n "
            ,alpha, beta, gamma, x0);
    
    for (int i=0; i<nbins; i++) { // Park & Petrosian 1995  
        
        double y = beta * pow(xarr[i],alpha);

        double y0 = beta * pow(x0,alpha);

        double nu = fabs((q-1+a)/(2*alpha));

        result[i] = 1/fabs(alpha) /(2*tau) * pow(xarr[i], 0.5*(1-q+a)) 
            * pow(x0, 0.5*(1-q-a));

        double arg = pow(xarr[i]*x0,alpha)/(2*tau*alpha*alpha);
        double e_arg = (pow(xarr[i],2*alpha)+pow(x0,2*alpha))
            /(4*tau*alpha*alpha);

        result[i] *= gsl_sf_bessel_Inu(nu, arg) * exp(-e_arg)*exp(-theta*tau);
    }

    printf("Output to %s \n", fout);

    if (!(fp = fopen(fout, "w")))
        printf("Can't open <%s> for writing ! \n",fout);
    
    for (int i=0; i<nbins; i++)
        fprintf(fp, "%g %g \n",xarr[i], result[i]);
   
    fclose(fp);
    free(xarr); 
    free(result);

    return EXIT_SUCCESS;
} 
#endif
