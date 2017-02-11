#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>


/* Implémentation de la méthode de Kuo et Mallik basée sur de
   l'adaptive rejection sampling
 */
double logPrior(double *theta, double sd, int p, int type, SEXP listevar)
{
    int i, pga;
    double pri;
    pga = length(listevar);
    
    pri = 0.0;
    for (i = pga; i < (pga+p+1); i++) {
	pri = pri + dnorm(theta[i], 0.0, sd, 1);
    }
    if (type == 3) {
	pri = pri + dgamma(1.0/theta[pga+p+1], 0.01, 100, 1);
    }
    return(pri);
}


double logVrais(double *theta, SEXP y, SEXP X, int type, SEXP listevar)
{
    int i, j, k, l, n, p, pga;
    double pred, like, *yr;
    
    n = length(y);
    p = length(X);
    pga = length(listevar);
    yr = REAL(y);
    
    
    like = 0.0;
    l=0;
    for (i = 0; i < n; i++) {
	
	pred = theta[pga]; /* intercept */
	
	for (j = 0; j < pga; j++) {
	    for (k = 0; k < length(VECTOR_ELT(listevar, j)); k++) {
		l = INTEGER(VECTOR_ELT(listevar, j))[k];
		pred = pred + REAL(VECTOR_ELT(X, l))[i] * 
		    theta[j] * 
		    theta[pga+1+l];
	    }
	}
	if (type == 1) { /* régression logistique */
	    like = like + (yr[i] * pred) - log(1 + exp(pred));
	}
	if (type == 2) { /* régression Poisson */
	    like = like + dpois(yr[i], exp(pred), 1);
	}
	if (type == 3) { /* régression Gaussienne  */
	    like = like + dnorm(yr[i], pred, sqrt(theta[pga+p+1]), 1);
	}
    }
    return(like);    
}


double logPoste(double *theta, double sd, SEXP y, SEXP X, int type, SEXP listevar)
{
    double pos;
    pos = logPrior(theta, sd, length(X), type, listevar) + 
	logVrais(theta, y, X, type, listevar);
    return(pos);
}





void Compteur(double k)
{
    Rprintf("%i\%\r", ((int) round(100*k)));
//    fflush(stdout);
    R_CheckUserInterrupt();
}


/* Doubling procedure décrite par Neal 2000 (fig. 3).
*/

SEXP doublingproc(double x0, int pp, double ys, int i,
		  double *theta, SEXP y, SEXP X, double sd, int type, SEXP listevar)
{
    
    double w, u2, lb, ub, flb, fub, v;
    int k2, ok, p;
    SEXP ulb;
    
    w = 100; /* largeur de la tranche */
    u2 = 0.0;
    p = length(X);
    
    PROTECT(ulb = allocVector(REALSXP,2));
    
    /*  Selon les recommandations de Neal 2003 (discussion), 
	on tire au sort la largeur de l'intervalle w
	dans une distribution avec une queue infinie  */
    GetRNGstate();
    w = rgamma(2, 2);
    u2 = unif_rand();
    PutRNGstate();
    
    /* On positionne l'intevalle aléatoirement autour de x0 */
    lb = x0 - w*u2;
    ub = lb+w;
    k2 = pp;

    /* cas où on bosse sur la variance résiduelle (type == 3, i.e. cas gaussien) */
    if (i == (p+p+1)) {
	/* cas où lb <0 */
	if (lb < 0.0000000001) {
	    lb = 0.0000000001;
	    ub = w;
	}
    }
    
    /* Puis on va augmenter w jusqu'à ce que les deux bornes 
       soient en dehors de la tranche */
    ok = 0;
    
    while (!ok) {	
	/* calcul de la posterior à la borne inf */
	theta[i] = lb;
	flb = logPoste(theta, sd, y, X, type, listevar);
	theta[i] = ub;
	fub = logPoste(theta, sd, y, X, type, listevar);
	
	/* si les deux bornes sont en dehors de la tranche, c'est ok */
	if ((k2>0)&&((ys < flb)||(ys < fub))) {
	    ok = 0;
	} else {
	    ok = 1;
	}
	if (!ok) {
	    k2--;
	    GetRNGstate();
	    v = unif_rand();
	    PutRNGstate();
	    
	    if (i == (p+p+1)) {
		if ((v < 0.5)&&(lb>0.000001)) {
		    lb = lb - (ub - lb);
		} else {
		    ub = ub + (ub - lb);
		}
	    } else {
		if (v < 0.5) {
		    lb = lb - (ub - lb);
		} else {
		    ub = ub + (ub - lb);
		}
	    }
	    
	}
    }

    REAL(ulb)[0] = lb;
    REAL(ulb)[1] = ub;
    UNPROTECT(1);
    return(ulb);
}



SEXP slisamplcoeur(double ub, double lb, double ys, double x0, int i,
		   SEXP y, SEXP X, double *theta, double sd, int type,
		   SEXP listevar)
{
    int ok, k2;
    double u, prop;
    SEXP ulb;
    PROTECT(ulb = allocVector(REALSXP,2));
    
    ok = 0;
    k2 = 0;
    
    while (!ok) {
	GetRNGstate();
	u = unif_rand();
	PutRNGstate();
	u = lb + u*(ub-lb);
	theta[i] = u;
	prop = logPoste(theta, sd, y, X, type, listevar);
	
	if (ys < prop) {
	    ok = 1;
	} else {
	    /* Une étape de shrinkage */
	    if (u < x0) {
		lb = u;
	    } else {
		ub = u;
	    }
	}
    }
    REAL(ulb)[0] = lb;
    REAL(ulb)[1] = ub;
    UNPROTECT(1);
    return(ulb);
}






/* Slice sampling */
/* Mise en oeuvre du gibbs sampling. Mêmes arguments que la posterior + 2:
   - tabr: le tableau destiné à contenir les valeurs de paramètres simulés
   - nsim: le nombre de simulations
   Remarque: theta contient les valeurs de départ
 */

SEXP GibbsSampling(SEXP y, SEXP X, double *theta,
		   int nsim, double sd, int type, SEXP listevar)
{
    int i, k, p, pp, ma, pga;
    double x0, ys, u, v, lb, ub, prop;
    SEXP resu, ulb;
    
    v = 0.0;
    pp = 5;
    ys = 0.0;
    lb=0.0;
    ub=0.0;
    p = length(X);
    
    /* listevar contient les "regroupements de variables", 
       i.e. les ensembles de variables pour lesquels le coefficient
       gamma est égal à 0. 
       C'est une liste contenant les numéros de variable (numérotés à partir de 0)
       associés
    */
    /* Le nombre de coefficients gamma est donc: */
    pga = length(listevar);
    
    /* Nombre total de coefficients: pga gamma + p beta + 1 intercept 
       (+ une variance si gaussien) */
    ma = pga+1+p;
    if (type==3)
	ma = pga+2+p;

    PROTECT(resu = allocVector(REALSXP, nsim * ma));
    
    
    /* au début, les valeurs initiales sont dans theta */
    Rprintf("Itération              ");

    for (k = 0; k < nsim; k++) { 
	
	/* Affichage du numéro d'itération */
	Compteur(((double)k)/((double) (nsim-1)));
	
	/* Boucle pour chaque paramètre */	
	for (i = 0; i < ma; i++) { 
	    
	    /* tirage au sort d'une valeur comprise entre 0 et f(x_0), 
	       qui déterminera la tranche */
	    x0 = theta[i];
	    ys = logPoste(theta, sd, y, X, type, listevar);
	    
	    GetRNGstate();
	    u = unif_rand();
	    PutRNGstate();
	    ys = ys + log(u);
	    
	    /* Premier cas: le coefficient est un gamma */
	    if (i < pga) {		
		/* Cas zero */
		theta[i] = 0;
		lb = logPoste(theta, sd, y, X, type, listevar);

		/* Cas 1 */
		theta[i] = 1;
		ub = logPoste(theta, sd, y, X, type, listevar);

		/* standardisation sur les possibles */
		prop=exp(lb)+exp(ub);
		lb=exp(lb)/prop;
		ub=exp(ub)/prop;
		
		/* Et tirage au sort */
		if (u<lb) {
		    theta[i] = 0;
		} else {
		    theta[i] = 1;
		}
		
		
	    } else {
		/* Sinon slice sampling classique */
		PROTECT(ulb = doublingproc(x0, pp, ys, i, theta, y, X, sd, type, 
					   listevar));
		lb = REAL(ulb)[0];
		ub = REAL(ulb)[1];
		UNPROTECT(1);
		
		/*  On tire au sort dans la tranche telle que ys < f(x_1) */
		PROTECT(ulb = slisamplcoeur(ub, lb, ys, x0, i,
					    y, X, theta, sd, type, listevar));
		lb = REAL(ulb)[0];
		ub = REAL(ulb)[1];
		UNPROTECT(1);
	    } 
	    
	    REAL(resu)[k+(i*nsim)] = theta[i];
	}
    }
    UNPROTECT(1);
    Rprintf("\n");
    return(resu);
}




/* On suppose que X ne contient pas l'intercept */
SEXP kuomallik(SEXP y, SEXP X, SEXP nsimr, SEXP thetar, SEXP sdr, SEXP typer, 
	       SEXP listevar)
{
    /* déclaration des variables */
    int nsim, p, type;
    double sd, *theta;
    SEXP nsimrv, name, resu, typev;
    
    /* Variables utiles */
    p = length(X); /* nombre de variables expli */
    PROTECT(name = getAttrib(X, R_NamesSymbol)); /* nom de ces variables*/
    sd = REAL(sdr)[0];
    PROTECT(nsimrv = coerceVector(nsimr, INTSXP));
    nsim = INTEGER(nsimrv)[0];
    theta = REAL(thetar);
    PROTECT(typev = coerceVector(typer, INTSXP));
    type = INTEGER(typev)[0];
    
    PROTECT(resu = GibbsSampling(y, X, theta,
				 nsim, sd, type, listevar));        
    UNPROTECT(4);
    return(resu);
}





/* *******************************************************************************
   The functions below implement the one-dimensional projection pursuit approach
   ******************************************************************************* */



/* Solid angle transform. Takes as argument a vector of length n-1 corresponding to a point
   in a hypercube with side [0,1], and return a normed vector (squared size = 1) in a 
   space of dimension n (i.e. a point on a hypersphere with radius 1 in dimension n).
   From Friedman and Steppel 1973, after correction, since the formulas given by these
   authors for the odd case are incorrect. 

   Output: a SEXP containing a vector
*/
SEXP SAT(SEXP etas)
{ 
    int n, i, j, i0, j0, ev;
    SEXP Xs;
    double *eta, *X, pr;

    /* Length of the resulting vector */
    n=length(etas)+1;

    /* We store the results in C tables */
    PROTECT(Xs = allocVector(REALSXP, n));
    X = REAL(Xs);
    eta = REAL(etas);
    
    /* initialization of X */
    for (i = 0; i < n; i++) {
	X[i] = 0.0;
    }
    
    /* is n odd or even? */
    ev = 1 - (n%2);
    
    
    if (ev) {
	/* If n is even */
        for (i0 = 0; i0 < (n/2-1); i0++) {
	    i = i0+1;
            pr = 1.0;
            if (i>1) {
                for (j0 = 0; j0 < (i-1); j0++) {
		    j = j0+1;
                    pr = pr*R_pow(eta[(2*j)-1], (1.0/(((double) (n-2*j)))));
		}
	    }
	    pr = pr*cos(asin(R_pow(eta[(2*i)-1],(1.0/(((double) (n-2*i)))))))*
		sin(2.0*M_PI*eta[2*i-2]);  
	    X[2*i-1] = pr;
	}
	pr = 1.0;
        for (j0 = 0; j0 < (n/2-1); j0++) {
	    j = j0+1;
            pr = pr*R_pow(eta[2*j-1],(1.0/((double) (n-2*j))));
	}
        X[n-1] = pr*sin(2.0*M_PI*eta[n-2]);
	for (i0 = 0; i0 < (n/2); i0++) {
	    i = i0+1;
            X[2*i-2] = X[2*i-1]/tan(2.0*M_PI*eta[2*i-2]);
        }
	
    } else {
        
	/* If n is odd */
	pr=1.0;
        for (j0 = 0; j0 < ((n-3)/2); j0++) {
	    j=j0+1;
            pr = pr*R_pow(eta[2*j-1],(1.0/((double) (n-2*j))));
        }
        X[n-3] = pr*sqrt(4.0*eta[n-3]-(4.0*R_pow(eta[n-3],2.0)))*cos(2.0*M_PI*eta[n-2]);
        X[n-2] = pr*sqrt(4.0*eta[n-3]-(4.0*R_pow(eta[n-3],2.0)))*sin(2.0*M_PI*eta[n-2]);
        X[n-1] = pr*(1.0-2.0*eta[n-3]);

        for (i0 = 0; i0 < ((n-3)/2); i0++) {
	    i = i0+1;
            pr = 1.0;
            if (i>1) {
                for (j0 = 0; j0 < (i-1); j0++) {
		    j=j0+1;
                    pr = pr*R_pow(eta[2*j-1],(1.0/((double) (n-2*j))));
		}
            }
            X[2*i-1] = pr*cos(asin(R_pow(eta[2*i-1],(1/((double) (n-2*i))))))*
		sin(2.0*M_PI*eta[2*i-2]);
        }
        for (i0 = 0; i0< ((n-3)/2); i0++) {
	    i = i0+1;
            X[2*i-2] = X[2*i-1]/tan(2.0*M_PI*eta[2*i-2]);
	}
    }
    
    UNPROTECT(1);
    return(Xs);
}


/* Inverse transform of the solid angle transform: Xs is a vector of length n corresponding
   to a point on a hypersphere with radius 1, and returns a vector corresponding to
   a point in a hypercube with side [0,1].
   Note that this transform may return "not a number" if at least one of the values 
   is exactly 0 (rare).

   Output: a SEXP containing a vector
*/
SEXP invSAT(SEXP Xs)
{
    int n, i, i0, ev, j0, j;
    SEXP etas;
    double *eta, *X, pr, z, uu;

    n = length(Xs);
    uu=0.0;
    j=1;
    j0=1;

    /* We store the results in C tables */
    PROTECT(etas = allocVector(REALSXP, n-1));
    X = REAL(Xs);
    eta = REAL(etas);

    for (i = 0; i < n-1; i++) {
	eta[i] = 0.0;
    }
    
    /* is n odd or even? */
    ev = 1 - (n%2);

    if (ev) {
	/* even n */
	for (i0 = 0; i0 <(n/2); i0++) {
	    i = i0+1;
            pr = atan2(X[2*i-1], X[2*i-2]);
	    if (pr < 0)
		pr = pr + 2.0*M_PI;
            eta[2*i-2] = pr/(2.0*M_PI);
        }
        z = 1.0;
	
        for (i0 = 0; i0 < ((n/2)-1); i0++) {
	    i = i0+1;
            if (i>1) {
                z = z * R_pow(eta[2*(i-1)-1], (1.0/((double) (n-2*(i-1)))));
	    }
            eta[2*i-1] = R_pow((sin(acos(X[2*i-1]/(z*(sin(2.0*M_PI*eta[2*i-2])))))), (n-2*i));
        }
    } else {
	/* odd n */
        for (i0 = 0; i0 < ((n-3)/2); i0 ++) {
	    i = i0+1;
            pr = atan2(X[2*i-1], X[2*i-2]);
	    if (pr < 0)
		pr = pr + 2.0*M_PI;
            eta[2*i-2] = pr/(2.0*M_PI);
        }
        z = 1.0;
	
        for (i0 = 0; i0 < (((n-3)/2)); i0++) {
	    i = i0+1;
            if (i>1) {
                z = z*R_pow(eta[2*(i-1)-1], (1.0/((double) (n-2*(i-1)))));
	    }
            eta[2*i-1] = R_pow(sin(acos(X[2*i-1]/(z*(sin(2.0*M_PI*eta[2*i-2]))))), 
			       ((double) (n-2*i)));
        }

	pr=1.0;
        for (j0 = 0; j0 < ((n-3)/2); j0++) {
	    j=j0+1;
            pr = pr*R_pow(eta[2*j-1],(1.0/((double) (n-2*j))));
        }

	eta[n-3] = (1.0 - X[n-1]/(pr))/2.0;
        uu = pr*sqrt((4.0*eta[n-3])-(4.0*R_pow(eta[n-3],2.0)));
        uu = atan2(X[n-2]/uu, X[n-3]/uu);
	if (uu < 0)
	    uu = uu + 2.0*M_PI;
	eta[n-2] = uu/(2.0*M_PI);

    }
    UNPROTECT(1);
    return(etas);
}




/* Converts a vector with values comprised between -infinity and +infinity to a vector where
   each value is comprised in the range [0,1] (using logit, if versO1 =1) and conversely 
   (using inverse logit, if vers01 = 0).

   Output: a SEXP containing a vector 
 */
SEXP logitVector(SEXP vec, int vers01)
{
    SEXP vecso;
    int n, i;
    
    n = length(vec);
    PROTECT(vecso = allocVector(REALSXP, n));
    if (vers01) {
	for (i = 0; i < n; i++) {
	    REAL(vecso)[i] = exp(REAL(vec)[i])/(1.0 + exp(REAL(vec)[i]));		
	}
    } else {
	for (i = 0; i < n; i++) {
	    REAL(vecso)[i] = log(REAL(vec)[i]/(1.0 - REAL(vec)[i]));		
	}
    }
    UNPROTECT(1);
    return(vecso);
}


/* Converts a point on a unit sphere in a n-dimensional space to a point in an
   unbounded Euclidean (n-1) dimensional space: this function combines the functions invSAT
   and logitVector (with vers01=0).

   Output: a SEXP containing a vector 
*/
SEXP axeVersEuc(SEXP x)
{
    SEXP xHypercube, xEuc;
    PROTECT(xHypercube = invSAT(x));
    PROTECT(xEuc = logitVector(xHypercube, 0));
    UNPROTECT(2);
    return(xEuc);
}

/* Converts a point in an in an unbounded Euclidean (n-1) dimensional space to a 
   point on a unit sphere in a n-dimensional space: this function combines the 
   functions logitVector (with vers01=0) and invSAT.

   Output: a SEXP containing a vector
*/
SEXP eucVersAxe(SEXP x)
{
    SEXP xHypercube, axe;
    PROTECT(xHypercube = logitVector(x, 1));
    PROTECT(axe = SAT(xHypercube));
    UNPROTECT(2);
    return(axe);
}


/*  This structure is used for the optimization. It is passed as the *ex object
    containing the parameters used by the function calculating the criterion.
    This structure contains the data.frame explored by projection pursuit, 
    the R function calculating the criterion, the list of parameters passed 
    to this function, and the environment in which this function will be used.
 */
typedef struct optimpp{
    SEXP df;
    SEXP fonc;
    SEXP pare;
    SEXP envir;
} optimpp, *OPTPP;



/* Calculates the covariance matrix associated to a 
   data.frame.
   
   output: a SEXP containing the matrix stored as a vector.
 */
SEXP matcovar(SEXP df)
{
    int n, p, i1, i2, j;
    SEXP mat;
    double *matr, *v1, *v2;
    
    p = length(df);
    n = length(VECTOR_ELT(df, 0));
    PROTECT(mat = allocVector(REALSXP, p*p));
    matr=REAL(mat);
    
    for (i1 = 0; i1 < p; i1++) {
	v1 = REAL(VECTOR_ELT(df, i1));
	for (i2 = 0; i2 < p; i2++) {
	    v2 = REAL(VECTOR_ELT(df, i2));	    
	    matr[i1+p*i2] = 0.0;
	    for (j = 0; j < n; j++) {
		matr[i1+p*i2] += (v1[j] * v2[j]);
	    }
	    matr[i1+p*i2] = matr[i1+p*i2]/((double) n);
	}
    }
    UNPROTECT(1);
    return(mat);
}



/* Wrapper function to find the eigenvectors and eigenvalues of
   a matrix. Use the LAPACK function dsyev.
   
   output: a SEXP corresponding to a list containing:
   (i) a SEXP storing the eigenvalues as a vector
   (ii) a SEXP storing a list of SEXP, each one being an 
   eigenvector.
  */
SEXP eigendecompo(SEXP mat)
{
    int m, info, i, j, nvp, k;
    SEXP valp, work, dfso, matc, valpfin, vecpfin, vecptmp;
    double wo, *workr;

    /* Prepare the data for the use of dsyev */
    const char jobz = 'V'; /* We need the eigenvectors *and* eigenvalues (otherwise, 
			      would be 'N' for values only */
    const char uplo = 'U'; /* The upper triangle of the symetric mat will be passed to the 
			      function (and not the lower -- it would be 'L' */
    const int n = length(mat); /* the matrix is n = pxp */
    const int p = (int) (sqrt((double) n)+0.2); /* p is the order of the matrix 
						   (number of cols). We use +0.2 to 
						   avoid rounding errors  */
    PROTECT(matc = allocVector(REALSXP, n)); /* Will store the matrix to be analysed */
    int lwork = -1;  /*  dimension of the workspace. According to the help of lapack:
			 LWORK >= max(1,3*N-1). We will use 3*n */
    lwork = 3*n;
    info = 0; /* For the output: checks whether successful exit (=0) */
    wo = 0;
    
    PROTECT(valp = allocVector(REALSXP, n)); /* output: will be used to store the eigenvalues */
    PROTECT(work = allocVector(REALSXP, lwork)); /* matrix used as workspace by dsyev
						    dimension lwork (we use lwork=3n)
						 */
    
    /* initialize valp and copy the values of mat in matc */
    for (i = 0; i < n; i++) {
	REAL(valp)[i]=0;
	REAL(matc)[i] = REAL(mat)[i];
    }

    F77_NAME(dsyev)(&jobz, &uplo,
		    &p, REAL(matc), &p, REAL(valp), 
		    REAL(work), &lwork, &info);
    

    /* Checks if the algorithm converged */
    if (info < 0)
	error("Illegal value in the calculation of eigenvalues");
    if (info > 0)
	error("The algorithm used to calculate the eigenvalues did not converge");
    
    
    
    /* How many positive eigenvalues ? */
    nvp = 0;
    for (i = 0; i < p; i++) {
	if (REAL(valp)[i] > 0.000000001) {
	    nvp++;
	}
    }

    /* Stores the results in the suitable vectors and data.frames */
    PROTECT(valpfin = allocVector(REALSXP, nvp));
    PROTECT(vecpfin = allocVector(VECSXP, nvp));
    k=0;
    for (i = 0; i < p; i++) {
	if (REAL(valp)[i] > 0.000000001) {
	    REAL(valpfin)[nvp-1-k] = REAL(valp)[i];
	    PROTECT(vecptmp = allocVector(REALSXP, p));
	    for (j = 0; j < p; j++) {
		REAL(vecptmp)[j] = REAL(matc)[i*p+j];
	    }
	    SET_VECTOR_ELT(vecpfin, (nvp-1-k), vecptmp);
	    UNPROTECT(1);
	    k++;
	} 
    }

    /* And stores everything in the output list */
    PROTECT(dfso = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dfso, 0, valpfin);
    SET_VECTOR_ELT(dfso, 1, vecpfin);
 
    /* output */
    UNPROTECT(6);
    return(dfso);    
}



/* Projects the data.frame df on the hyperplane orthogonal to the vector u
   (expected to be normed) and returns the result. More Precisely, let X be the
   matrix corresponding to the data.frame, this function calculates:
   Y = X(I - uu^t)
   where I is the nxn identity matrix (and n is the number of columns in X)
   
   Output: a SEXP storing the data.frame Y
*/
SEXP orthoProjDf(SEXP df, SEXP u)
{
    int i, j, k, n, nl;
    double tmp;
    SEXP uut, dfso, colo;

    n = length(u);
    nl = length(VECTOR_ELT(df, 0));
    PROTECT(uut = allocVector(REALSXP, n*n));
    PROTECT(dfso = allocVector(VECSXP, n));

    /* Calculation of uu^t */
    for (i = 0; i<n; i++) {
	for (j = 0; j<n; j++) {
	    REAL(uut)[j+n*i] = REAL(u)[i]*REAL(u)[j];
	}
    }
    
    /* Calculation of the data.frame */
    for (i = 0; i < n; i++) {
	PROTECT(colo = allocVector(REALSXP, nl));
	for (j = 0; j < nl; j++) {
	    tmp = 0.0;
	    for (k = 0; k< n; k++) {
		tmp += REAL(VECTOR_ELT(df, k))[j]*REAL(uut)[k+i*n];
	    }
	    REAL(colo)[j] = REAL(VECTOR_ELT(df, i))[j]-tmp;
	}
	SET_VECTOR_ELT(dfso, i, colo);
	UNPROTECT(1);
    }	
    UNPROTECT(2);
    return(dfso);
}




/* Évaluates the function fonc with the parameters stored in the list par in the environment
   envir, based on the vector storing the coordinates of the rows of the data.frame df
   on the direction pointed by the vector u. The vector u contains the coordinates of a point
   on a unit sphere in a n-dimensional space, and is calculated from the SEXP x (a vector in
   the (n-1)-dimensional Euclidean space).
   
   output: a double containing the value of the function
*/
double evalfonc(SEXP x, SEXP df, SEXP fonc, SEXP par, SEXP envir)
{
    int n, nl, i, j;
    double *xr, *dfi, *vecr, out;
    SEXP xSAT, vec, resu;

    if(!isEnvironment(envir)) error("'env' should be an environment");
    
    /* Transforming a vector from the dimension n-1 to a vector storing the 
       coordinates of a point on a unit sphere in a n-dimensional space */
    PROTECT(xSAT = eucVersAxe(x));
    
    /* Prepare the objects to store the output */ 
    n = length(xSAT);
    nl = length(VECTOR_ELT(df, 0));
    PROTECT(vec = allocVector(REALSXP, nl));
    xr = REAL(xSAT);   
    
    /* projection of df on x */
    vecr = REAL(vec);
    for (i = 0; i < nl; i++) {
	vecr[i] = 0.0;
	for (j = 0; j < n; j++) {
	    vecr[i] += xr[j] * REAL(VECTOR_ELT(df, j))[i];
	}
    }
    
    /* evaluation of the function on vec */
    defineVar(install("x"), vec, envir);
    defineVar(install("par"), par, envir);
    PROTECT(resu = eval(fonc, envir));
    out = REAL(resu)[0];    
    UNPROTECT(3);

    /* return the result */
    return(out);
}




/* Same function as evalfonc, but the vector xSAT is already a normed vector in a 
   n-dimensional space.

   output: a double containing the value of the function   
*/
double evalfoncSuraxe(SEXP xSAT, SEXP df, SEXP fonc, SEXP par, SEXP envir)
{
    int n, nl, i, j;
    double *xr, *dfi, *vecr, out;
    SEXP vec, resu;

    if(!isEnvironment(envir)) error("'env' should be an environment");
        
    /* Préparation of the structures to store the results */ 
    n = length(xSAT);
    nl = length(VECTOR_ELT(df, 0));
    PROTECT(vec = allocVector(REALSXP, nl));
    xr = REAL(xSAT);   
    
    /* projection on df on x */
    vecr = REAL(vec);
    for (i = 0; i < nl; i++) {
	vecr[i] = 0.0;
	for (j = 0; j < n; j++) {
	    vecr[i] += xr[j] * REAL(VECTOR_ELT(df, j))[i];
	}
    }
    
    /* evaluation of the function on vec */
    defineVar(install("x"), vec, envir);
    defineVar(install("par"), par, envir);
    PROTECT(resu = eval(fonc, envir));
    out = REAL(resu)[0];    
    UNPROTECT(2);

    return(out);
}




/* Function to be optimized when searching the value of the axis 
   that maximizes the criterion. This function has the correct arguments for use 
   with the minimization function nmmin(): 
   n is the dimension of the vector,
   par contains the starting value of the vector. This vector has length n. The function will 
   transform this unconstrained vector in a n-space into a unit-normed vector
   in a n-dimensional space
   ex contains the parameters associated to the function itself (i.e. the R function, 
   the data.frame, etc.). It is an "unclassed" object of class optimpp (see above).

   output: a double storing the value of the function at the current value of the vector par.
 */
static double optimfonc(int n, double *par, void *ex)
{
    SEXP x, x2, dforth, df, fonc, pare, envir;
    int i;
    OPTPP pal = (OPTPP) ex;
    double so;
    
    /* We get all the required elements from *ex */
    PROTECT(x = allocVector(REALSXP, n));
    df = pal->df;
    fonc = pal->fonc;
    pare = pal->pare;
    envir = pal->envir;
    
    /* We store the current value for x */
    for (i = 0; i < n; i++) {
	REAL(x)[i] = par[i];
    }
    
    /* And we evaluate the function */
    so = -evalfonc(x, df, fonc, pare, envir);

    UNPROTECT(1);

    return(so);
}




/* This function performs the optimization. The arguments are:
   - the starting value for the vector x (unconstrained vector in n-1 dimensional space)
   - the data.frame df
   - the R function fonc to be optimized
   - the list containing the parameters of the fonction
   - the environnement envir
   - the maximum number of iterations maxitr
   - whether the optimisation should be verbose (1) or silent (0)

   output: a SEXP containing the optimal vector in the unconstrained n-1 dimensional space
 */
SEXP optimisewrap(SEXP x, SEXP df, SEXP fonc, 
		  SEXP pare, SEXP envir, SEXP maxitr, SEXP verbo)
{
    int i, fail, fncount, maxit;
    OPTPP monobj = (OPTPP) malloc(sizeof(optimpp));
    double *xin, *xpr, Fmin, abstol, intol, su;
    SEXP xp, xlo, dfso, maxitv;
    
    /* checks that maxitr is actually an integer */
    PROTECT(maxitv = coerceVector(maxitr, INTSXP));
    maxit = INTEGER(maxitv)[0];
    UNPROTECT(1);
    
    /* the object optimpp used by the function optimfonc (argument *ex) */
    monobj->df = df;
    monobj->fonc = fonc;
    monobj->pare = pare;
    monobj->envir = envir;
    
    /* prepare the elements required by the function nmmin */
    PROTECT(xp = allocVector(REALSXP, length(x))); /* Will store the result of optimization */
    PROTECT(xlo = allocVector(REALSXP, length(x))); /* Will store the starting value of the 
						       vector x */
    for (i = 0; i < length(x); i++)
	REAL(xlo)[i] = REAL(x)[i];
    xin = REAL(xlo);
    xpr = REAL(xp);

    Fmin = 0.0; /* Final value of the function for the optimal axis */
    abstol = -1000000000000; /* absolute convergence tolerance. We do not care */
    intol = 1e-10; /* relative convergence tolerance */
    fncount = 0; /* output: number of iterations of the algorithm */
    fail = 0; /* result of the function: convergence or not? */
        
    /* Optimization with the Nelder-Mead algorithm */
    nmmin(length(x), xin, xpr, &Fmin, optimfonc,
	  &fail, abstol, intol, (void *)monobj,
	  1.0, 0.5, 2.0, 0,
	  &fncount, maxit);

    /* If the optimisation fails */
	if (fail!=0)
	    Rprintf("Minimization failed. \nCode returned by the function optim (Nelder-Mead algo): %i (cf. ?optim)\n", fail);
    
    /* Clean the workspace */
    free(monobj);
    UNPROTECT(2);
    
    return(xp);
}




/* The preoptimisation approach recommended by Friedman 1987, p. 56.
   Because, for a given criterion, there can be many local maximum, Friedman
   recommends to first perform an imprecise preliminary search of good candidate axes, 
   and then to find the vector that maximize the criterion using a 
   classical maximization approach. 
   The arguments are the following:
   - axe is the n-vector (note that it is expected that this vector is normed before the 
   use of the function). It is a point on the unit sphere in the n-dimensional space
   - df is the data.frame to be explored
   - fonc is the R function to calculate the criterion
   - par is the list of parameters to be passed to the R function
   - envir is the R environment in which the R function will be used.
   
   output: a SEXP storing the candidate n-vector
 */
SEXP preoptimFriedman1(SEXP axe, SEXP df, 
		       SEXP fonc, SEXP par, SEXP envir)
{
    double I0, *eir, fp, fm, *alphao, *alphar, f, s, tmp, fm1, su;
    int conti, i, nc, j;
    SEXP ei, alphac, alphalo;
    
    conti = 1;
    nc = length(df);

    /* local copy of the axis */
    PROTECT(alphac = allocVector(REALSXP, nc));
    PROTECT(alphalo = allocVector(REALSXP, nc));
    alphao = REAL(alphac);
    alphar = REAL(alphalo);
    
    for (i = 0; i < nc; i++) {
	alphar[i] = REAL(axe)[i];
	alphao[i] = alphar[i];
    }
    
    /* repeat until convergence */
    while (conti) {
	
	/* calculation of the function for the starting value */
	I0 = evalfoncSuraxe(alphalo, df, fonc, par, envir);
	fm1 = I0;
	
	/* Then, for each dimension */
	for (i = 0; i < nc; i++) {
	    
	    /* calculation of f+ */
	    for (j = 0; j < nc; j++) {
		alphao[j] = alphar[j];
	    }
	    alphao[i] = alphar[i] + 1;
	    for (j = 0; j < nc; j++) {
		alphao[j] = (alphao[j]/sqrt(2.0))/(sqrt(1.0 + alphar[i]));
	    }
	    fp = evalfoncSuraxe(alphac, df, fonc, par, envir);
	    
	    /* calculation of f- */
	    for (j = 0; j < nc; j++) {
		alphao[j] = alphar[j];
	    }
	    alphao[i] = alphar[i] - 1;
	    for (j = 0; j < nc; j++) {
		alphao[j] = (alphao[j]/sqrt(2.0))/(sqrt(1.0 - alphar[i]));
	    }
	    fm = evalfoncSuraxe(alphac, df, fonc, par, envir);
	    
	    /* selection of the direction leading to the best increase */
	    if (fp > fm) {
		f = fp;
		s = 1.0;
	    } else {
		f = fm;
		s = -1.0;
	    }
	    if (f > fm1) {
		tmp = alphar[i];
		alphar[i] = alphar[i] + s;
		for (j = 0; j < nc; j++) {
		    alphar[j] = (alphar[j]/sqrt(2.0))/(sqrt(1.0 +s*tmp));
		}		
		fm1 = f;
	    }
	}
	
	/* Check the convergence */
	f = evalfoncSuraxe(alphalo, df, fonc, par, envir);
	if (fabs(I0 - f) < 1e-8)
	    conti = 0;
	
	I0 = f;
    }

    /* On rajoute un bruit sur les valeurs nulles, pas bien gérées par la SAT */
    for (i = 0; i < length(alphalo); i++) {
	if (fabs(REAL(alphalo)[i])< 0.0001) {
	    REAL(alphalo)[i] = 0.0001;
	}
    }
    
    /* Restandardisation */
    su = 0.0;
    for (i = 0; i < length(alphalo); i++) {
	su += R_pow(REAL(alphalo)[i], 2.0);
    }
    for (i = 0; i < length(alphalo); i++) {
	REAL(alphalo)[i] = REAL(alphalo)[i]/sqrt(su);
    }

    
    UNPROTECT(2);
    return(alphalo);
}




/* This function combines the preoptimisation approach of Friedman and the optimisation 
   itself (with the Nelder-Mead algorithm). 
   This function has the same arguments as optimizewrap (except for x, which is 
   called axe here). 

   Note that the function preoptimFriedman1 only accepts a normed vector
   in a n-dimensional space (as the present function: axe is expected to be normed), 
   whereas optimisewrap only accepts an unconstrained vector
   in a (n-1)-dimensional space. This function use the functions axeVersEuc and eucVersAxe to
   transform to and from unconstrained space.

   This function is the core of the projection pursuit approach.

   output: a SEXP storing the coordinates of optimal axis
 */
SEXP trouvaxe(SEXP df, SEXP axe, SEXP fonc, SEXP par, 
	      SEXP envir, SEXP maxitr, SEXP verbo)
{
    SEXP axe2, versnmu, resu, axeres, axeresfin;
    int nc, nl, i, j;
      
    nc = length(df);
    nl = length(VECTOR_ELT(df, 0));
    
    /* Friedman preoptimisation approach */
    PROTECT(axe2 = preoptimFriedman1(axe, df, 
				     fonc, par, envir));
    
    
    /* The minimization. (before: from n to (n-1); after: from (n-1) to n)     
     */
    PROTECT(versnmu = axeVersEuc(axe2));
    PROTECT(axeres = optimisewrap(versnmu, df, fonc, par, envir, 
				  maxitr, verbo));
    
    PROTECT(axeresfin = eucVersAxe(axeres));
    
    UNPROTECT(4);
    return(axeresfin);
}



/* This function samples a point on a unit sphere in a n-dimensional space, where n is the number
   of column in the data.frame df passed as argument.

   output: a SEXP storing the coordinates of the point in a vector
*/
SEXP initrand(SEXP df)
{
    int nc, i, ok;
    double su;
    SEXP vec;

    nc = length(df);
    PROTECT(vec = allocVector(REALSXP, nc));
    
    /* sample a direction on a sphere. We therefore eliminate all vectors with
       a size > 1
     */
    ok = 0;
    GetRNGstate();
    while (!ok) {
	for (i = 0; i < nc; i++) {
	    REAL(vec)[i] = (unif_rand()*2.0)-1.0;
	}
	su = 0.0;
	for (i = 0; i < nc; i++) {
	    su += R_pow(REAL(vec)[i],2.0);
	}
	if (su < 1)
	    ok = 1;
    }
    PutRNGstate();

    /* Because the SAT transform does not manage the zero values, we add a small noise in the 
       unlikely case where one of the coordinates would be equal to 0 */
    for (i = 0; i < length(vec); i++) {
	if (fabs(REAL(vec)[i])< 0.0001) {
	    REAL(vec)[i] = 0.0001;
	}
    }
    
    /* standardization of the vector */
    su = 0.0;
    for (i = 0; i < nc; i++) {
	su += R_pow(REAL(vec)[i], 2.0);
    }
    for (i = 0; i < nc; i++) {
	REAL(vec)[i] = REAL(vec)[i]/sqrt(su);
    }


    
    UNPROTECT(1);
    return(vec);
}




/* This function selects the principal axis that maximizes the function fonc. Prior to the 
   projection pursuit, a possibility currently unused but programmed could be to select as
   a starting value for the axis the principal axis of the data.frame df for which the function 
   fonc is maximized (instead of a random vector as with initrand). The arguments are:
   - the data.frame df
   - the R function fonc
   - the list of parameters passed to the function fonc
   - the environment envir where this function will be used
   
   output: the SEXP storing the vector maximizing the criterion.
 */
SEXP initpca(SEXP df, SEXP fonc, SEXP par, SEXP envir)
{
    SEXP sigma, vpp, vecp, cr, vec, hyperc;
    double *crr, m1;
    int n, i, im1;
    
    /* We perform the PCA and we store the principal axes in vecp */
    PROTECT(sigma = matcovar(df));
    PROTECT(vpp = eigendecompo(sigma));
    PROTECT(vecp = VECTOR_ELT(vpp, 1));
    n = length(vecp);

    /* We calculate the value of the criterion for each principal axis 
       (after conversion from a sphere  in n dimension to an unconstrained
       (n-1)-dimensional space with the SAT transform)
     */
    PROTECT(cr = allocVector(REALSXP, n));
    crr = REAL(cr);
    for (i = 0; i < n; i++) {
	PROTECT(hyperc = axeVersEuc(VECTOR_ELT(vecp, i)));
	crr[i] = evalfonc(hyperc, df, fonc, par, envir);
	UNPROTECT(1);
    }
    
    /* Selection of the maximal value */
    m1=crr[0];
    im1=0;
    for (i = 1; i < n; i++) {
	if (crr[i] > m1) {
	    m1 = crr[i];
	    im1 = i;
	}
    }
    
    /* Because the SAT transform does not manage the zero values, we add a small noise in the 
       case where one of the coordinates would be equal to 0 */
    vec = VECTOR_ELT(vecp, im1);
    for (i = 0; i < length(vec); i++) {
	if (fabs(REAL(vec)[i])< 0.0001) {
	    REAL(vec)[i] = 0.0001;
	}
    }
    
    UNPROTECT(5);
    return(vec);
}





/* This function performs the matrix multiplication XU. Here, X and U 
   are actually data.frames stored in SEXP objects
   
   output: a SEXP storing the result
 */
SEXP XU(SEXP X, SEXP U)
{
    int nc1, nl1, nc2, i, j, k;
    SEXP resu, vectmp, uij;
    double *vtmp, *uijr;
    
    nc1 = length(X);
    nl1 = length(VECTOR_ELT(X,0));
    nc2 = length(U);
    
    PROTECT(resu = allocVector(VECSXP, nc2));
    
    for (j = 0; j < nc2; j++) {
	PROTECT(uij = VECTOR_ELT(U, j));
	uijr = REAL(uij);
	PROTECT(vectmp = allocVector(REALSXP, nl1));
	vtmp = REAL(vectmp);
	
	for (i = 0; i < nl1; i++) {
	    vtmp[i] = 0.0;
	    for (k = 0; k < nc1; k++) {
		vtmp[i] += REAL(VECTOR_ELT(X, k))[i] * uijr[k];
	    }
	}
	SET_VECTOR_ELT(resu, j, vectmp);
	UNPROTECT(2);
    }
    UNPROTECT(1);
    return(resu);
}





/* This function performs a PCA of the data.frame df. 
   IF sphere = 1, this function returns the normed row scores (component l1 in dudi.pca).
   If sphere = 0, the function returns the non-normed scores (component li)
*/
SEXP pca(SEXP df, SEXP sphere)
{
    SEXP mat, ov, diago, resu, spherev, listesortie;
    double *matr, *vectmp, *valpro;
    int i, j, spherei, nvp, nl;
    
    PROTECT(ov = allocVector(INTSXP, 1));
    PROTECT(listesortie = allocVector(VECSXP, 2));
    INTEGER(ov)[0] = 0;
    PROTECT(spherev = coerceVector(sphere, INTSXP));
    spherei = INTEGER(spherev)[0];
    UNPROTECT(1);
    nvp = 0;
    nl = 0;
    
    /* Inertia matrix */
    PROTECT(mat = matcovar(df));
    
    /* eigendecomposition */
    PROTECT(diago = eigendecompo(mat));
    
    /* projection of the rows */
    PROTECT(resu = XU(df, VECTOR_ELT(diago,1)));
	
    /* sphering, if asked */
    if (spherei) {
	nvp = length(VECTOR_ELT(diago, 0));
	nl = length(VECTOR_ELT(df, 0));
	
	valpro = REAL(VECTOR_ELT(diago, 0));
	for (i = 0; i < nvp; i++) {
	    vectmp = REAL(VECTOR_ELT(resu, i));
	    for (j = 0; j < nl; j++) {
		vectmp[j] = vectmp[j]/sqrt(valpro[i]);
	    }
	}
    }
        
    SET_VECTOR_ELT(listesortie, 0, resu);
    SET_VECTOR_ELT(listesortie, 1, diago);
    UNPROTECT(5);

    return(listesortie);
}






/* This function is used by the function ProjectionPursuit to sample nrepr random
   directions in the n-dimensional space defined by the columns of the data.frame df,
   and to measure the value of the criterion programmed in the R function foo in 
   these directions.
   par and envir are respectively the list of arguments needed by the R function foo
   and the environment in which this function will be evaluated
 */
SEXP DistriCritere(SEXP df, SEXP nrepr, SEXP foo, SEXP par, SEXP envir)
{
    SEXP nrepv, axe, resu;
    int nrep, r, i, ok;
    double *axer, su, *resur;
    
    PROTECT(nrepv=coerceVector(nrepr, INTSXP));
    nrep = INTEGER(nrepv)[0];
    UNPROTECT(1);
    PROTECT(axe=allocVector(REALSXP, length(df)));
    PROTECT(resu=allocVector(REALSXP, nrep));
    axer = REAL(axe);
    resur = REAL(resu);
    
    for (r=0; r < nrep; r++) {
	ok = 0;
	while (!ok) {
	    /* Generates a random direction */
	    for (i = 0; i < length(df); i++) {
		GetRNGstate();
		axer[i] = (unif_rand()*2.0) - 1.0;
		PutRNGstate();
	    }
	    /* checks that the length does not exceed 1 */
	    su = 0.0;
	    for (i = 0; i < length(df); i++) {
		su += R_pow(axer[i], 2.0);
	    }
	    if (su < 1.0)
		ok = 1;
	}
	/* restandardisation */
	for (i = 0; i < length(df); i++) {
	    axer[i] = axer[i]/sqrt(su);
	}
	/* evaluation of the function */
	resur[r] = evalfoncSuraxe(axe, df, foo, par, envir);
    }
    UNPROTECT(2);
    return(resu);
}





/* The projection pursuit approach as such. This function takes the following arguments:
   - df is the data.frame to be explored
   - foo is the R function to be maximised
   - par is the list of arguments to be passed to the function foo
   - envir is the environment where this function will be evaluated
   - nrepr is the number of repetitions passed to the function DistriCritere
     (number of random directions in the n-dimensional space where the function will
      be evaluated)
   - maxitr is the maximum number of iterations of the Nelder-Mead algorithm
   - naxesr is the number of orthogonal axes required.
   - if initr=0 the starting value is calculated with initpca. If
     initr>0 (default in R), the starting direction is randomly sampled.

     output: a SEXP corresponding to a list with one element per returned axis (i.e. with
     naxesr elements), each element being a list itself with the following elements:
     - The coordinates of the optimal axis
     - The results of the pca performed on the "current" table to eliminate superfluous
       dimensions (see below)
     - The results of the simulations carried out with the function DistriCritere
     - The observed value of the criterion on the optimal axis
 */
SEXP ProjectionPursuit(SEXP df, SEXP foo, SEXP par, SEXP envir, 
		       SEXP nrepr, SEXP maxitr, SEXP naxesr, 
		       SEXP initr)
{
    int nrep, init,  i, naxes,j;
    SEXP nrepv, dfc, dfc2, initv, axe, axes, noax, axecur, simus, liso, verbo, respca, obs;
    SEXP sphere;

    /* Check that nrep, naxesr are actually integers (maxitr are
       checked in other functions -- initr is checked below)
     */
    PROTECT(verbo = allocVector(INTSXP, 1));
    INTEGER(verbo)[0] = 1;
    PROTECT(nrepv = coerceVector(nrepr, INTSXP));
    PROTECT(noax = coerceVector(naxesr, INTSXP));
    PROTECT(initv = coerceVector(initr, INTSXP));
    init = INTEGER(initv)[0];
    nrep = INTEGER(nrepv)[0];
    naxes = INTEGER(noax)[0];
    UNPROTECT(4);    

    /* prepare the output */
    PROTECT(liso = allocVector(VECSXP, naxes));
    PROTECT(sphere = allocVector(INTSXP, 1));
    INTEGER(sphere)[0] = 0;
    dfc = df;
    dfc2 = df;

    /* A loop with as many repetitions as desired axes. For each step:
       1. We get the row scores of the PCA of the table (because the table with 
          n columns may store coordinates in a (n-1)-dimensional space -- see #5:
	  if this is step >1, the table has been projected orthogonally to the 
	  optimal vector found at the previous step -- the PCA then returns a table
	  with n-1 columns).
       2. Initialization of the search (pca, rand, initr)
       3. Search in itself
       4. Use of the function DistriCritere to sample random directions and calculate the
          criterion (and comparison with the value corresponding to the optimal axis)
       5. Projection of the table orthogonally to this optimal axis
     */    

    for (i = 0; i < naxes; i++) {
	/* Preliminary PCA to eliminate superfluous dimensions */
	PROTECT(respca=pca(dfc2, sphere));
	PROTECT(obs = allocVector(REALSXP, 1));
	REAL(obs)[0]=0.0;
	dfc = VECTOR_ELT(respca, 0);
	
	/* initialization */
	if (init == 0) {
	    PROTECT(axe = initpca(dfc, foo, par, envir));
	} 
	if (init == 1) {
	    PROTECT(axe = initrand(dfc));
	}
		
	/* Search of the best axis */
	PROTECT(axecur = trouvaxe(dfc, axe, foo, par, envir, maxitr, verbo));
	PROTECT(axes = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(axes, 0, axecur);
	SET_VECTOR_ELT(axes, 1, respca);

	/* Districritere */
	PROTECT(simus=DistriCritere(dfc, nrepr, foo, par, envir));
	SET_VECTOR_ELT(axes, 2, simus);
	
	/* Evaluation on the optimal direction */
	REAL(obs)[0]=evalfoncSuraxe(axecur, dfc, foo, par, envir);
	SET_VECTOR_ELT(axes, 3, obs);

	/* Store the results and clean the memory */
	SET_VECTOR_ELT(liso, i, axes);
	if (i==0) {
	    UNPROTECT(6);
	} else {
	    UNPROTECT(7);
	}

	/* Projection orthogonally to this axis */
	if (i<(naxes-1)) {
	    PROTECT(dfc2=orthoProjDf(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(liso,i), 1),0), 
				     VECTOR_ELT(VECTOR_ELT(liso,i),0)));
	}
    }
    
    UNPROTECT(2);
    return(liso);
}












/* Same function as ProjectionPursuit, but the axis is searched using essaisr
   random starting points and the best axis is returned
*/
SEXP ProjectionPursuitViolent(SEXP df, SEXP foo, SEXP par, SEXP envir, 
			      SEXP nrepr, SEXP maxitr, SEXP naxesr, 
			      SEXP essaisr)
{
    int nrep, init, i, naxes,j, nessais, k;
    SEXP nrepv, dfc, dfc2, initv, axe, axenow, axes, noax, axecur;
    SEXP simus, liso, verbo, respca, obs, essaisi, sphere;
    double val, valnow;

    
    /* Check that nrep, naxesr, essaisi are actually integers (maxitr is
       checked in other functions)
     */
    PROTECT(verbo = allocVector(INTSXP, 1));
    INTEGER(verbo)[0] = 1;
    PROTECT(nrepv = coerceVector(nrepr, INTSXP));
    PROTECT(noax = coerceVector(naxesr, INTSXP));
    PROTECT(essaisi = coerceVector(essaisr, INTSXP));
    nrep = INTEGER(nrepv)[0];
    naxes = INTEGER(noax)[0];
    nessais = INTEGER(essaisi)[0];
    UNPROTECT(4);
    

    /* prepare the output */
    PROTECT(liso = allocVector(VECSXP, naxes));
    PROTECT(sphere = allocVector(INTSXP, 1));
    INTEGER(sphere)[0] = 0;
    dfc = df;
    dfc2 = df;

    
    for (i = 0; i < naxes; i++) {
	/* Preliminary PCA to eliminate superfluous dimensions */
	PROTECT(respca=pca(dfc2, sphere));
	PROTECT(obs = allocVector(REALSXP, 1));
	REAL(obs)[0]=0.0;
	dfc = VECTOR_ELT(respca, 0);
	
	/* initialization */
	PROTECT(axe = initrand(dfc));
	PROTECT(axenow = initrand(dfc));
	val=evalfoncSuraxe(axe, dfc, foo, par, envir);
	valnow=evalfoncSuraxe(axenow, dfc, foo, par, envir);
	if (valnow > val) {
	    for (j = 0; j < length(axe); j++) {
		REAL(axe)[j] = REAL(axenow)[j];
		val = valnow;
	    }
	}

	/* axe contains the best axis -- axenow contains the current axis: */
	for (j = 0; j < nessais; j++) {
	    PROTECT(axenow = initrand(dfc));
	    PROTECT(axecur = trouvaxe(dfc, axenow, foo, par, envir, maxitr, verbo));
	    valnow = evalfoncSuraxe(axecur, dfc, foo, par, envir);
	    if (valnow > val) {
		for (k = 0; k < length(axe); k++) {
		    REAL(axe)[k] = REAL(axecur)[k];
		    val = valnow;
		}
 	    }
	    UNPROTECT(2);
	}
	
	
	PROTECT(axes = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(axes, 0, axe);
	SET_VECTOR_ELT(axes, 1, respca);

	/* Districritere */
	PROTECT(simus=DistriCritere(dfc, nrepr, foo, par, envir));
	SET_VECTOR_ELT(axes, 2, simus);
	
	/* Evaluation on the optimal direction */
	REAL(obs)[0]=evalfoncSuraxe(axe, dfc, foo, par, envir);
	SET_VECTOR_ELT(axes, 3, obs);

	/* Store the results and clean the memory */
	SET_VECTOR_ELT(liso, i, axes);
	if (i==0) {
	    UNPROTECT(6);
	} else {
	    UNPROTECT(7);
	}

	/* Projection orthogonally to this axis */
	if (i<(naxes-1)) {
	    PROTECT(dfc2=orthoProjDf(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(liso,i), 1),0), 
				     VECTOR_ELT(VECTOR_ELT(liso,i),0)));
	}
    }
    
    UNPROTECT(2);
    return(liso);
}





/* Standardize a vector to unit length. vec is the pointer to the table, and n is the length
   of the table
   
   no output: the vector vec is standardized
 */
void StandardizeVec(double *vec, int n)
{
    double su;
    int i;
    
    su=0.0;
    for (i = 0; i < n; i++) {
	su += R_pow(vec[i], 2.0);
    }
    su = sqrt(su);
    for (i = 0; i < n; i++) {
	vec[i] = vec[i]/su;
    }
    return;
}


/* Sample a random plane in a C-dimensional space, where C is the number of elemnts/columns
   in the list/data.frame df.
   
   output: a list with two elements, each element being one of the two vectors defining 
   the plane
 */
SEXP initplan(SEXP df)
{
    SEXP sorties, axe1, axe2, axe1b;
    double costheta, su;
    int i;
    
    /* generates two random axes */
    PROTECT(axe1 = initrand(df));
    PROTECT(axe2 = initrand(df));
    PROTECT(axe1b = allocVector(REALSXP, length(axe2)));
    PROTECT(sorties = allocVector(VECSXP, 2));
    
    costheta = 0;
    for (i = 0; i < length(df); i++) {
	costheta += REAL(axe1)[i]*REAL(axe2)[i];
    }
    for (i = 0; i < length(df); i++) {
	REAL(axe1b)[i] = REAL(axe2)[i]-REAL(axe1)[i]*costheta;
    }

    StandardizeVec(REAL(axe1b), length(axe1b));
    SET_VECTOR_ELT(sorties, 0, axe1);
    SET_VECTOR_ELT(sorties, 1, axe1b);
    UNPROTECT(4);
    return(sorties);    
}




/* Évaluates the function fonc with the parameters stored in the list par in the environment
   envir, based on the vector storing the coordinates of the rows of the data.frame df
   on the plane defined by the two directions defined by the orthogonal vectors u1 and u2, 
   where u1 and u2 are respectively the first and second elements of the unnamed list plan.
   
   output: a double containing the value of the function
*/
double evalfoncSurplan(SEXP plan, SEXP df, SEXP fonc, SEXP par, SEXP envir)
{
    int n, nl, i, j;
    double *xr1, *xr2, *dfi, *vecr1, *vecr2, out;
    SEXP vec1, vec2, resu, xSAT1, xSAT2, dfproj;

    if(!isEnvironment(envir)) error("'env' should be an environment");
    
    /* Préparation of the structures to store the results */ 
    xSAT1 = VECTOR_ELT(plan,0);
    xSAT2 = VECTOR_ELT(plan,1);
    n = length(xSAT1);
    nl = length(VECTOR_ELT(df, 0));

    PROTECT(vec1 = allocVector(REALSXP, nl));
    PROTECT(vec2 = allocVector(REALSXP, nl));
    PROTECT(dfproj = allocVector(VECSXP, 2));

    xr1 = REAL(xSAT1);   
    xr2 = REAL(xSAT2);   
    
    /* projection on df on x */
    vecr1 = REAL(vec1);
    vecr2 = REAL(vec2);
    for (i = 0; i < nl; i++) {
	vecr1[i] = 0.0;
	vecr2[i] = 0.0;
	for (j = 0; j < n; j++) {
	    vecr1[i] += xr1[j] * REAL(VECTOR_ELT(df, j))[i];
	    vecr2[i] += xr2[j] * REAL(VECTOR_ELT(df, j))[i];
	}
    }
    
    SET_VECTOR_ELT(dfproj,0,vec1);
    SET_VECTOR_ELT(dfproj,1,vec2);

    /* evaluation of the function on vec */
    defineVar(install("x"), dfproj, envir);
    defineVar(install("par"), par, envir);
    PROTECT(resu = eval(fonc, envir));
    out = REAL(resu)[0];    
    UNPROTECT(4);

    return(out);
}


/* Just a wrapper for the previous function, to test the function in R:
   same argument as evalfoncSurplan, but the output is now a SEXP
 */
SEXP wrapevfsp(SEXP plan, SEXP df, SEXP fonc, SEXP par, SEXP envir)
{
    SEXP so;
    PROTECT(so=allocVector(REALSXP,1));
    REAL(so)[0]=evalfoncSurplan(plan, df, fonc, par, envir);
    UNPROTECT(1);
    return(so);
}


/* This function implements the algorithm of Posse (1995). The arguments are:
   df: the data.frame containing the data
   foo: the R function implementing the criterion
   par: the R object containing the data needed to calculate the criterion (may be NULL)
   envir: the R environment where the criterion will be calculated
   conr: the value of c (see the help page of posse1995)
   maxitr: the value of maxit (see the help page of posse1995)
   tolr: the value of tol (see the help page of posse1995)
   nrep: the value of nrep (see the help page of posse1995)
   maxhalf: the value of maxhalf (see the help page of posse1995)

   output: a SEXP corresponding to a list with two elements:
   - A list with the two vectors defining the plane found by the algorithm
   - The results of the simulations carried out to measure the importance of the
     criterion on the found plane. This element is a list with two elements:
     - The observed value of the criterion on the optimal axis
     - The nrep simulated values of the criterion on randomly sampled planes
 */

SEXP trouveplan(SEXP df, SEXP foo, SEXP par, SEXP envir,
		SEXP conr, SEXP maxitr, SEXP tolr, SEXP nrepr, SEXP maxhalf)
{

    SEXP v, a1, a2, b1, b2, alpha, beta, plane2;
    SEXP cond, maxiti, told, nrepi, maxhalfd, plane, alphatmp, betatmp, simus, obs, resim, liso;
    int i, j, cont, maxit, m, ite, nrep;
    double c, tol, a1b, a2b, fp, fm, f, g, half;

    /* Checks that the format is correct for the arguments */
    PROTECT(cond = coerceVector(conr, REALSXP));
    PROTECT(told = coerceVector(tolr, REALSXP));
    PROTECT(maxhalfd = coerceVector(maxhalf, INTSXP));
    PROTECT(maxiti = coerceVector(maxitr, INTSXP));
    PROTECT(nrepi = coerceVector(nrepr, INTSXP));
    c=REAL(cond)[0];
    tol=REAL(told)[0];
    m=INTEGER(maxhalfd)[0];
    maxit=INTEGER(maxiti)[0];
    nrep = INTEGER(nrepi)[0];
    UNPROTECT(5);
    fp=0;
    fm=0;
    f=0;
    g=0;
    half=0;
    
    /* Starts with a random plane */
    PROTECT(plane=initplan(df));
    
    /* allocates memory for the working vectors */
    PROTECT(alpha = allocVector(REALSXP, length(df)));
    PROTECT(beta = allocVector(REALSXP, length(df)));
    PROTECT(a1 = allocVector(REALSXP, length(df)));
    PROTECT(a2 = allocVector(REALSXP, length(df)));
    PROTECT(b1 = allocVector(REALSXP, length(df)));
    PROTECT(b2 = allocVector(REALSXP, length(df)));
    PROTECT(alphatmp = allocVector(REALSXP, length(df)));
    PROTECT(betatmp = allocVector(REALSXP, length(df)));
    PROTECT(simus = allocVector(REALSXP, nrep));
    PROTECT(obs = allocVector(REALSXP, 1));
    PROTECT(resim = allocVector(VECSXP, 2));
    PROTECT(liso = allocVector(VECSXP, 2));
    


    /* Copy the current plane in alpha and beta */
    for (i = 0; i < length(df); i++) {
	REAL(alpha)[i]=REAL(VECTOR_ELT(plane,0))[i];
	REAL(beta)[i]=REAL(VECTOR_ELT(plane,1))[i];
    }
    
    /* starts the loop */
    cont=1;
    ite=0;
    while (cont) {
 	/* Random vector v */
	PROTECT(v = initrand(df));
	
	/* a1 and a2 */
	for (i = 0; i < length(df); i++) {
	    REAL(a1)[i] = REAL(alpha)[i] + c*REAL(v)[i];
	    REAL(a2)[i] = REAL(alpha)[i] - c*REAL(v)[i];
	}
	StandardizeVec(REAL(a1), length(a1));
	StandardizeVec(REAL(a2), length(a1));
	
	/* a1^t beta et a2^t beta */
	a1b = 0.0;
	a2b = 0.0;
	for (i = 0; i < length(df); i++) {
	    a1b += REAL(a1)[i]*REAL(beta)[i];
	    a2b += REAL(a2)[i]*REAL(beta)[i];
	}
	
	/* b1 and b2 */
	for (i = 0; i < length(df); i++) {
	    REAL(b1)[i] = REAL(beta)[i] - a1b * REAL(a1)[i];
	    REAL(b2)[i] = REAL(beta)[i] - a2b * REAL(a2)[i];
	}
	StandardizeVec(REAL(b1), length(b1));
	StandardizeVec(REAL(b2), length(b1));
	
	/* plane + */
	SET_VECTOR_ELT(plane,0,a1);
	SET_VECTOR_ELT(plane,1,b1);
	fp=evalfoncSurplan(plane, df, foo, par, envir);
	SET_VECTOR_ELT(plane,0,a2);
	SET_VECTOR_ELT(plane,1,b2);
	fm=evalfoncSurplan(plane, df, foo, par, envir);

	if (fp > fm) {
	    f = fp;
	    for (i = 0; i < length(df); i++) {
		REAL(alphatmp)[i] = REAL(a1)[i];
		REAL(betatmp)[i] = REAL(b1)[i];
	    }
	} else {
	    f = fm;
	    for (i = 0; i < length(df); i++) {
		REAL(alphatmp)[i] = REAL(a2)[i];
		REAL(betatmp)[i] = REAL(b2)[i];
	    }
	}
	SET_VECTOR_ELT(plane,0,alpha);
	SET_VECTOR_ELT(plane,1,beta);
	g=evalfoncSurplan(plane, df, foo, par, envir);
	
	if (f > g) {
	    for (i = 0; i < length(df); i++) {
		REAL(alpha)[i] = REAL(alphatmp)[i];
		REAL(beta)[i] = REAL(betatmp)[i];
	    }
	    half = 0;
	} else {
	    half = half+1;
	}
	if (half>m) {
	    c = c/2.0;
	    half = 0;
	}

	if (c <tol)
	    cont = 0;
	ite++;
	if (ite > maxit)
	    error("Maximum number of iterations reached");
	UNPROTECT(1);
    }
    SET_VECTOR_ELT(plane,0,alpha);
    SET_VECTOR_ELT(plane,1,beta);
    
    /* randtest */
    REAL(obs)[0] = g;
    
    /* measure of the criteria on random planes */
    for (i = 0; i < nrep; i++) {
	PROTECT(plane2 = initplan(df));
	REAL(simus)[i]=evalfoncSurplan(plane2, df, foo, par, envir);
	UNPROTECT(1);
    }

    SET_VECTOR_ELT(resim,0,obs);
    SET_VECTOR_ELT(resim,1,simus);
    SET_VECTOR_ELT(liso,0,plane);
    SET_VECTOR_ELT(liso,1,resim);
    

    UNPROTECT(13);
    return(liso);
}



/* **************************************************
   The functions below implement the MaxEnt approach
   ************************************************** */



/* 
   Original implementation of maxent, as described in the paper of Phillips et al. in 
   Ecological Modelling. I do not use here the exact same algorithm as in their paper, 
   but a more classical BFGS algorithm (quasi-Newton approach).

   The result returned by this approach is identical as the Poisson regression 
   described by Renner et Warton. Note that this result is identical to the one that would 
   be obtained with  the function glm, and this, despite the warnings caused by the fact 
   that the response variable is not a Poisson variable
 */



/* structure required for the optimization, for the use of the function vmmin (implementing
   the BfGS algorithm: this structure is intended to store the data required for the 
   calculation of the penalized likelihood (the penalty coefficients beta, 
   the response y, and the features used for prediction (data.frame) */
typedef struct {
    double *bbeta;
    SEXP featuresp;
    double *y;
} Pouroptim;



/* The following function calculates the core of the likelihood for MaxEnt. Let us recall
   that the distribution that we want to fit is a Gibbs distribution:
   
   exp(lambda%*%x)/Z_lambda
   
   where Z_lambda is a constante, depending on the parameters in lambda, allowing this
   distribution to integrate to 1 on the whole domain (which is a problem, since this domain
   is continous, we are constrained to discretize it by sampling context points in it, and 
   calculating Z_lambda by "summing" the function exp(lambda%*%x) over all these context 
   points.
   
   Actually, the function below calculates:
   (i) either -1/m *(\sum_i lambda %*% x_i), where i loops over all presence points
   (ii) or Z_lambda, calculated as \sum_j exp(lambda %*% x_j), where j loops over all points
   (context and presence pooled)

   Parameters:
   * featuresp: SEXP data.frame containing the features
   * N: Number of points (number of rows in SEXP)
   * m: Number of presence points
   * P: number of features
   * y: pointer to a table of double describing the type of points in featuresp
   (y=1 is a presence, y=0 is a context point)
   * lambdar: pointer to a table containing the coefficients of the model
   * lg: if 1, then we calculate (i) above, otherwise we calculate (ii)
   
*/

double calcCoeurVrais(SEXP featuresp, int N, int m, int P, double *y, double *lambdar, int lg)
{
    double vrais, cl;
    int i, j;
    
    vrais = 0.0;
    for (i = 0; i < N; i++) {
	if (y[i] >0.5) {
	    cl = 0.0;
	    for (j=0; j<P; j++) {
		cl += REAL(VECTOR_ELT(featuresp, j))[i]*lambdar[j];
	    }
	    if (lg ==1) {
		vrais += cl;
	    } else {
		vrais += exp(cl);
	    }
	}
    }
    if (lg == 1) {
	vrais = -(vrais/(((double) m)));
    } else {
	vrais = log(vrais);
    }
    return(vrais);
}


/* This function calculates the penalized likelihood:
 * bbeta: value of the coefficient used for the LASSO regularization
 * lambdar: pointer towarrd the table containing the coefficients
 * featuresp: SEXP containing the features
 * y: response variable (cf. previous function).
 */

double calcVraisPen(double *bbeta, double *lambdar, SEXP featuresp, double *y)
{
    int m, i, j, P, N;
    double vrais, cl, vrais2;
    SEXP ybisr;

    m=0;
    N=length(VECTOR_ELT(featuresp,0));
    for (i = 0; i < N; i++) {
	if (y[i] > 0.5)
	    m++;
    }
    P=length(featuresp);
    
    /* Calculation of the likelihood */
    vrais = calcCoeurVrais(featuresp, N, m, P, y, lambdar, 1);

    /* renormalisation (calculation of the renormalisation constant on the context points) */
    PROTECT(ybisr = allocVector(REALSXP, N));
    for (i = 0; i < N; i++) {
	REAL(ybisr)[i] = 1;
    }
    
    vrais2 = calcCoeurVrais(featuresp, N, m, P, REAL(ybisr), lambdar, 0);
    vrais=vrais + vrais2;
    
    /* Regularisation l1 */
    for (j=0; j<P; j++) {
	vrais += fabs(lambdar[j])*bbeta[j];
    }
    
    /* output */
    UNPROTECT(1);
    return(vrais);
}


/* Numerical calculation of the gradient vector. 
   Same arguments as the previous function
*/
SEXP calcGradient(double *bbeta, double *lambdar, SEXP featuresp, double *y)
{
    double vr1, vr2, h;
    SEXP gradient;
    int i;
    
    PROTECT(gradient = allocVector(REALSXP, length(featuresp)));
    h = 1e-5;

    for (i = 0; i < length(featuresp); i++) {
	lambdar[i] = lambdar[i]+h;
	vr1 = calcVraisPen(bbeta, lambdar, featuresp, y);
	lambdar[i] = lambdar[i]-2*h;
	vr2 = calcVraisPen(bbeta, lambdar, featuresp, y);
	REAL(gradient)[i] = (vr1-vr2)/(2.0*h);
	lambdar[i] = lambdar[i]+h;
    }
    UNPROTECT(1);
    return(gradient);
}



/* Function having the correct format for use with vmmin -- calculates the 
   penalized likelihood
*/
double vraispenpropt(int n, double *par, void *ex)
{
    double res;
    Pouroptim *ex2 = (Pouroptim *) ex;
        
    res = calcVraisPen(ex2->bbeta, par, ex2->featuresp, ex2->y);
    
    return(res);
}


/* Function calculating the gradient vector, having the correct format for use with
   vmmin
*/
void gradpropt(int n, double *par, double *gr, void *ex)
{
    SEXP res;
    int i;
    Pouroptim *ex2 = (Pouroptim *) ex;

    res = calcGradient(ex2->bbeta, par, ex2->featuresp, ex2->y);
    
    for (i = 0; i < length(ex2->featuresp); i++) {
	gr[i] = REAL(res)[i];
    }
}



/* Wrapper for the function vmmin, for the BFGS optimization:
   * lambdaInit: SEXP containing the initial values of lambda
   * bbeta: vector containing the value of the beta coefficients (LASSO penalization)
   * featuresp: SEXP containing the features
   * maxIter: maximum number of iterations of the BFGS algorithm
   * stopCrit: criterion to detect convergence
   * verb: verbose or not?
   * y: vector containing the "response variable" (1=presence/0=context points).
*/
SEXP EstimMaxVraisPen(SEXP lambdaInit, double *bbeta, SEXP featuresp,
		      int maxIter, double stopCrit, int verb, double *y)
{
    
    Pouroptim ex;
    double *lambdar, vrais1, vraisPenInit;
    int *mask, i, nREPORT;
    int fncount, grcount, convergence;
    SEXP lambda, resultats, nbits, convergences;
    
    /* Required by the algorithm bgfs */
    double abstol = log(0.0);
    double reltol = stopCrit;
    nREPORT = 1;
    fncount = 0;
    grcount = 0;
    convergence = 0;
    ex.bbeta = bbeta;
    ex.featuresp = featuresp;
    ex.y = y;

    /* The vector of parameters */
    PROTECT(resultats = allocVector(VECSXP, 3));
    PROTECT(nbits = allocVector(INTSXP, 1));
    PROTECT(convergences = allocVector(INTSXP, 1));
    PROTECT(lambda = allocVector(REALSXP, length(featuresp)));
    lambdar = REAL(lambda);
    for (i = 0; i < length(featuresp); i++)
	lambdar[i] = REAL(lambdaInit)[i];
    
    /* Starting penalized likelihood */
    vrais1 = calcVraisPen(bbeta, lambdar, featuresp, y);
    vraisPenInit = vrais1;
    
    mask = (int *) R_alloc(length(lambdaInit), sizeof(int));
    for (i = 0; i < length(lambdaInit); i++) 
	mask[i] = 1;

    vmmin(length(lambdaInit), lambdar, &vrais1,
	  vraispenpropt, gradpropt, maxIter, verb,
	  mask, abstol, reltol, nREPORT,
	  &ex, &fncount, &grcount, &convergence);
    /* The arguments of the function vmmin are the following:
       - number of parameters
       - Starting values of the parameters
       - Value of the likelihood that should be returned by the function
       - function for the calculation of the penalized likelihood
       - function for the calculation of the gradient
       - maximum number of iterations
       - verbose?
       - a vector indicating which elements of the vector of parameters should 
       be changed in this minimisation (mask=1), and which elements should not be taken
       into account (mask=0). Here, all the elements are equal to 1
       - absolute tolerance (the function should not be lower than this value)
       - relative tolerance (used to identify the convergence)
       - show a report every how many iterations?
       - data used by the likelihood and gradient functions (as object ex)
       - as output: the number of times that optim calls the penalized likelihood function
       - as output: the number of times that optim calls the gradient function
       - as output: whether the optimization converged
    */
    
    if (convergence!=0)
	Rprintf("Minimization failed. \nCode returned by the function optim (Bfgs algo): %i (cf. ?optim)\n", convergence);
    INTEGER(nbits)[0] = fncount;
    INTEGER(convergences)[0] = convergence;
    SET_VECTOR_ELT(resultats, 0, lambda);
    SET_VECTOR_ELT(resultats, 1, nbits);
    SET_VECTOR_ELT(resultats, 2, convergences);
    
    UNPROTECT(4);
    return(resultats);
}



/* Interface to R. Wrapper for the function EstimMaxVraisPen. see previously for the
   arguments */
SEXP maxentBFGS(SEXP featuresp, SEXP bbeta, SEXP mi, SEXP sc, SEXP ve, SEXP lambda, SEXP y)
{
    SEXP lambdafin;
    int i;

    PROTECT(lambdafin = EstimMaxVraisPen(lambda, REAL(bbeta), featuresp,
					 INTEGER(mi)[0], REAL(sc)[0], INTEGER(ve)[0], REAL(y)));
    UNPROTECT(1);
    return(lambdafin);
    
}




/* 
   Implementation of maxent based on Dudik and Phillips 2004. We rely here on a 
   sequential updating algorithm: the values of lambda are changed individually
   when they allow to reduce the penalized likelihood.
   
   We obtain the same result as the BFGS algorithm, only if we do not have any intercept:
   this algorithm requires features scaled between 0 and 1, which is not possible 
   with an intercept (see Dudik and Phillips, op.cit).
 */


/* Calculation of Z_lambda (coefficient of standardization of the Gibbs distribution)
   for a given value of lambda, with featuresp the feature dataframe */

double Z_lambda(SEXP featuresp, double *lambda)
{
    double Z,cl;
    int i,j;
    Z=0.0;
    for (i = 0; i < length(VECTOR_ELT(featuresp,0)); i++) {
	cl = 0.0;
	for (j = 0; j < length(featuresp); j++) {
	    cl+= (REAL(VECTOR_ELT(featuresp,j))[i] * lambda[j]);
	}
	Z+=exp(cl);	
    }
    return(Z);    
}

/* Expectation of a feature k in the dataframe of features featuresp, calculated with the
   fitted Gibbs distribution (Z is the coefficient of standardization of the distribution) */

double E_q_lambda_fk(SEXP featuresp, double *lambda, int k, double Z)
{

    double eql, cl;
    int i,j;

    eql=0.0;    
    for (i = 0; i < length(VECTOR_ELT(featuresp,0)); i++) {
	cl = 0.0;
	for (j = 0; j < length(featuresp); j++) {
	    cl+=REAL(VECTOR_ELT(featuresp,j))[i] * lambda[j];
	}
	eql+=(exp(cl)/Z)*REAL(VECTOR_ELT(featuresp,k))[i];	
    }
    return(eql);
}


/* Mean of a feature k calculated from the sample of observations
   (featurespy is a dataframe containing the value of the features for each observation
   of the species). */

double pi_tilde_fk(SEXP featurespy, int k)
{
    double ep;
    int i, j;
    
    ep=0.0;    
    for (i = 0; i < length(VECTOR_ELT(featurespy,0)); i++) {
	ep+=REAL(VECTOR_ELT(featurespy,k))[i];	
    }
    ep=ep/((double) length(VECTOR_ELT(featurespy,0)));
    return(ep);
}

/* Calculation of the possible value of delta when lambda + delta is >0
   The arguments are the mean of the feature measured on the sample (ep),
   the vector of the penalization coefficients betaj, the expectation of the
   fitted Gibbs distribution (eql), and the number of the feature considered
*/
double delta_lambdapos(double ep, double *betaj, double eql, int j)
{
    double res;
    res = log(((ep - betaj[j])*(1.0 - eql))/((1.0-ep+betaj[j])*eql));
    return(res);
}


/* Calculation of the possible value of delta when lambda + delta is < 0
   The arguments are the mean of the feature measured on the sample (ep),
   the vector of the penalization coefficients betaj, the expectation of the
   fitted Gibbs distribution (eql), and the number of the feature considered
*/
double delta_lambdaneg(double ep, double *betaj, double eql, int j)
{
    double res;
    res = log(((ep + betaj[j])*(1.0 - eql))/((1.0-ep-betaj[j])*eql));
    return(res);
}

/* The value of Fj that we want to minimize at each step of the algorithm;
   equation (12) of the paper of Dudik et Phillips.
   We give the value of delta, the mean of the feature, and its expectation
   from the Gibbs distribution fitted at time t, the penalisation betaj corresponding
   to this feature, and the coefficient lambdaj corresponding to this feature.
*/
double Fjf(double delta, double ep, double eq, double betaj, double lambdaj)
{
    double res;

    res = -delta*ep + log(1.0 + (exp(delta) - 1.0)*eq) + 
	betaj*(fabs(lambdaj+delta) - fabs(lambdaj));
    return(res);
}

/* Calculation of the likelihood of the sample
   We give the data.frame featuresp containing the features, the number m of species 
   observations, the vector y indicating (0/1) the context points/presence, The vector 
   of coefficients lambdar, and the standardization coefficient of the Gibbs distribution 
   associated to the value of the coefficients lambdar
 */
double vraisemblance(SEXP featuresp, int m, double *y, double *lambdar, double Zl)
{
    double vrais, cl;
    int i, j;
    
    vrais = 0.0;
    for (i = 0; i < length(VECTOR_ELT(featuresp,0)); i++) {
	if (y[i] >0.5) {
	    cl = 0.0;
	    for (j=0; j< length(featuresp); j++) {
		cl += REAL(VECTOR_ELT(featuresp, j))[i]*lambdar[j];
	    }
	    vrais += cl;
	}
    }
    vrais = -(vrais/(((double) m))) + log(Zl);
    return(vrais);
}


/* This function calculates the penalized likelihood:
 * betaj: value of the coefficient serving for the LASSO penalization
 * lambdar: pointer to the table of coefficients
 * featuresp: SEXP containing the features
 * y: response variable (as before).
 * m: number of presences of the species in the dataset
 * Zl: standardisation coefficient of the Gibbs distribution for the value of lambdar
 */

double vraisemblance_penalisee(double *betaj, double *lambdar, SEXP featuresp, 
			       double *y, int m, double Zl)
{
    int i, j, N;
    double vrais, cl, vrais2;
    
    /* Calculation of the likelihood */
    vrais = vraisemblance(featuresp, m, y, lambdar, Zl);
    
    /* Regularisation l1 */
    for (j=0; j<length(featuresp); j++) {
	vrais += fabs(lambdar[j])*betaj[j];
    }
    
    /* output */
    return(vrais);
}



/* Main algorithm of maxent. It corresponds exactly to the implementation of figure 1
   in Dudik and Phillips (2004). This procedure changes the values in *lambda as output.

   featuresp is the dataframe containing the features on all points
   featurespy is the dataframe containing the features on the presence points
   betaj is the vector containing the coefficients for the LASSO penalization
   lambda is the vector containing the coefficients associated to the features
   (0 as input, estimations as output)
   y is the vector indicating which points in featuresp are presences (1)
   or context points (0)
   vraiscrit is the criterion used to determine convergence (based on the
   likelihood difference between two successive itérations).
   verbose = 1 indicates the changes of the penalized likelihood with time
   maxIter is selfexplanatory
*/

void algo_principal(SEXP featuresp, SEXP featurespy, double *betaj, 
		    double *lambda, double *y, double vraiscrit, int verbose, int maxIter,
		    int *nbit, int *convergencearenv)
{
    SEXP deltar,Fjr, epr, eqr;
    double *delta,  *Fj, Fjmin, Zl, *eq, *ep, lamj;
    double delta0, delta1, deltap, deltan, Fj1, Fj0, vraisprev, vraisnew;
    int j, convergence, achanger, m, i, it;
    
    /* Memory Allocation */
    PROTECT(deltar = allocVector(REALSXP, length(featuresp)));
    PROTECT(Fjr = allocVector(REALSXP, length(featuresp)));
    PROTECT(epr = allocVector(REALSXP, length(featuresp)));
    PROTECT(eqr = allocVector(REALSXP, length(featuresp)));
    delta = REAL(deltar);
    Fj = REAL(Fjr);
    ep = REAL(epr);
    eq = REAL(eqr);
    
    /* Calculation of m */
    m=0;
    for (i = 0; i < length(VECTOR_ELT(featuresp,0)); i++) {
	if (y[i] > 0.5)
	    m++;
    }

    /* Calculation of ep once and for all */
    for (j = 0; j < length(featuresp); j++) {
	ep[j] = pi_tilde_fk(featurespy, j);
    }

    /* Initialisation */
    Zl = Z_lambda(featuresp, lambda);
    vraisnew = vraisemblance_penalisee(betaj, lambda, featuresp, 
				       y, m, Zl);
    vraisprev = 0.0;
    it = 0;
    
    /* Start the loop */
    convergence = 0;
    *convergencearenv = 0;
    *nbit = 0; /* Number of iterations */
    while (!convergence) {
	R_CheckUserInterrupt();
	(*nbit)++;
	/* calculation of Fj for each variable */
	for (j = 0; j < length(featuresp); j++) {
	    
	    /* Calculation of the expectation of q estimated and of the expectation 
	       of the estimation p for the feature j */
	    eq[j] = E_q_lambda_fk(featuresp, lambda, j, Zl);

	    /* Calculation of the different possible values for delta */
	    lamj = lambda[j];
	    delta0 = -lamj;
	    deltap = delta_lambdapos(ep[j], betaj, eq[j], j);

	    /* By default, we set delta equal do -lambda_j, and we define F_j 
	       equal to the corresponding F_j */
	    delta[j] = delta0;
	    Fj0 = Fjf(delta0, ep[j], eq[j], betaj[j], lamj);
	    Fj[j] = Fj0;
	    
	    /* Identification of the minimum delta */
	    if ((lamj + deltap)>=0) {
		Fj1 = Fjf(deltap, ep[j], eq[j], betaj[j], lamj);
		delta1 = deltap;
	    } else {
		deltan = delta_lambdaneg(ep[j], betaj, eq[j], j);
		Fj1 = Fjf(deltan, ep[j], eq[j], betaj[j], lamj);
		delta1 = deltan;
	    }
	    if (Fj1 < Fj0) {
		delta[j] = delta1;
		Fj[j] = Fj1;
	    }
	}
	
	/* Identification of the minimum couple (delta,j) */
	achanger = 0;
	Fjmin = Fj[0];
	for (j = 1; j < length(featuresp); j++) {
	    if (Fj[j] < Fjmin) {
		Fjmin = Fj[j];
		achanger = j;
	    }
	}
	
	/* Update */
	lambda[achanger] = lambda[achanger]+delta[achanger];
	
	/* Calculation of Z_lambda and likelihood */
	Zl = Z_lambda(featuresp, lambda);
	vraisprev = vraisnew;
	vraisnew = vraisemblance_penalisee(betaj, lambda, featuresp, 
					   y, m, Zl);
	it++;
	if (verbose)
	    Rprintf("Iteration %i: %f\n", it, vraisnew);
	/* Tests convergence */
	if (fabs(vraisprev-vraisnew) < vraiscrit) {
	    convergence=1;
	    (*convergencearenv) = 1;
	}
	if ((*nbit)>maxIter) {
	    Rprintf("Maximum number of iterations reached.\n");
	    convergence=1;
	    *convergencearenv = 0;
	}
    }

    UNPROTECT(4);
    return;

}



/* Interface to R. Wrapper for the function algo_principal. See previously for the arguments */
SEXP maxentSequential(SEXP featuresp, SEXP featurespy, SEXP bbeta, SEXP yr, 
		      SEXP vraiscrit, SEXP verbose, SEXP mi)
{
    SEXP lambdafin, resultats, nbits, convergences;
    double *betaj, *y, *lambda;
    int j, nbit, convergence;
    
    nbit = 0;
    convergence = 0;
    betaj = REAL(bbeta);
    y = REAL(yr);
    PROTECT(lambdafin = allocVector(REALSXP, length(featuresp)));
    PROTECT(nbits = allocVector(INTSXP, 1));
    PROTECT(convergences = allocVector(INTSXP, 1));
    PROTECT(resultats = allocVector(VECSXP, 3));
    lambda = REAL(lambdafin);
    
    for (j = 0; j < length(lambdafin); j++) {
	lambda[j] = 0.0;
    }
    algo_principal(featuresp, featurespy, betaj, lambda, y, REAL(vraiscrit)[0], 
		   INTEGER(verbose)[0], INTEGER(mi)[0], &nbit, &convergence);

    INTEGER(nbits)[0] = nbit;
    INTEGER(convergences)[0] = convergence;
    SET_VECTOR_ELT(resultats, 0, lambdafin);
    SET_VECTOR_ELT(resultats, 1, nbits);
    SET_VECTOR_ELT(resultats, 2, convergences);

    UNPROTECT(4);
    return(resultats);
    
}


SEXP vraisemblance_penaliseer(SEXP betajs, SEXP lambdas, SEXP featuresp, 
			      SEXP ys)
{
    double *betaj, *lambdar, *y, Zl;
    int m, i;
    SEXP resu;

    PROTECT(resu = allocVector(REALSXP, 2));
    
    betaj = REAL(betajs);
    lambdar = REAL(lambdas);
    y = REAL(ys);
    m=0;
    for (i = 0; i < length(VECTOR_ELT(featuresp,0)); i++) {
	if (y[i] > 0.5)
	    m++;
    }
    
    Zl = Z_lambda(featuresp, lambdar);
    REAL(resu)[0] = vraisemblance_penalisee(betaj, lambdar, featuresp, 
					    y, m, Zl);
    REAL(resu)[1] = vraisemblance(featuresp, m, y, lambdar, Zl);
    UNPROTECT(1);
    
    return(resu);
}

