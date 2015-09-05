/* Compression Alg. for CRe spectra. 
 * We first find a Cubic Hermite spline that fits
 * the normalised spectrum in log-log space 
 * Then we binary compress the knots further 
 * This if full of hacks, use the existing debug infrastructure to 
 * make pathological spectra fit, good luck */

#include "common.h"
#include "compress.h"

#define N N_SPEC_BINS   // Length of spectrum 

#define SIZEOF_FULLKNOT sizeof(struct full_point) 
#define SIZEOF_HALFKNOT sizeof(struct half_point)
#define MAXKNOTS (SPECSIZE_BYTES/SIZEOF_HALFKNOT) // can't be more nodes 

#define KNOT_MIN_DIST 4 // minimum distance between two knots 

//#define DEBUG_COMPRESSION // runs only on particle IPART! 
//#define IPART 5 // ID

void check_knots(const struct Knot *K, double *spectrum);

static int iii = -1;

static float log_p[N] = { 0 }, log_dp[N] = { 0 }, dL[N] = { 0 }, dR[N] = { 0 },
             ddL[N] = { 0 }, ddR[N] = { 0 };

static size_t Nknots, KnotMemSize; // count mem size & no 

#pragma omp threadprivate(Nknots,KnotMemSize,log_p,log_dp,dL,dR,ddL,ddR,iii)

void Compress(const float *input_spectrum, float *nCRe, char *spectrum)
{   
    Assert(spectrum != NULL, "Can't compress into NULL pointer");
    Assert(input_spectrum != NULL, "Can't compress from NULL pointer");
    
	double np[N_SPEC_BINS] = { 0 };
	
	for (int i=0; i<N_SPEC_BINS; i++) 
		np[i] = pow(10, input_spectrum[i]);

    *nCRe = compute_spec_parameters(np);

    if (*nCRe == 0)
        return ;

    Nknots = KnotMemSize = 0;

    struct Knot K[MAXKNOTS] = { { 0 } };

    find_first_knots(np, K);

    double curve[N] = { -DBL_MAX }; // holds reconstructed spectrum
    int err_idx = 0;
    
	for (;;) {

        update_control_points(K);    

		draw_curve(K, curve);    
    
        double max_err = find_max_error(K, np, curve, &err_idx);

        add_knot(np, K, err_idx, KNOT_FULL);

        if (KnotMemSize + SIZEOF_FULLKNOT > SPECSIZE_BYTES)
            break; // the STOP knot is already accounted for
    }
    
	update_control_points(K);

//	check_knots(K, np);

    compress_knots_binary(K, spectrum);

#ifdef DEBUG_COMPRESSION
    
	if (Ipart == IPART) {
        char name[80];
	    FILE *fd;
        sprintf(name, "./spec_uncompressed_%03i", Snap.SnapNum);
        fd = fopen(name, "w");

		double curve2[N] = { -DBL_MAX };

        draw_curve(K, curve2);    

        for (int i = 0; i < N; i++)
            fprintf(fd, "%d %g %g %g \n", i, log_p[i], curve2[i]*nCRe[0], 
		    		input_spectrum[i]);

        fclose(fd);

        sprintf(name, "./spec_knots_%03i", Snap.SnapNum);
        fd = fopen(name, "w");
	    fprintf(fd,"#i idx, P[1], Mleft[0] Mleft[1] Mright[0] Mright[1] nCRe\n");
        for (int i = 0; i < Nknots; i++)
            fprintf(fd, "%d %i %g %g %g %g %g %g \n", 
                i, K[i].idx, K[i].P[1], K[i].Mleft[0], K[i].Mleft[1], 
                K[i].Mright[0], K[i].Mright[1], *nCRe);

        fclose(fd);
	
	    Uncompress(*nCRe,spectrum,np );
    }
#endif

    return;
}

void Uncompress(const float nCRe, const char *spectrum, double *out_spectrum)
{ 
    Assert(out_spectrum != NULL, "Can't uncompress into NULL");
    Assert(spectrum != NULL, "Can't uncompress from NULL");

    if (nCRe == 0)
        return;

	for (int i=0; i<N; i++)
		dR[i] = dL[i] = ddR[i] = ddL[i] = 0;

    struct Knot K[MAXKNOTS] = { 0 };
    
    int nKnots = uncompress_knots_binary(spectrum, K);

    double np[N_SPEC_BINS] = { -FLT_MAX };

    draw_curve(K, np); // fill raw spectrum

#ifdef DEBUG_COMPRESSION
	
	//check_knots(K, np);
    
	if (Ipart == IPART) {
        char name[128];
        FILE *fd;
        sprintf(name, "./spec_compressed_%03i", Snap.SnapNum);
        fd = fopen(name, "w");
        for (int i = 0; i < N; i++)
            fprintf(fd, "%d %g %g  %g \n", i, log_p[i], pow(10, np[i]*nCRe) , np[i]);

        fclose(fd);
    }
#endif

    for (int i = 0; i < N_SPEC_BINS; i++)  
        out_spectrum[i] = pow(10, np[i] * nCRe);

	return;
}

void Setup_compression()  
{
    for (int i=0; i<N_SPEC_BINS; i++ )
        log_p[i] = log10(p[i]/m_e/c);

    for (int i=1; i<N_SPEC_BINS; i++ )
        log_dp[i] = log_p[i]-log_p[i-1];
    	
	rprintf("Using Hermite spline compression of spectrum\n"
			"  Max total size          : %d bytes \n"
			"  Full knot size          : %zu bytes \n"
			"  Half knot size          : %zu bytes \n"
			"  Max number of knots     : %zu \n"
			"  Typical number of knots : %g \n",
			SPECSIZE_BYTES, SIZEOF_FULLKNOT, SIZEOF_HALFKNOT, MAXKNOTS, 
			5.0 + floor((SPECSIZE_BYTES-5.0*SIZEOF_HALFKNOT ) 
                / SIZEOF_FULLKNOT));
	return;
}

/* find integral, take log10, clean, find idx of maximum, normalise 
 * & calc derivatives */
double compute_spec_parameters(double *np)
{
    int i, global_max_idx = 0;

    double integral = 0;

    for (i = LowIdx; i < HighIdx; i++) { // in the boundary regions is garbage 

        integral += np[i] * dp[i];

        np[i] = fmax(log10(np[i]), -DBL_MAX);

        if ( ! isfinite(np[i]) ) 
            np[i] = 0;

        if (np[i] > np[global_max_idx])
            global_max_idx = i;
    }

	if (global_max_idx == LowIdx) // leave room for second derivative
		global_max_idx++;

    const double norm = np[global_max_idx];
    
    Assert(isfinite(norm) || norm == 0, "Spec norm is not finite or 0");
    
    for (i = LowIdx-1; i < HighIdx; i++) { // apply norm, compute derivatives
        
        np[i+1] /= norm; // normalise, global max == 1 

        dL[i+1] = (np[i+1]-np[i])/(0.5*(log_dp[i]+log_dp[i+1]));
        dR[i-1] = dL[i];

        ddL[i+1] = (np[i+1]-2*np[i]+np[i-1])/(log_dp[i]+log_dp[i+1]);
        ddR[i-1] = ddL[i];
    }
 
    for (i = 0; i < LowIdx; i++) { // treat boundaries
        np[i] =  -FLT_MAX;
        dL[i] = dR[i] = ddL[i] = ddR[i] = 0;
    }

    for (i = HighIdx; i < N; i++) { 
        np[i] =  -FLT_MAX;
        dL[i] = dR[i] = ddL[i] = ddR[i] = 0;
    }

#ifdef DEBUG_COMPRESSION
    if (Ipart == IPART) 
	    printf("Setting  maxidx=%d, norm=%g integral=%g \n", 
            global_max_idx, norm, integral );
#endif

    return norm;
}

/* find Knot points from start, stop, min/max, max global 
 * we leave one cell space to the boundary at start and stop so the 
 * derivaties are well defined */
void find_first_knots(const double *np, struct Knot *K)
{
    const double threshold = 0.1;

    int lastKnot = 0;
    int nMinMax = 0;
    int iMindL = 1;
    
    for (int i = LowIdx; i < HighIdx; i++) {
        
        enum knot_types type = KNOT_EMPTY;

        if (Nknots == 0) { // find start knot first

            if (np[i] > threshold || np[i] == 1 ) {
                
				type = KNOT_START;
	
				ddR[i] = ddR[i+1]; // extrapolate
				ddL[i] = ddL[i+1];

			} else        
                continue; 
	    }

        int d_last = i - lastKnot;

		if (dL[i]*dR[i] <= 0 && d_last >= KNOT_MIN_DIST && i != HighIdx-1) {

            type = KNOT_MINMAX;
                
            nMinMax++;

        } else if (np[i+1] < threshold || dL[i+1] < -4 || i == HighIdx-1) {

            type = KNOT_STOP;

			ddL[i] = ddL[i-1]; // extrapolate
			ddR[i] = ddR[i-1];

		}
		    
        if (fabs(dL[i]) < 0.01 && type == KNOT_MINMAX) 
            iMindL = i; // add other minmax knot here if needed

        if (type == KNOT_EMPTY) 
            continue; // nothing to add here
            
        add_knot(np, K, i, type);
    
        if (type == KNOT_STOP)
            break;    

        lastKnot = i;
    }
            
    if ( (nMinMax % 2) == 0  && nMinMax != 0) 
        add_knot(np, K, iMindL, KNOT_MINMAX); // we miss a MINMAX knot

    return;
}

void add_knot(const double * np, struct Knot * K, int idx, 
		enum knot_types type)
{
    int thisKnot = Nknots;

    switch (type) { 
	
        case KNOT_EMPTY:
            
            Assert(0, "Can't add empty knot");

            break;

        case KNOT_START:  
        case KNOT_STOP:
            
            KnotMemSize += SIZEOF_HALFKNOT;

            break;

        case KNOT_MINMAX: // impose assumptions
            
            dL[idx] = dR[idx] = 0;
            
             for (thisKnot = 0; thisKnot < Nknots; thisKnot++)
               if (K[thisKnot].idx > idx )              
                    break;

            size_t nBytes = (Nknots-thisKnot) * sizeof(*K); // make space 
            memmove(&(K[thisKnot+1]),&(K[thisKnot]), nBytes);

            KnotMemSize += SIZEOF_HALFKNOT;
            
            break;

        case KNOT_FULL: // insert 
            
            for (thisKnot = 0; thisKnot < Nknots; thisKnot++)
               if (K[thisKnot].idx > idx )              
                    break;

            nBytes = (Nknots-thisKnot) * sizeof(*K); // make space 
            memmove(&(K[thisKnot+1]),&(K[thisKnot]), nBytes);
            
            K[thisKnot].Mleft[0] = 0; // mark for recalculation of M 
            K[thisKnot].Mright[0] = 0;

            KnotMemSize += SIZEOF_FULLKNOT;
            
            break;

        default:

            Assert(0, "The bug you just found is awesome");
            
            break;
    }

	if (np[idx] == 1)
		K[thisKnot].Is_Global_Max = true;

    K[thisKnot].type = type;
    
    K[thisKnot].idx = idx;
    
    K[thisKnot].P[0] = log_p[idx];
    K[thisKnot].P[1] = np[idx];

    Nknots++;

#ifdef DEBUG_COMPRESSION
    if (Ipart == IPART) 
        printf("Adding Knot type %d, idx %d, p %g, np %g dL %g dR %g "
                "ddL %g ddR %g %d \n", type, idx, log_p[idx], np[idx], 
                dL[idx], dR[idx], ddL[idx], ddR[idx], K[thisKnot].Is_Global_Max);
#endif

    return;
}

/* fills missing / updates changed control points */
void update_control_points(struct Knot *K) 
{
    int next = 0;

    for (int this = 0; this < Nknots-1; this++) {
        
        next++;

        if( (K[next].Mleft[0]) && (K[this].Mright[0]))
            continue;       // nothing to do 

        int thisIdx = K[this].idx;
        int nextIdx = K[next].idx;
        
        /* x from second derivative, y from first */
        K[this].Mright[0] = -1.0/6.0 * ddL[nextIdx] 
            - 1.0/3.0*ddR[thisIdx] - K[this].P[0] + K[next].P[0];        
        K[next].Mleft[0] = 1.0/6.0 * ddR[thisIdx] 
            + 1.0/3.0*ddL[nextIdx] - K[this].P[0] + K[next].P[0];

        K[this].Mright[1] = K[this].Mright[0] * dR[thisIdx];
        K[next].Mleft[1] = K[next].Mleft[0] * dL[nextIdx];

#ifdef DEBUG_COMPRESSION
    if (Ipart == IPART) {
        printf(" Cor Knot %d M0 %d Mr0=%g Mr1=%g ddR=%g dR=%g \n", 
				this,thisIdx, K[this].Mright[0],K[this].Mright[1], 
				ddR[thisIdx], dR[thisIdx]);
        printf(" Cor Knot %d M1 %d Ml0=%g Ml1=%g ddL=%g dL=%g \n", 
				next,nextIdx, K[next].Mleft[0],K[next].Mleft[1],
				ddL[nextIdx], dL[nextIdx]);
    }
#endif
    }

    return;
}

/* find y(x), cubic Hermite spline */
void draw_curve(const struct Knot *K, double *Curve) 
{
    double c0[2], c1[2], c2[2], c3[2];
    
	for (int i=0; i < K[0].idx; i++)
		Curve[i] = -DBL_MAX;
    
	for (int i = 0; i < Nknots-1; i++) { // def. Hermitian polynome (cspline)
        
        c0[0] = 2*K[i].P[0] - 2*K[i+1].P[0] + K[i].Mright[0] + K[i+1].Mleft[0];
        c0[1] = 2*K[i].P[1] - 2*K[i+1].P[1] + K[i].Mright[1] + K[i+1].Mleft[1];

        c1[0] = -3*K[i].P[0]+3*K[i+1].P[0]-2*K[i].Mright[0] - K[i+1].Mleft[0];
        c1[1] = -3*K[i].P[1]+3*K[i+1].P[1]-2*K[i].Mright[1] - K[i+1].Mleft[1];

        c2[0] = K[i].Mright[0];
        c2[1] = K[i].Mright[1];

        c3[0] = K[i].P[0];
        c3[1] = K[i].P[1];

#ifdef DEBUG_COMPRESSION
    if (Ipart == IPART) 
        printf("Curve: %d@%d | P=%g,%g | Mr=%g,%g Ml=%g,%g | "
				"c0=%g c1=%g c2=%g c3=%g\n",
				i, K[i].idx, K[i].P[0],K[i].P[1],K[i].Mright[0],K[i].Mright[1]
				, K[i].Mleft[0],K[i].Mleft[1] ,c0[1], c1[1], c2[1], c3[1]  );
#endif

        for (int j = K[i].idx; j < K[i+1].idx+1; j++) {

			float t = 0, q = 0;
			int it = 0;

			do { // Newton-Raphson, damn it's fast 

				q = c0[0]*t*t*t + c1[0]*t*t + c2[0]*t + c3[0];
				t -= (q - log_p[j]) / (3*c0[0]*t*t + 2 * c1[0]*t + c2[0]);

			} while ( fabs(q-log_p[j]) / log_p[j] > 1e-2 && it++ < 50 );
	
			if(it > 49) 
				printf("Newton-Raphson failed ID=%d %d %g %g %g\n", 
						Ipart, j, log_p[j], q, t); 

			Curve[j] = c0[1]*t*t*t + c1[1]*t*t + c2[1]*t + c3[1];
        }
    }

	for (int i = K[Nknots-1].idx+1; i < N; i++)
		Curve[i] = -DBL_MAX;
    
	return;
}

/* returns max error between curve and data and its index */
float find_max_error(const struct Knot *K, const double *np, 
        double *curve, int *max_err_idx)
{
    float  dErr[N] = { 0 }, err[N] = { 0 }, max_err = 0;

    for (int i = K[0].idx; i < K[Nknots-1].idx; i++)
        err[i] = fabs(curve[i]-np[i]) / np[i] * log_p[i];

    for (int i = 0; i < Nknots; i++) {

		int jMin = K[i].idx + KNOT_MIN_DIST;
		int jMax = K[i+1].idx - KNOT_MIN_DIST;

    	for (int j = jMin; j < jMax; j++) {

        	//dErr[j] = (err[j+1]-err[j]) / (0.5 * (log_dp[j]+log_dp[j+1]));

        	//if ((dErr[j] * dErr[j-1] < 0) && (err[j] > max_err) ) {
        	if ( err[j] > max_err ) {
            	*max_err_idx = j;

            	max_err = err[j];
        	}
		}
        
        if (max_err == 0 && i == Nknots-1) { // refine the last bin

            *max_err_idx = jMin + (jMax-jMin) * 0.5;

            max_err = err[*max_err_idx];
        }
    }

    return max_err;
}

void check_knots(const struct Knot *K, double *spectrum)
{
    double curve[N] = { -DBL_MAX }; 
    
	draw_curve(K, curve);    

	int max_idx = K[Nknots-1].idx;
	int min_idx = 20;

	if (max_idx == 0)
		printf("CHECK: END, ID=%d N=%d \n", Ipart, Nknots-1);

	double err_max = 0;
	int i_err_max = -1;

	for (int i = min_idx; i < max_idx; i++) {
		
		double err = (curve[i]-spectrum[i])/spectrum[i];

		if (fabs(err) > 0.10) {
			err_max = err;
			i_err_max = i;
		}
	}

	if (err_max > 0.1)
		printf("CHECK: ID=%zu, idx=%d, err=%g, curve=%g, spec=%g \n", 
				Ipart, i_err_max, err_max, curve[i_err_max], spectrum[i_err_max]);
	
	for(int i = 0; i < Nknots; i++) {

		if (K[i].P[1] > 8 || K[i].P[1] < -8)
			printf("CHECK: P1=%g ID=%d idx=%d type=%d \n", 
					K[i].P[1], Ipart, K[i].idx);

		if (K[i].Mleft[1] > 8 || K[i].Mleft[1] < -8)
			printf("CHECK: Mleft1=%g ID=%d idx=%d type=%d\n", 
					K[i].Mleft[1], Ipart, K[i].idx, K[i].type);

		if (K[i].Mright[1] > 8 || K[i].Mright[1] < -8 && K[i].type != 2)
			printf("CHECK: Mright1=%g ID=%d idx=%d type=%d\n", 
					K[i].Mright[1], Ipart, K[i].idx, K[i].type);

		if (K[i].Mleft[0] > 2 || K[i].Mleft[0] < -2)
			printf("CHECK: Mleft0=%g ID=%d idx=%d type=%d\n", 
					K[i].Mleft[0], Ipart, K[i].idx, K[i].type);
		
		if ( (K[i].Mright[0] > 2 || K[i].Mright[0] < -2) && K[i].type != 2)
			printf("CHECK: Mright0=%g ID=%d idx=%d type=%d\n", 
					K[i].Mright[0], Ipart, K[i].idx, K[i].type);
	}

	return ;
}

void compress_knots_binary(const struct Knot *K, char *spectrum)
{
    char pkg[SPECSIZE_BYTES] = { 0 }; // data package
    struct half_point hp = { 0 };
    struct full_point fp = { 0 };

    size_t offset = 0;
	float dxdy = 0;
	
	hp.x = K[0].idx; // KNOT_START 
    hp.y = compressFloat_16bit(K[0].P[1]);
	
	hp.M.xyLR[0] = compressFloat_8bit(K[0].Mright[0]);
	hp.M.xyLR[1] = compressFloat_8bit(K[0].Mright[1]);

#ifdef DEBUG_COMPRESSION
    	if (Ipart == IPART) 
      		printf("Compress START :  x=%d y=%g, MR0=%g, MR1=%g \n", 
            	hp.x, K[0].P[1], K[0].Mright[0], K[0].Mright[1]);
#endif

	memcpy(pkg+offset, &hp.x, sizeof(hp) ); // offset == 0 obviously

    offset += SIZEOF_HALFKNOT;

    for (int i = 1; i < Nknots-1; i++) { // middle knots

        int8_t idx = K[i].idx;

        iii = K[i].idx;
            
        if (K[i].type == KNOT_FULL) { // KNOT_FULL
            
            fp.x = 0x80 | idx; // mark first bit to set KNOT_FULL
            fp.y = compressFloat_16bit(K[i].P[1]);
            
            fp.xleft = compressFloat_8bit(K[i].Mleft[0]);
            fp.yleft = compressFloat_16bit(K[i].Mleft[1]);
        
            fp.xright = compressFloat_8bit(K[i].Mright[0]);
            fp.yright = compressFloat_16bit(K[i].Mright[1]);

            memcpy(pkg+offset, &fp.x, sizeof(fp) );

            offset += SIZEOF_FULLKNOT;

#ifdef DEBUG_COMPRESSION
    	if (Ipart == IPART) 
            printf("Compress FULL:  x=%d y=%g xl=%g, yl=%g xr=%g, yr=%g \n", 
                    0x7F & fp.x,K[i].P[1],K[i].Mleft[0],K[i].Mleft[1],
                    K[i].Mright[0], K[i].Mright[1]);
#endif

        } else { // KNOT_MINMAX

            hp.x = K[i].idx;
            hp.y = compressFloat_16bit(K[i].P[1]);

            hp.M.xLR[0] = compressFloat_8bit(K[i].Mleft[0]);
            hp.M.xLR[1] = compressFloat_8bit(K[i].Mright[0]);

            memcpy(pkg+offset, &hp.x, sizeof(hp) );

            offset += SIZEOF_HALFKNOT;

#ifdef DEBUG_COMPRESSION
    		if (Ipart == IPART) 
        		printf("Compress MINMAX:  x=%d y=%g, xL=%g, xR=%g \n", 
                    hp.x,K[i].P[1],K[i].Mleft[0],K[i].Mright[0]);
#endif
        }
    }

    hp.x = K[Nknots-1].idx; // KNOT_STOP
    hp.y = compressFloat_16bit(K[Nknots-1].P[1]);

  	hp.M.xyLR[0] = compressFloat_8bit(K[Nknots-1].Mleft[0]);
	hp.M.xyLR[1] = compressFloat_8bit(K[Nknots-1].Mleft[1]);

	memcpy(pkg+offset, &hp.x, sizeof(hp) );
    offset += SIZEOF_HALFKNOT;

#ifdef DEBUG_COMPRESSION
    if (Ipart == IPART) 
    	printf("Compress STOP:  x=%d y=%g Ml0=%g Ml1=%g\n", 
				K[Nknots-1].idx,K[Nknots-1].P[1], K[0].Mleft[0],K[0].Mleft[1]);
#endif
        
    memcpy(spectrum, pkg, SPECSIZE_BYTES); // move to package

    return;
}

int uncompress_knots_binary(const char *spectrum, struct Knot *K)
{
    struct half_point hp = { 0 };
    struct full_point fp = { 0 };

    const char *src = spectrum;
        
    memcpy(&hp, src, SIZEOF_HALFKNOT); // KNOT_START
    src += SIZEOF_HALFKNOT;    

    K[0].type = KNOT_START;

    size_t idx = K[0].idx = hp.x;

    K[0].P[0] = log_p[idx];
    K[0].P[1] = uncompressFloat_16bit(hp.y);

	if (K[0].P[1] == 1)
		K[0].Is_Global_Max = true;

	K[0].Mright[0] = uncompressFloat_8bit(hp.M.xyLR[0]);
	K[0].Mright[1] = uncompressFloat_8bit(hp.M.xyLR[1]);

	Nknots = 1;           
    KnotMemSize = SIZEOF_HALFKNOT;

#ifdef DEBUG_COMPRESSION
    if (Ipart == IPART) 
        printf("DeCompress START, x=%d y=%g Mr[0]=%g MR[1] = %g \n", 
            K[0].idx, K[0].P[1], K[0].Mright[0], K[0].Mright[1] );
#endif

	int nMinMaxKnots = 0;

    for (int i = 1; i < MAXKNOTS; i++) {
        
        if (*src & 0x80) { // marker bit set
            
			K[i].type = KNOT_FULL;
            
            memcpy(&fp, src, SIZEOF_FULLKNOT);

            src += SIZEOF_FULLKNOT;
            KnotMemSize += SIZEOF_FULLKNOT;

            idx = K[i].idx = 0x7F & fp.x; // remove marker bit

            K[i].P[0] = log_p[idx];
            K[i].P[1] = uncompressFloat_16bit(fp.y);

            K[i].Mleft[0] = uncompressFloat_8bit(fp.xleft);
            K[i].Mleft[1] = uncompressFloat_16bit(fp.yleft);

            K[i].Mright[0] = uncompressFloat_8bit(fp.xright);
            K[i].Mright[1] = uncompressFloat_16bit(fp.yright);

#ifdef DEBUG_COMPRESSION
    if (Ipart == IPART) 
            printf("DeCompress FULL, x=%d y=%g xL=%g yL=%g, xR=%g yR=%g \n", 
                K[i].idx, K[i].P[1],K[i].Mleft[0],K[i].Mleft[1],K[i].Mright[0],
                K[i].Mright[1]);
#endif
        } else { // marker bit not set 
            
			K[i].type = KNOT_MINMAX;

            memcpy(&hp, src, SIZEOF_HALFKNOT);
            
            src += SIZEOF_HALFKNOT;
            KnotMemSize += SIZEOF_HALFKNOT;

            idx = K[i].idx = hp.x;

            K[i].P[0] = log_p[idx];
            K[i].P[1] = uncompressFloat_16bit(hp.y);

            K[i].Mleft[0] = uncompressFloat_8bit(hp.M.xLR[0]);
            K[i].Mleft[1] = 0;           

            K[i].Mright[0] = uncompressFloat_8bit(hp.M.xLR[1]);
            K[i].Mright[1] = 0;

			nMinMaxKnots++;
        
#ifdef DEBUG_COMPRESSION
    if (Ipart == IPART) 
            printf("DeCompress MINMAX, x=%d y=%g xL=%g yL=%g, xR=%g yR=%g \n", 
                K[i].idx, K[i].P[1],K[i].Mleft[0],K[i].Mleft[1],K[i].Mright[0],
                K[i].Mright[1]);
#endif
        } 
        
        Nknots++;           

        if (nMinMaxKnots % 2 
			&& KnotMemSize + SIZEOF_FULLKNOT+SIZEOF_HALFKNOT > SPECSIZE_BYTES ) 
            break; // only KNOT_STOP possible, exit loop
    }
    
    memcpy(&hp, src, SIZEOF_HALFKNOT); // KNOT_STOP

	if (hp.x == 0 ) {// we overshot, if START==MINMAX 
		
		src -= SIZEOF_HALFKNOT;

		memcpy(&hp, src, SIZEOF_HALFKNOT);
	}

    K[Nknots].type = KNOT_STOP;
    
    K[Nknots].idx = idx = hp.x;

    K[Nknots].P[0] = log_p[idx];
    K[Nknots].P[1] = uncompressFloat_16bit(hp.y);
     
	K[Nknots].Mleft[0] = uncompressFloat_8bit(hp.M.xyLR[0]);
	K[Nknots].Mleft[1] = uncompressFloat_8bit(hp.M.xyLR[1]);

#ifdef DEBUG_COMPRESSION
    if (Ipart == IPART) 
        printf("Decompress STOP, x=%d y=%g  ML0=%g ML1=%g \n", 
            K[Nknots].idx, K[Nknots].P[1],K[Nknots].Mleft[0],K[Nknots].Mleft[1]);
#endif

    Nknots++;           

    return Nknots;
}

/* float to fixed point float conversions, see Wikipedia. */ 
uint16_t compressFloat_16bit(float input) 
{
#ifdef COMPRESSION_RANGE_WARNING 
    if (input > 8 || input < -8) {

        rprintf("16bit Compression, input out of range %g %d %d\n", 
                input, Ipart, P[Ipart].ID );   
        Flag = 1;
    }
#endif

    input = fmax(-8, input);
    input = fmin(7.99, input); // rolls over at 8 and becomes negative

    return roundf( (input+8) * 4096); // * 2^12
}   

float uncompressFloat_16bit(uint16_t input)
{
    return input * 0.000244140625 - 8; // * 2^-12
}

uint8_t compressFloat_8bit(float input)
{
#ifdef COMPRESSION_RANGE_WARNING
    if (input > 3 || input < -1) {

        rprintf("8bit Compression, input out of range "
                "num=%g ipart=%d ID=%d idx=%d  \n",
                input , Ipart, P[Ipart].ID, iii);   
        Flag=1;
    }
#endif

    input = fmax(-1, input);
    input = fmin(2.99, input); // rolls over at 2 and becomes negative

    return roundf( (input+1) * 64);  // 2^6
}   

float uncompressFloat_8bit(uint8_t input)
{   
    return input * 0.015625 -1; // 2^-6
}

float constrainFloat_16bit(float input) // convert to next representable value
{
    return uncompressFloat_16bit(compressFloat_16bit(input));
}

float constrainFloat_8bit(float input)
{
    return uncompressFloat_8bit(compressFloat_8bit(input));
}

#undef N 
#undef MAXKNOTS 
#undef CUTOFF 
