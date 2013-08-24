/*
 * codelets4.h
 *
 *  Created on: Feb 21, 2013
 *      Author: USER
 */

#ifndef CODELETS4_H_
#define CODELETS4_H_

#include "codelets3.h"

void inline radix3_dit_rec(dblitr opr, dblitr opi,dblitr ipr, dblitr ipi, int sgn,int N,int l) {

	if (N == 1) {
		*opr = *ipr;
		*opi = *ipi;
		//return;
	} else {
		int m = N / 3;

		radix3_dit_rec(opr,opi,ipr,ipi,sgn,m,l*3);
		radix3_dit_rec(opr+m,opi+m,ipr+l,ipi+l,sgn,m,l*3);
		radix3_dit_rec(opr+2*m,opi+2*m,ipr+2*l,ipi+2*l,sgn,m,l*3);

		double PI2 = 6.28318530717958647692528676655900577;
		double theta =  -1.0 * sgn * PI2/N;

		double S =  sin(theta);
		double C =  cos(theta);

		double wlr =  1.0;
		double wli =  0.0;
		
		double wl2r =  1.0;
		double wl2i =  0.0;

		double tau0r,tau0i,tau1r,tau1i,tau2r,tau2i;
		double ar,ai,br,bi,cr,ci;
		int tkm1,tkm2;

		for (int k = 0; k < m; ++k) {
			ar = *(opr+k);
			ai = *(opi+k);
			tkm1 = k + m;
			tkm2 = tkm1 + m;

			br = *(opr+tkm1)*wlr - *(opi+tkm1)*wli;
			bi = *(opi+tkm1)*wlr + *(opr+tkm1)*wli;

			cr = *(opr+tkm2)*wl2r - *(opi+tkm2)*wl2i;
			ci = *(opi+tkm2)*wl2r + *(opr+tkm2)*wl2i;
			
			tau0r = br + cr;
			tau0i = bi + ci;
			
			tau1r = sgn * 0.86602540378 * (br - cr);
			tau1i = sgn * 0.86602540378 * (bi - ci);
			
			tau2r = ar - tau0r * 0.5000000000;
			tau2i = ai - tau0i * 0.5000000000; 										
								

			*(opr+k) = ar + tau0r ;
			*(opi+k) = ai + tau0i;
			
			*(opr+tkm1) = tau2r + tau1i;
			*(opi+tkm1) = tau2i - tau1r;
			
			*(opr+tkm2) = tau2r - tau1i;
			*(opi+tkm2) = tau2i + tau1r;			

			double temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			wl2r = wlr*wlr - wli*wli;
			wl2i = 2.0*wlr*wli;			

		}

	}
}

void inline radix4_dit_rec(dblitr opr, dblitr opi,dblitr ipr, dblitr ipi, int sgn,int N,int l) {
	if (N == 1) {
		*opr = *ipr;
		*opi = *ipi;
		//return;
	} else {
		int m = N / 4;

		radix4_dit_rec(opr,opi,ipr,ipi,sgn,m,l*4);
		radix4_dit_rec(opr+m,opi+m,ipr+l,ipi+l,sgn,m,l*4);
		radix4_dit_rec(opr+2*m,opi+2*m,ipr+2*l,ipi+2*l,sgn,m,l*4);
		radix4_dit_rec(opr+3*m,opi+3*m,ipr+3*l,ipi+3*l,sgn,m,l*4);

		double PI2 = 6.28318530717958647692528676655900577;
		double theta =  -1.0 * sgn * PI2/N;

		double S =  sin(theta);
		double C =  cos(theta);

		double wlr =  1.0;
		double wli =  0.0;
		
		double wl2r =  1.0;
		double wl2i =  0.0;
		
		double wl3r =  1.0;
		double wl3i =  0.0;

		double tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
		double ar,ai,br,bi,cr,ci,dr,di;
		int tkm1,tkm2,tkm3;

		for (int k = 0; k < m; ++k) {
			ar = *(opr+k);
			ai = *(opi+k);
			tkm1 = k + m;
			tkm2 = tkm1 + m;
			tkm3 = tkm2 + m;

			br = *(opr+tkm1)*wlr - *(opi+tkm1)*wli;
			bi = *(opi+tkm1)*wlr + *(opr+tkm1)*wli;

			cr = *(opr+tkm2)*wl2r - *(opi+tkm2)*wl2i;
			ci = *(opi+tkm2)*wl2r + *(opr+tkm2)*wl2i;
			
			dr = *(opr+tkm3)*wl3r - *(opi+tkm3)*wl3i;
			di = *(opi+tkm3)*wl3r + *(opr+tkm3)*wl3i;
			
			tau0r = ar + cr;
			tau0i = ai + ci;

			tau1r = ar - cr;
			tau1i = ai - ci;

			tau2r = br + dr;
			tau2i = bi + di;

			tau3r = sgn* (br - dr);
			tau3i = sgn* (bi - di);							
								

			*(opr+k) = tau0r + tau2r ;
			*(opi+k) = tau0i + tau2i;
			
			*(opr+tkm1) = tau1r + tau3i;
			*(opi+tkm1) = tau1i - tau3r;
			
			*(opr+tkm2) = tau0r - tau2r;
			*(opi+tkm2) = tau0i - tau2i;		
			
			*(opr+tkm3) = tau1r - tau3i;
			*(opi+tkm3) = tau1i + tau3r;		

			double temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			wl2r = wlr*wlr - wli*wli;
			wl2i = 2.0*wlr*wli;		
			
			wl3r = wl2r*wlr - wli*wl2i;
			wl3i= wl2r*wli + wl2i*wlr;	

		}

	}

}

void inline radix5_dit_rec(dblitr opr, dblitr opi,dblitr ipr, dblitr ipi, int sgn,int N,int l) {
	if (N == 1) {
		*opr = *ipr;
		*opi = *ipi;
		//return;
	} else {
		int m = N / 5;

		radix5_dit_rec(opr,opi,ipr,ipi,sgn,m,l*5);
		radix5_dit_rec(opr+m,opi+m,ipr+l,ipi+l,sgn,m,l*5);
		radix5_dit_rec(opr+2*m,opi+2*m,ipr+2*l,ipi+2*l,sgn,m,l*5);
		radix5_dit_rec(opr+3*m,opi+3*m,ipr+3*l,ipi+3*l,sgn,m,l*5);
		radix5_dit_rec(opr+4*m,opi+4*m,ipr+4*l,ipi+4*l,sgn,m,l*5);

		double PI2 = 6.28318530717958647692528676655900577;
		double theta =  -1.0 * sgn * PI2/N;

		double S =  sin(theta);
		double C =  cos(theta);

		double wlr =  1.0;
		double wli =  0.0;
		
		double wl2r =  1.0;
		double wl2i =  0.0;
		
		double wl3r =  1.0;
		double wl3i =  0.0;

		double wl4r =  1.0;
	    double wl4i =  0.0;

		double c1 = 0.30901699437;
		double c2 = -0.80901699437;
		double s1 = 0.95105651629;
		double s2 = 0.58778525229;

		double tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
		double tau4r,tau4i,tau5r,tau5i;
		double ar,ai,br,bi,cr,ci,dr,di,er,ei;
		
		int tkm1,tkm2,tkm3,tkm4;

		for (int k = 0; k < m; ++k) {
			ar = *(opr+k);
			ai = *(opi+k);
			tkm1 = k + m;
			tkm2 = tkm1 + m;
			tkm3 = tkm2 + m;
			tkm4 = tkm3 + m;

			br = *(opr+tkm1)*wlr - *(opi+tkm1)*wli;
			bi = *(opi+tkm1)*wlr + *(opr+tkm1)*wli;

			cr = *(opr+tkm2)*wl2r - *(opi+tkm2)*wl2i;
			ci = *(opi+tkm2)*wl2r + *(opr+tkm2)*wl2i;
			
			dr = *(opr+tkm3)*wl3r - *(opi+tkm3)*wl3i;
			di = *(opi+tkm3)*wl3r + *(opr+tkm3)*wl3i;
			
			er = *(opr+tkm4)*wl4r - *(opi+tkm4)*wl4i;
			ei = *(opi+tkm4)*wl4r + *(opr+tkm4)*wl4i;
			
			tau0r = br + er;
			tau0i = bi + ei;

			tau1r = cr + dr;
			tau1i = ci + di;

			tau2r = br - er;
			tau2i = bi - ei;

			tau3r = cr - dr;
			tau3i = ci - di;					
								

			*(opr+k) = ar + tau0r + tau1r;
			*(opi+k) = ai + tau0i + tau1i;
			
			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;

			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);			
			
			*(opr+tkm1) = ar + tau4r + tau5i;
			*(opi+tkm1) = ai + tau4i - tau5r;
			
			*(opr+tkm4) = ar + tau4r - tau5i;
			*(opi+tkm4) = ai + tau4i + tau5r;	
			
			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;

			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);			
			
			*(opr+tkm2) = ar + tau4r + tau5i;
			*(opi+tkm2) = ai + tau4i - tau5r;

			*(opr+tkm3) = ar + tau4r - tau5i;
			*(opi+tkm3) = ai + tau4i + tau5r;	

			double temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			wl2r = wlr*wlr - wli*wli;
			wl2i = 2.0*wlr*wli;		
			
			wl3r = wl2r*wlr - wli*wl2i;
			wl3i= wl2r*wli + wl2i*wlr;	
			
			wl4r = wl2r*wl2r - wl2i*wl2i;
			wl4i = 2.0*wl2r*wl2i;

		}

	}

}

void inline radix7_dit_rec(dblitr opr, dblitr opi,dblitr ipr, dblitr ipi, int sgn,int N,int l) {
	if (N == 1) {
		*opr = *ipr;
		*opi = *ipi;
		//return;
	} else {
		int m = N / 7;

		radix7_dit_rec(opr,opi,ipr,ipi,sgn,m,l*7);
		radix7_dit_rec(opr+m,opi+m,ipr+l,ipi+l,sgn,m,l*7);
		radix7_dit_rec(opr+2*m,opi+2*m,ipr+2*l,ipi+2*l,sgn,m,l*7);
		radix7_dit_rec(opr+3*m,opi+3*m,ipr+3*l,ipi+3*l,sgn,m,l*7);
		radix7_dit_rec(opr+4*m,opi+4*m,ipr+4*l,ipi+4*l,sgn,m,l*7);
		radix7_dit_rec(opr+5*m,opi+5*m,ipr+5*l,ipi+5*l,sgn,m,l*7);
		radix7_dit_rec(opr+6*m,opi+6*m,ipr+6*l,ipi+6*l,sgn,m,l*7);

		double PI2 = 6.28318530717958647692528676655900577;
		double theta =  -1.0 * sgn * PI2/N;

		double S =  sin(theta);
		double C =  cos(theta);

		double wlr =  1.0;
		double wli =  0.0;
		
		double wl2r =  1.0;
		double wl2i =  0.0;
		
		double wl3r =  1.0;
		double wl3i =  0.0;

		double wl4r =  1.0;
	    double wl4i =  0.0;
	    
		double wl5r =  1.0;
	    double wl5i =  0.0;
	    
		double wl6r =  1.0;
	    double wl6i =  0.0;	    	    	    

		double c1 = 0.62348980185;
		double c2 = -0.22252093395;
		double c3 = -0.9009688679;
		double s1 = 0.78183148246;
		double s2 = 0.97492791218;
		double s3 = 0.43388373911;

		double tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
		double tau4r,tau4i,tau5r,tau5i;
		double tau6r,tau6i,tau7r,tau7i;
		double ar,ai,br,bi,cr,ci,dr,di,er,ei;
		double fr,fi,gr,gi;
		
		int tkm1,tkm2,tkm3,tkm4,tkm5,tkm6;

		for (int k = 0; k < m; ++k) {
			ar = *(opr+k);
			ai = *(opi+k);
			tkm1 = k + m;
			tkm2 = tkm1 + m;
			tkm3 = tkm2 + m;
			tkm4 = tkm3 + m;
			tkm5 = tkm4 + m;
			tkm6 = tkm5 + m;

			br = *(opr+tkm1)*wlr - *(opi+tkm1)*wli;
			bi = *(opi+tkm1)*wlr + *(opr+tkm1)*wli;

			cr = *(opr+tkm2)*wl2r - *(opi+tkm2)*wl2i;
			ci = *(opi+tkm2)*wl2r + *(opr+tkm2)*wl2i;
			
			dr = *(opr+tkm3)*wl3r - *(opi+tkm3)*wl3i;
			di = *(opi+tkm3)*wl3r + *(opr+tkm3)*wl3i;
			
			er = *(opr+tkm4)*wl4r - *(opi+tkm4)*wl4i;
			ei = *(opi+tkm4)*wl4r + *(opr+tkm4)*wl4i;
			
			fr = *(opr+tkm5)*wl5r - *(opi+tkm5)*wl5i;
			fi = *(opi+tkm5)*wl5r + *(opr+tkm5)*wl5i;
			
			gr = *(opr+tkm6)*wl6r - *(opi+tkm6)*wl6i;
			gi = *(opi+tkm6)*wl6r + *(opr+tkm6)*wl6i;
			
			tau0r = br + gr;
			tau0i = bi + gi;

			tau1r = cr + fr;
			tau1i = ci + fi;

			tau2r = dr + er;
			tau2i = di + ei;

			tau3r = br - gr;
			tau3i = bi - gi;

			tau4r = cr - fr;
			tau4i = ci - fi;

			tau5r = dr - er;
			tau5i = di - ei;				
								

			*(opr+k) = ar + tau0r + tau1r + tau2r;
			*(opi+k) = ai + tau0i + tau1i + tau2i;
			
			tau6r = ar + c1 * tau0r + c2 * tau1r + c3 * tau2r;
			tau6i = ai + c1 * tau0i + c2 * tau1i + c3 * tau2i;

			tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
			tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);		
			
			*(opr+tkm1) = tau6r - tau7i;
			*(opi+tkm1) = tau6i + tau7r;
			
			*(opr+tkm6) = tau6r + tau7i;
			*(opi+tkm6) = tau6i - tau7r;	
			
			tau6r = ar + c2 * tau0r + c3 * tau1r + c1 * tau2r;
			tau6i = ai + c2 * tau0i + c3 * tau1i + c1 * tau2i;

			tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
			tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

			*(opr+tkm2) = tau6r - tau7i;
			*(opi+tkm2) = tau6i + tau7r;

			*(opr+tkm5) = tau6r + tau7i;
			*(opi+tkm5) = tau6i - tau7r;
			
			tau6r = ar + c3 * tau0r + c1 * tau1r + c2 * tau2r;
			tau6i = ai + c3 * tau0i + c1 * tau1i + c2 * tau2i;

			tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
			tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

			*(opr+tkm3) = tau6r - tau7i;
			*(opi+tkm3) = tau6i + tau7r;

			*(opr+tkm4) = tau6r + tau7i;
			*(opi+tkm4) = tau6i - tau7r;

			double temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			wl2r = wlr*wlr - wli*wli;
			wl2i = 2.0*wlr*wli;		
			
			wl3r = wl2r*wlr - wli*wl2i;
			wl3i= wl2r*wli + wl2i*wlr;	
			
			wl4r = wl2r*wl2r - wl2i*wl2i;
			wl4i = 2.0*wl2r*wl2i;
			
			wl5r = wl3r*wl2r - wl2i*wl3i;
			wl5i= wl3r*wl2i + wl3i*wl2r;

			wl6r = wl3r*wl3r - wl3i*wl3i;
			wl6i = 2.0*wl3r*wl3i;

		}

	}

}

void inline radix8_dit_rec(dblitr opr, dblitr opi,dblitr ipr, dblitr ipi, int sgn,int N,int l) {
	if (N == 1) {
		*opr = *ipr;
		*opi = *ipi;
		//return;
	} else {
		int m = N / 8;

		radix8_dit_rec(opr,opi,ipr,ipi,sgn,m,l*8);
		radix8_dit_rec(opr+m,opi+m,ipr+l,ipi+l,sgn,m,l*8);
		radix8_dit_rec(opr+2*m,opi+2*m,ipr+2*l,ipi+2*l,sgn,m,l*8);
		radix8_dit_rec(opr+3*m,opi+3*m,ipr+3*l,ipi+3*l,sgn,m,l*8);
		radix8_dit_rec(opr+4*m,opi+4*m,ipr+4*l,ipi+4*l,sgn,m,l*8);
		radix8_dit_rec(opr+5*m,opi+5*m,ipr+5*l,ipi+5*l,sgn,m,l*8);
		radix8_dit_rec(opr+6*m,opi+6*m,ipr+6*l,ipi+6*l,sgn,m,l*8);
		radix8_dit_rec(opr+7*m,opi+7*m,ipr+7*l,ipi+7*l,sgn,m,l*8);

		double PI2 = 6.28318530717958647692528676655900577;
		double theta =  -1.0 * sgn * PI2/N;

		double S =  sin(theta);
		double C =  cos(theta);

		double wlr =  1.0;
		double wli =  0.0;
		
		double wl2r =  1.0;
		double wl2i =  0.0;
		
		double wl3r =  1.0;
		double wl3i =  0.0;

		double wl4r =  1.0;
	    double wl4i =  0.0;
	    
		double wl5r =  1.0;
	    double wl5i =  0.0;
	    
		double wl6r =  1.0;
	    double wl6i =  0.0;	    	    	
	    
	    double wl7r =  1.0;
	    double wl7i =  0.0;        

		double c1 = 0.70710678118654752440084436210485;
		//T c2 = -0.70710678118654752440084436210485;
		double s1 = 0.70710678118654752440084436210485;
		//T s2 = 0.70710678118654752440084436210485;

		double tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
		double tau4r,tau4i,tau5r,tau5i;
		double tau6r,tau6i,tau7r,tau7i;
		double tau8r,tau8i,tau9r,tau9i;
		double ar,ai,br,bi,cr,ci,dr,di,er,ei;
		double fr,fi,gr,gi,hr,hi;
		double temp1r,temp1i,temp2r,temp2i;
		
		int tkm1,tkm2,tkm3,tkm4,tkm5,tkm6,tkm7;

		for (int k = 0; k < m; ++k) {
			ar = *(opr+k);
			ai = *(opi+k);
			tkm1 = k + m;
			tkm2 = tkm1 + m;
			tkm3 = tkm2 + m;
			tkm4 = tkm3 + m;
			tkm5 = tkm4 + m;
			tkm6 = tkm5 + m;
			tkm7 = tkm6 + m;

			br = *(opr+tkm1)*wlr - *(opi+tkm1)*wli;
			bi = *(opi+tkm1)*wlr + *(opr+tkm1)*wli;

			cr = *(opr+tkm2)*wl2r - *(opi+tkm2)*wl2i;
			ci = *(opi+tkm2)*wl2r + *(opr+tkm2)*wl2i;
			
			dr = *(opr+tkm3)*wl3r - *(opi+tkm3)*wl3i;
			di = *(opi+tkm3)*wl3r + *(opr+tkm3)*wl3i;
			
			er = *(opr+tkm4)*wl4r - *(opi+tkm4)*wl4i;
			ei = *(opi+tkm4)*wl4r + *(opr+tkm4)*wl4i;
			
			fr = *(opr+tkm5)*wl5r - *(opi+tkm5)*wl5i;
			fi = *(opi+tkm5)*wl5r + *(opr+tkm5)*wl5i;
			
			gr = *(opr+tkm6)*wl6r - *(opi+tkm6)*wl6i;
			gi = *(opi+tkm6)*wl6r + *(opr+tkm6)*wl6i;
			
			hr = *(opr+tkm7)*wl7r - *(opi+tkm7)*wl7i;
			hi = *(opi+tkm7)*wl7r + *(opr+tkm7)*wl7i;
			
			tau0r = ar + er;
			tau0i = ai + ei;

			tau1r = br + hr;
			tau1i = bi + hi;

			tau2r = dr + fr;
			tau2i = di + fi;

			tau3r = cr + gr;
			tau3i = ci + gi;

			tau4r = ar - er;
			tau4i = ai - ei;

			tau5r = br - hr;
			tau5i = bi - hi;

			tau6r = dr - fr;
			tau6i = di - fi;

			tau7r = cr - gr;
			tau7i = ci - gi;				
								

			*(opr+k) = tau0r + tau1r + tau2r + tau3r;
			*(opi+k) = tau0i + tau1i + tau2i + tau3i;
			
			*(opr+tkm4) = tau0r - tau1r - tau2r + tau3r;
			*(opi+tkm4) = tau0i - tau1i - tau2i + tau3i;
			
			temp1r = tau1r - tau2r;
			temp1i = tau1i - tau2i;

			temp2r = tau5r + tau6r;
			temp2i = tau5i + tau6i;

			tau8r =  tau4r + c1 * temp1r;
			tau8i =  tau4i + c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r - tau7r);
			tau9i = sgn * ( -s1 * temp2i - tau7i);	
			
			*(opr+tkm1) = tau8r - tau9i;
			*(opi+tkm1) = tau8i + tau9r;

			*(opr+tkm7) = tau8r + tau9i;
			*(opi+tkm7) = tau8i - tau9r;	
			
			tau8r = tau0r - tau3r;
			tau8i = tau0i - tau3i;

			tau9r = sgn * ( -tau5r + tau6r);
			tau9i = sgn * ( -tau5i + tau6i);

			*(opr+tkm2) = tau8r - tau9i;
			*(opi+tkm2) = tau8i + tau9r;

			*(opr+tkm6) = tau8r + tau9i;
			*(opi+tkm6) = tau8i - tau9r;
			
			tau8r = tau4r - c1 * temp1r;
			tau8i = tau4i - c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r + tau7r);
			tau9i = sgn * ( -s1 * temp2i + tau7i);

			*(opr+tkm3) = tau8r - tau9i;
			*(opi+tkm3) = tau8i + tau9r;

			*(opr+tkm5) = tau8r + tau9i;
			*(opi+tkm5) = tau8i - tau9r;

			double temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			wl2r = wlr*wlr - wli*wli;
			wl2i = 2.0*wlr*wli;		
			
			wl3r = wl2r*wlr - wli*wl2i;
			wl3i= wl2r*wli + wl2i*wlr;	
			
			wl4r = wl2r*wl2r - wl2i*wl2i;
			wl4i = 2.0*wl2r*wl2i;
			
			wl5r = wl3r*wl2r - wl2i*wl3i;
			wl5i= wl3r*wl2i + wl3i*wl2r;

			wl6r = wl3r*wl3r - wl3i*wl3i;
			wl6i = 2.0*wl3r*wl3i;
			
			wl7r = wl3r*wl4r - wl4i*wl3i;
			wl7i= wl3r*wl4i + wl3i*wl4r;

		}

	}

}

void inline split_radix_rec2(dblitr r, dblitr i, int sgn,int N) {
	if (N == 1) {
		return;
	} else if (N == 2) {

		double taur =  *(r+1);
		double taui =  *(i+1);
			
		*(r+1) = *r - taur; 
		*(i+1) = *i - taui; 
			
		*r = *r + taur; 
		*i = *i + taui; 

		N=1;
	} else {
		int m = N/2;
		int p = N/4;
		
		double PI2 = 6.28318530717958647692528676655900577;
		double theta =  -1.0 * sgn * PI2/N;
		
		double S =  sin(theta);
		double C =  cos(theta);
		
		double PI6 =  3.0*6.28318530717958647692528676655900577;
		double theta3 =  -1.0 * sgn * PI6/N;
		
		double S3 =  sin(theta3);
		double C3 =  cos(theta3);
		
		double wlr =  1.0;
		double wli =  0.0;
		
		//T wl2r = (T) 1.0;
		//T wl2i = (T) 0.0;
		
		double wl3r =  1.0;
		double wl3i =  0.0;
		
		double tau1r,tau1i,tau2r,tau2i;
		double ur,ui,vr,vi;
		
		for (int j = 0; j < p; j++) {
			wlr =  cos(theta*j);
		    wli =  sin(theta*j);
			
			wl3r =  cos(theta3*j);
		    wl3i =  sin(theta3*j);
			
			int index1 = j+m;
			int index2 = index1+p;
			int index3 = j+p;
			
			tau1r = *(r+index1)*wlr - *(i+index1)*wli;
			tau1i = *(i+index1)*wlr + *(r+index1)*wli;
						
			tau2r = *(r+index2)*wl3r - *(i+index2)*wl3i;;
			tau2i = *(i+index2)*wl3r + *(r+index2)*wl3i;
						
			ur = tau1r + tau2r;
			ui = tau1i + tau2i;
						
			vr = sgn* (tau2r - tau1r);
			vi = sgn* (tau2i - tau1i);
						
			*(r+index2) = *(r+index3) - vi;
			*(i+index2) = *(i+index3) + vr;
						
			*(r+index1) = *(r+j) - ur;
			*(i+index1) = *(i+j) - ui;
						
			*(r+index3) = *(r+index3) + vi;
			*(i+index3) = *(i+index3) - vr;
						
			*(r+j) = *(r+j) + ur;
			*(i+j) = *(i+j) + ui;
			
			double temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			temp = wl3r;
			wl3r = C3 * wl3r - S3 * wl3i;
			wl3i = S3 * temp + C3 * wl3i;
			 
			
		}
		
		split_radix_rec2(r,i,sgn,m);
		split_radix_rec2(r+m,i+m,sgn,p);
		split_radix_rec2(r+m+p,i+m+p,sgn,p);
	}
	
	
}



#endif 
