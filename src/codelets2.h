/*
 * codelets2.h
 *
 *  Created on: Feb 12, 2013
 *      Author: USER
 */

#ifndef CODELETS2_H_
#define CODELETS2_H_

#include "codelets.h"

template <typename T>
void inline radix2_mixed_dit_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 2;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T taur,taui;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index1 = k*L+j;
			int index = index1+Ls;
			taur = wlr*x.re[index] - wli*x.im[index];
			taui = wlr*x.im[index] + wli*x.re[index];

			x.re[index] = x.re[index1] - taur;
			x.im[index] = x.im[index1] - taui;

			x.re[index1] = x.re[index1] + taur;
			x.im[index1] = x.im[index1] + taui;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;
	}

}

template <typename T>
void inline radix3_mixed_dit_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 3;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i;
	T ar,ai,br,bi,cr,ci;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];



			tau0r = br + cr;
			tau0i = bi + ci;

			tau1r = sgn * 0.86602540378 * (br - cr);
			tau1i = sgn * 0.86602540378 * (bi - ci);

			x.re[index] = ar + tau0r;
			x.im[index] = ai + tau0i;

			tau2r = ar - tau0r * 0.5000000000;
			tau2i = ai - tau0i * 0.5000000000;

			x.re[index1] = tau2r + tau1i;
			x.im[index1] = tau2i - tau1r;

			x.re[index2] = tau2r - tau1i;
			x.im[index2] = tau2i + tau1r;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

	}

}

template <typename T>
void inline radix4_mixed_dit_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 4;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T ar,ai,br,bi,cr,ci,dr,di;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];

			dr = wl3r*x.re[index3] - wl3i*x.im[index3];
			di = wl3r*x.im[index3] + wl3i*x.re[index3];

			tau0r = ar + cr;
			tau0i = ai + ci;

			tau1r = ar - cr;
			tau1i = ai - ci;

			tau2r = br + dr;
			tau2i = bi + di;

			tau3r = sgn* (br - dr);
			tau3i = sgn* (bi - di);

			x.re[index] = tau0r + tau2r;
			x.im[index] = tau0i + tau2i;

			x.re[index1] = tau1r + tau3i;
			x.im[index1] = tau1i - tau3r;

			x.re[index2] = tau0r - tau2r;
			x.im[index2] = tau0i - tau2i;

			x.re[index3] = tau1r - tau3i;
			x.im[index3] = tau1i + tau3r;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;
	}

}

template <typename T>
void inline radix5_mixed_dit_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 5;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T c1 = 0.30901699437;
	T c2 = -0.80901699437;
	T s1 = 0.95105651629;
	T s2 = 0.58778525229;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];

			dr = wl3r*x.re[index3] - wl3i*x.im[index3];
			di = wl3r*x.im[index3] + wl3i*x.re[index3];

			er = wl4r*x.re[index4] - wl4i*x.im[index4];
			ei = wl4r*x.im[index4] + wl4i*x.re[index4];

			tau0r = br + er;
			tau0i = bi + ei;

			tau1r = cr + dr;
			tau1i = ci + di;

			tau2r = br - er;
			tau2i = bi - ei;

			tau3r = cr - dr;
			tau3i = ci - di;

			x.re[index] = ar + tau0r + tau1r;
			x.im[index] = ai + tau0i + tau1i;

			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;

			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

			x.re[index1] = ar + tau4r + tau5i;
			x.im[index1] = ai + tau4i - tau5r;

			x.re[index4] = ar + tau4r - tau5i;
			x.im[index4] = ai + tau4i + tau5r;

			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;

			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

			x.re[index2] = ar + tau4r + tau5i;
			x.im[index2] = ai + tau4i - tau5r;

			x.re[index3] = ar + tau4r - tau5i;
			x.im[index3] = ai + tau4i + tau5r;


		}
		T temp = wlr;
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

template <typename T>
void inline radix2_mixed_dit_inplace(fft_data<T> &x,fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 2;
	int r = n / L;

	T wlr,wli;

	T taur,taui;

	for (int j = 0; j < Ls; j++) {
		int ind = r*j;
		wlr = wl.re[ind];
		wli = wl.im[ind];
		for (int k =0; k < r; k++) {
			int index1 = k*L+j;
			int index = index1+Ls;
			taur = wlr*x.re[index] - wli*x.im[index];
			taui = wlr*x.im[index] + wli*x.re[index];

			x.re[index] = x.re[index1] - taur;
			x.im[index] = x.im[index1] - taui;

			x.re[index1] = x.re[index1] + taur;
			x.im[index1] = x.im[index1] + taui;


		}

	}

}


template <typename T>
void inline radix3_mixed_dit_inplace(fft_data<T> &x,fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 3;
	int r = n / L;


	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i;
	T ar,ai,br,bi,cr,ci;
	T wlr,wli,wl2r,wl2i;

	for (int j = 0; j < Ls; j++) {
		int ind = r*j;
		wlr = wl.re[ind];
		wli = wl.im[ind];
		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];



			tau0r = br + cr;
			tau0i = bi + ci;

			tau1r = sgn * 0.86602540378 * (br - cr);
			tau1i = sgn * 0.86602540378 * (bi - ci);

			x.re[index] = ar + tau0r;
			x.im[index] = ai + tau0i;

			tau2r = ar - tau0r * 0.5000000000;
			tau2i = ai - tau0i * 0.5000000000;

			x.re[index1] = tau2r + tau1i;
			x.im[index1] = tau2i - tau1r;

			x.re[index2] = tau2r - tau1i;
			x.im[index2] = tau2i + tau1r;


		}
	}

}

template <typename T>
void inline radix4_mixed_dit_inplace(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 4;
	int r = n / L;

	T wlr,wli,wl2r,wl2i,wl3r,wl3i;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T ar,ai,br,bi,cr,ci,dr,di;

	for (int j = 0; j < Ls; j++) {
		int ind = r*j;
		wlr = wl.re[ind];
		wli = wl.im[ind];
		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];

			dr = wl3r*x.re[index3] - wl3i*x.im[index3];
			di = wl3r*x.im[index3] + wl3i*x.re[index3];

			tau0r = ar + cr;
			tau0i = ai + ci;

			tau1r = ar - cr;
			tau1i = ai - ci;

			tau2r = br + dr;
			tau2i = bi + di;

			tau3r = sgn* (br - dr);
			tau3i = sgn* (bi - di);

			x.re[index] = tau0r + tau2r;
			x.im[index] = tau0i + tau2i;

			x.re[index1] = tau1r + tau3i;
			x.im[index1] = tau1i - tau3r;

			x.re[index2] = tau0r - tau2r;
			x.im[index2] = tau0i - tau2i;

			x.re[index3] = tau1r - tau3i;
			x.im[index3] = tau1i + tau3r;


		}


	}

}

template <typename T>
void inline radix5_mixed_dit_inplace(fft_data<T> &x,fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 5;
	int r = n / L;


	T c1 = 0.30901699437;
	T c2 = -0.80901699437;
	T s1 = 0.95105651629;
	T s2 = 0.58778525229;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;

		wl4r = wl2r*wl2r - wl2i*wl2i;
		wl4i = 2.0*wl2r*wl2i;
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];

			dr = wl3r*x.re[index3] - wl3i*x.im[index3];
			di = wl3r*x.im[index3] + wl3i*x.re[index3];

			er = wl4r*x.re[index4] - wl4i*x.im[index4];
			ei = wl4r*x.im[index4] + wl4i*x.re[index4];

			tau0r = br + er;
			tau0i = bi + ei;

			tau1r = cr + dr;
			tau1i = ci + di;

			tau2r = br - er;
			tau2i = bi - ei;

			tau3r = cr - dr;
			tau3i = ci - di;

			x.re[index] = ar + tau0r + tau1r;
			x.im[index] = ai + tau0i + tau1i;

			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;

			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

			x.re[index1] = ar + tau4r + tau5i;
			x.im[index1] = ai + tau4i - tau5r;

			x.re[index4] = ar + tau4r - tau5i;
			x.im[index4] = ai + tau4i + tau5r;

			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;

			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

			x.re[index2] = ar + tau4r + tau5i;
			x.im[index2] = ai + tau4i - tau5r;

			x.re[index3] = ar + tau4r - tau5i;
			x.im[index3] = ai + tau4i + tau5r;


		}

	}

}

template <typename T>
void inline radix2_mixed_dif_inplace(fft_data<T> &x,int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 2;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T taur,taui;
	T xlsr,xlsi;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index1 = k*L+j+Ls;
			int index2 = k*L + j;
			taur = x.re[index1];
			taui = x.im[index1];

			xlsr = x.re[index2] - taur;
			xlsi = x.im[index2] - taui;

			x.re[index1] = wlr * xlsr - wli * xlsi;
			x.im[index1] = wlr * xlsi + wli * xlsr;

			x.re[index2] = x.re[index2] + taur;
			x.im[index2] = x.im[index2] + taui;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;
	}

}

template <typename T>
void inline radix3_mixed_dif_inplace(fft_data<T> &x, int q,int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 3;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T tau0r,tau0i,tau1r,tau1i;
	T br,bi,xlsr,xlsi;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;

			tau0r = x.re[index1] + x.re[index2];
			tau0i = x.im[index1] + x.im[index2];

			tau1r = sgn * 0.86602540378 * (x.re[index1] - x.re[index2]);
			tau1i = sgn * 0.86602540378 * (x.im[index1] - x.im[index2]);


			xlsr = x.re[index] - tau0r * 0.5000000000;
			xlsi = x.im[index] - tau0i * 0.5000000000;

			x.re[index]+= tau0r;
			x.im[index]+= tau0i;

			br = xlsr + tau1i;
			bi = xlsi - tau1r;

			x.re[index1] = wlr*br - wli*bi;
			x.im[index1] = wlr*bi + wli*br;

			br = xlsr - tau1i;
			bi = xlsi + tau1r;

			x.re[index2] = wl2r*br - wl2i*bi;
			x.im[index2] = wl2r*bi + wl2i*br;

		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

	}

}

template <typename T>
void inline radix4_mixed_dif_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 4;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T xlsr,xlsi;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;

			tau0r = x.re[index] + x.re[index2];
			tau0i = x.im[index] + x.im[index2];

			tau1r = x.re[index] - x.re[index2];
			tau1i = x.im[index] - x.im[index2];

			tau2r = x.re[index1] + x.re[index3];
			tau2i = x.im[index1] + x.im[index3];

			tau3r = sgn * (x.re[index1] - x.re[index3]);
			tau3i = sgn * (x.im[index1] - x.im[index3]);

			x.re[index] = tau0r + tau2r;
			x.im[index] = tau0i + tau2i;

			xlsr = tau1r + tau3i;
			xlsi = tau1i - tau3r;

			x.re[index1] = wlr * xlsr - wli * xlsi;
			x.im[index1] = wlr * xlsi + wli * xlsr;

			xlsr = tau0r - tau2r;
			xlsi = tau0i - tau2i;

			x.re[index2] = wl2r * xlsr - wl2i * xlsi;
			x.im[index2] = wl2r * xlsi + wl2i * xlsr;

			xlsr = tau1r - tau3i;
			xlsi = tau1i + tau3r;

			x.re[index3] = wl3r * xlsr - wl3i * xlsi;
			x.im[index3] = wl3r * xlsi + wl3i * xlsr;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;
	}

}

template <typename T>
void inline radix5_mixed_dif_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 5;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T c1 = 0.30901699437;
	T c2 = -0.80901699437;
	T s1 = 0.95105651629;
	T s2 = 0.58778525229;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T br,bi,cr,ci,dr,di,er,ei;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;

			tau0r = x.re[index1] + x.re[index4];
			tau0i = x.im[index1] + x.im[index4];

			tau1r = x.re[index2] + x.re[index3];
			tau1i = x.im[index2] + x.im[index3];

			tau2r = x.re[index1] - x.re[index4];
			tau2i = x.im[index1] - x.im[index4];

			tau3r = x.re[index2] - x.re[index3];
			tau3i = x.im[index2] - x.im[index3];

			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;

			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

			br = x.re[index] + tau4r + tau5i;
			bi = x.im[index] + tau4i - tau5r;

			er = x.re[index] + tau4r - tau5i;
			ei = x.im[index] + tau4i + tau5r;

			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;

			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

			cr = x.re[index] + tau4r + tau5i;
			ci = x.im[index] + tau4i - tau5r;

			dr = x.re[index] + tau4r - tau5i;
			di = x.im[index] + tau4i + tau5r;

			x.re[index]+= tau0r + tau1r;
			x.im[index]+= tau0i + tau1i;

			x.re[index1] = wlr*br - wli*bi;
			x.im[index1] = wlr*bi + wli*br;

			x.re[index2] = wl2r*cr - wl2i*ci;
			x.im[index2] = wl2r*ci + wl2i*cr;

			x.re[index3] = wl3r*dr - wl3i*di;
			x.im[index3] = wl3r*di + wl3i*dr;

			x.re[index4] = wl4r*er - wl4i*ei;
			x.im[index4] = wl4r*ei + wl4i*er;


		}
		T temp = wlr;
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

template <typename T>
void inline radix2_mixed_dif_inplace(fft_data<T> &x,fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 2;
	int r = n / L;

	T wlr,wli;
	T taur,taui;
	T xlsr,xlsi;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];
		for (int k =0; k < r; k++) {
			int index1 = k*L+j+Ls;
			int index2 = k*L + j;
			taur = x.re[index1];
			taui = x.im[index1];

			xlsr = x.re[index2] - taur;
			xlsi = x.im[index2] - taui;

			x.re[index1] = wlr * xlsr - wli * xlsi;
			x.im[index1] = wlr * xlsi + wli * xlsr;

			x.re[index2] = x.re[index2] + taur;
			x.im[index2] = x.im[index2] + taui;


		}

	}

}


template <typename T>
void inline radix3_mixed_dif_inplace(fft_data<T> &x,fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 3;
	int r = n / L;


	T tau0r,tau0i,tau1r,tau1i;
	T br,bi,xlsr,xlsi;
	T wlr,wli,wl2r,wl2i;

	for (int j = 0; j < Ls; j++) {
		int ind = r*j;
		wlr = wl.re[ind];
		wli = wl.im[ind];
		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;

			tau0r = x.re[index1] + x.re[index2];
			tau0i = x.im[index1] + x.im[index2];

			tau1r = sgn * 0.86602540378 * (x.re[index1] - x.re[index2]);
			tau1i = sgn * 0.86602540378 * (x.im[index1] - x.im[index2]);


			xlsr = x.re[index] - tau0r * 0.5000000000;
			xlsi = x.im[index] - tau0i * 0.5000000000;

			x.re[index]+= tau0r;
			x.im[index]+= tau0i;

			br = xlsr + tau1i;
			bi = xlsi - tau1r;

			x.re[index1] = wlr*br - wli*bi;
			x.im[index1] = wlr*bi + wli*br;

			br = xlsr - tau1i;
			bi = xlsi + tau1r;

			x.re[index2] = wl2r*br - wl2i*bi;
			x.im[index2] = wl2r*bi + wl2i*br;

		}

	}

}

template <typename T>
void inline radix4_mixed_dif_inplace(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 4;
	int r = n / L;

	T wlr,wli,wl2r,wl2i,wl3r,wl3i;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T xlsr,xlsi;

	for (int j = 0; j < Ls; j++) {
		int ind = r*j;
		wlr = wl.re[ind];
		wli = wl.im[ind];

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;

			tau0r = x.re[index] + x.re[index2];
			tau0i = x.im[index] + x.im[index2];

			tau1r = x.re[index] - x.re[index2];
			tau1i = x.im[index] - x.im[index2];

			tau2r = x.re[index1] + x.re[index3];
			tau2i = x.im[index1] + x.im[index3];

			tau3r = sgn * (x.re[index1] - x.re[index3]);
			tau3i = sgn * (x.im[index1] - x.im[index3]);

			x.re[index] = tau0r + tau2r;
			x.im[index] = tau0i + tau2i;

			xlsr = tau1r + tau3i;
			xlsi = tau1i - tau3r;

			x.re[index1] = wlr * xlsr - wli * xlsi;
			x.im[index1] = wlr * xlsi + wli * xlsr;

			xlsr = tau0r - tau2r;
			xlsi = tau0i - tau2i;

			x.re[index2] = wl2r * xlsr - wl2i * xlsi;
			x.im[index2] = wl2r * xlsi + wl2i * xlsr;

			xlsr = tau1r - tau3i;
			xlsi = tau1i + tau3r;

			x.re[index3] = wl3r * xlsr - wl3i * xlsi;
			x.im[index3] = wl3r * xlsi + wl3i * xlsr;


		}

	}

}

template <typename T>
void inline radix5_mixed_dif_inplace(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 5;
	int r = n / L;

	T c1 = 0.30901699437;
	T c2 = -0.80901699437;
	T s1 = 0.95105651629;
	T s2 = 0.58778525229;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T br,bi,cr,ci,dr,di,er,ei;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;

		wl4r = wl2r*wl2r - wl2i*wl2i;
		wl4i = 2.0*wl2r*wl2i;

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;

			tau0r = x.re[index1] + x.re[index4];
			tau0i = x.im[index1] + x.im[index4];

			tau1r = x.re[index2] + x.re[index3];
			tau1i = x.im[index2] + x.im[index3];

			tau2r = x.re[index1] - x.re[index4];
			tau2i = x.im[index1] - x.im[index4];

			tau3r = x.re[index2] - x.re[index3];
			tau3i = x.im[index2] - x.im[index3];

			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;

			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

			br = x.re[index] + tau4r + tau5i;
			bi = x.im[index] + tau4i - tau5r;

			er = x.re[index] + tau4r - tau5i;
			ei = x.im[index] + tau4i + tau5r;

			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;

			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

			cr = x.re[index] + tau4r + tau5i;
			ci = x.im[index] + tau4i - tau5r;

			dr = x.re[index] + tau4r - tau5i;
			di = x.im[index] + tau4i + tau5r;

			x.re[index]+= tau0r + tau1r;
			x.im[index]+= tau0i + tau1i;

			x.re[index1] = wlr*br - wli*bi;
			x.im[index1] = wlr*bi + wli*br;

			x.re[index2] = wl2r*cr - wl2i*ci;
			x.im[index2] = wl2r*ci + wl2i*cr;

			x.re[index3] = wl3r*dr - wl3i*di;
			x.im[index3] = wl3r*di + wl3i*dr;

			x.re[index4] = wl4r*er - wl4i*ei;
			x.im[index4] = wl4r*ei + wl4i*er;


		}

	}

}

template <typename T>
void inline sh_mixed_radix2_dit(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 2;
	int r = n / L;
	int rs = n / Ls;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T taur,taui;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index1 = j*rs + k ;
			int index = index1 +r;
			taur = wlr*y.re[index] - wli*y.im[index];
			taui = wlr*y.im[index] + wli*y.re[index];

			x.re[(j+Ls)*r+k] = y.re[index1] - taur;
			x.im[(j+Ls)*r+k] = y.im[index1] - taui;

			x.re[j*r+k] = y.re[index1] + taur;
			x.im[j*r+k] = y.im[index1] + taui;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;
	}

}

template <typename T>
void inline sh_mixed_radix3_dit(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 3;
	int r = n / L;
	int rs = n / Ls;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i;
	T ar,ai,br,bi,cr,ci;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = j*rs+k;
			int index1 = index+r;
			int index2 = index1+r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];



			tau0r = br + cr;
			tau0i = bi + ci;

			tau1r = sgn * 0.86602540378 * (br - cr);
			tau1i = sgn * 0.86602540378 * (bi - ci);

			int indexo = j*r+k;
			int indexo1 = indexo+Ls*r;
			int indexo2 = indexo1+Ls*r;

			x.re[indexo] = ar + tau0r;
			x.im[indexo] = ai + tau0i;

			tau2r = ar - tau0r * 0.5000000000;
			tau2i = ai - tau0i * 0.5000000000;

			x.re[indexo1] = tau2r + tau1i;
			x.im[indexo1] = tau2i - tau1r;

			x.re[indexo2] = tau2r - tau1i;
			x.im[indexo2] = tau2i + tau1r;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

	}

}

template <typename T>
void inline sh_mixed_radix4_dit(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 4;
	int r = n / L;
	int rs = n / Ls;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T ar,ai,br,bi,cr,ci,dr,di;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = j*rs + k ;
			int index1 = index +r;
			int index2 = index1 +r;
			int index3 = index2 +r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];

			dr = wl3r*y.re[index3] - wl3i*y.im[index3];
			di = wl3r*y.im[index3] + wl3i*y.re[index3];

			tau0r = ar + cr;
			tau0i = ai + ci;

			tau1r = ar - cr;
			tau1i = ai - ci;

			tau2r = br + dr;
			tau2i = bi + di;

			tau3r = sgn* (br - dr);
			tau3i = sgn* (bi - di);

			int indexo = j*r+k;
			int indexo1 = indexo+Ls*r;
			int indexo2 = indexo1+Ls*r;
			int indexo3 = indexo2+Ls*r;

			x.re[indexo] = tau0r + tau2r;
			x.im[indexo] = tau0i + tau2i;

			x.re[indexo1] = tau1r + tau3i;
			x.im[indexo1] = tau1i - tau3r;

			x.re[indexo2] = tau0r - tau2r;
			x.im[indexo2] = tau0i - tau2i;

			x.re[indexo3] = tau1r - tau3i;
			x.im[indexo3] = tau1i + tau3r;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;
	}

}

template <typename T>
void inline sh_mixed_radix5_dit(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 5;
	int r = n / L;
	int rs = n / Ls;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T c1 = 0.30901699437;
	T c2 = -0.80901699437;
	T s1 = 0.95105651629;
	T s2 = 0.58778525229;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;

	fft_data<T> y = x;
	int lsr = Ls*r;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = j*rs+k;
			int index1 = index+r;
			int index2 = index1+r;
			int index3 = index2+r;
			int index4 = index3+r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];

			dr = wl3r*y.re[index3] - wl3i*y.im[index3];
			di = wl3r*y.im[index3] + wl3i*y.re[index3];

			er = wl4r*y.re[index4] - wl4i*y.im[index4];
			ei = wl4r*y.im[index4] + wl4i*y.re[index4];

			tau0r = br + er;
			tau0i = bi + ei;

			tau1r = cr + dr;
			tau1i = ci + di;

			tau2r = br - er;
			tau2i = bi - ei;

			tau3r = cr - dr;
			tau3i = ci - di;

			int indexo = j*r+k;
			int indexo1 = indexo+lsr;
			int indexo2 = indexo1+lsr;
			int indexo3 = indexo2+lsr;
			int indexo4 = indexo3+lsr;


			x.re[indexo] = ar + tau0r + tau1r;
			x.im[indexo] = ai + tau0i + tau1i;

			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;

			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

			x.re[indexo1] = ar + tau4r + tau5i;
			x.im[indexo1] = ai + tau4i - tau5r;

			x.re[indexo4] = ar + tau4r - tau5i;
			x.im[indexo4] = ai + tau4i + tau5r;

			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;

			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

			x.re[indexo2] = ar + tau4r + tau5i;
			x.im[indexo2] = ai + tau4i - tau5r;

			x.re[indexo3] = ar + tau4r - tau5i;
			x.im[indexo3] = ai + tau4i + tau5r;


		}
		T temp = wlr;
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

template <typename T>
void inline sh_mixed_radix2_dit(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 2;
	int r = n / L;
	int rs = n / Ls;

	T wlr,wli;

	T taur,taui;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];
		for (int k =0; k < r; k++) {
			int index1 = j*rs + k ;
			int index = index1 +r;
			taur = wlr*y.re[index] - wli*y.im[index];
			taui = wlr*y.im[index] + wli*y.re[index];

			x.re[(j+Ls)*r+k] = y.re[index1] - taur;
			x.im[(j+Ls)*r+k] = y.im[index1] - taui;

			x.re[j*r+k] = y.re[index1] + taur;
			x.im[j*r+k] = y.im[index1] + taui;


		}

	}

}

template <typename T>
void inline sh_mixed_radix3_dit(fft_data<T> &x,fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 3;
	int r = n / L;
	int rs = n / Ls;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i;
	T ar,ai,br,bi,cr,ci;
	T wlr,wli,wl2r,wl2i;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];
		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;
		for (int k =0; k < r; k++) {
			int index = j*rs+k;
			int index1 = index+r;
			int index2 = index1+r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];

			tau0r = br + cr;
			tau0i = bi + ci;

			tau1r = sgn * 0.86602540378 * (br - cr);
			tau1i = sgn * 0.86602540378 * (bi - ci);

			int indexo = j*r+k;
			int indexo1 = indexo+Ls*r;
			int indexo2 = indexo1+Ls*r;

			x.re[indexo] = ar + tau0r;
			x.im[indexo] = ai + tau0i;

			tau2r = ar - tau0r * 0.5000000000;
			tau2i = ai - tau0i * 0.5000000000;

			x.re[indexo1] = tau2r + tau1i;
			x.im[indexo1] = tau2i - tau1r;

			x.re[indexo2] = tau2r - tau1i;
			x.im[indexo2] = tau2i + tau1r;


		}

	}

}

template <typename T>
void inline sh_mixed_radix4_dit(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 4;
	int r = n / L;
	int rs = n / Ls;

	T wlr,wli,wl2r,wl2i,wl3r,wl3i;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T ar,ai,br,bi,cr,ci,dr,di;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		int ind = r*j;
		wlr = wl.re[ind];
		wli = wl.im[ind];

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;

		for (int k =0; k < r; k++) {
			int index = j*rs + k ;
			int index1 = index +r;
			int index2 = index1 +r;
			int index3 = index2 +r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];

			dr = wl3r*y.re[index3] - wl3i*y.im[index3];
			di = wl3r*y.im[index3] + wl3i*y.re[index3];

			tau0r = ar + cr;
			tau0i = ai + ci;

			tau1r = ar - cr;
			tau1i = ai - ci;

			tau2r = br + dr;
			tau2i = bi + di;

			tau3r = sgn* (br - dr);
			tau3i = sgn* (bi - di);

			int indexo = j*r+k;
			int indexo1 = indexo+Ls*r;
			int indexo2 = indexo1+Ls*r;
			int indexo3 = indexo2+Ls*r;

			x.re[indexo] = tau0r + tau2r;
			x.im[indexo] = tau0i + tau2i;

			x.re[indexo1] = tau1r + tau3i;
			x.im[indexo1] = tau1i - tau3r;

			x.re[indexo2] = tau0r - tau2r;
			x.im[indexo2] = tau0i - tau2i;

			x.re[indexo3] = tau1r - tau3i;
			x.im[indexo3] = tau1i + tau3r;


		}

	}

}

template <typename T>
void inline sh_mixed_radix5_dit(fft_data<T> &x,fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 5;
	int r = n / L;
	int rs = n / Ls;

	T c1 = 0.30901699437;
	T c2 = -0.80901699437;
	T s1 = 0.95105651629;
	T s2 = 0.58778525229;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;

	fft_data<T> y = x;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i;
	int lsr = Ls*r;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;

		wl4r = wl2r*wl2r - wl2i*wl2i;
		wl4i = 2.0*wl2r*wl2i;
		for (int k =0; k < r; k++) {
			int index = j*rs+k;
			int index1 = index+r;
			int index2 = index1+r;
			int index3 = index2+r;
			int index4 = index3+r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];

			dr = wl3r*y.re[index3] - wl3i*y.im[index3];
			di = wl3r*y.im[index3] + wl3i*y.re[index3];

			er = wl4r*y.re[index4] - wl4i*y.im[index4];
			ei = wl4r*y.im[index4] + wl4i*y.re[index4];

			tau0r = br + er;
			tau0i = bi + ei;

			tau1r = cr + dr;
			tau1i = ci + di;

			tau2r = br - er;
			tau2i = bi - ei;

			tau3r = cr - dr;
			tau3i = ci - di;


			int indexo = j*r+k;
			int indexo1 = indexo+lsr;
			int indexo2 = indexo1+lsr;
			int indexo3 = indexo2+lsr;
			int indexo4 = indexo3+lsr;


			x.re[indexo] = ar + tau0r + tau1r;
			x.im[indexo] = ai + tau0i + tau1i;

			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;

			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

			x.re[indexo1] = ar + tau4r + tau5i;
			x.im[indexo1] = ai + tau4i - tau5r;

			x.re[indexo4] = ar + tau4r - tau5i;
			x.im[indexo4] = ai + tau4i + tau5r;

			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;

			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

			x.re[indexo2] = ar + tau4r + tau5i;
			x.im[indexo2] = ai + tau4i - tau5r;

			x.re[indexo3] = ar + tau4r - tau5i;
			x.im[indexo3] = ai + tau4i + tau5r;


		}

	}

}

template <typename T>
void inline sh_mixed_radix2_dif(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 2;
	int r = n / L;
	//int rs = n / Ls;
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T taur,taui;
	T xlsr,xlsi;

	fft_data<T> y = x;

	for (int j =0; j < Ls; j++) {
		for (int k = 0; k < r; k++) {
			int indexo2 = k*Ls + j ;
			int indexo1 = (k+r) * Ls + j;
			int index1 = k*L+j+Ls;
			int index2 = k*L + j;
			taur = y.re[index1];
			taui = y.im[index1];

			xlsr = y.re[index2] - taur;
			xlsi = y.im[index2] - taui;

			x.re[indexo1] = wlr * xlsr - wli * xlsi;
			x.im[indexo1] = wlr * xlsi + wli * xlsr;

			x.re[indexo2] = y.re[index2] + taur;
			x.im[indexo2] = y.im[index2] + taui;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;
	}

}

template <typename T>
void inline sh_mixed_radix3_dif(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 3;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T tau0r,tau0i,tau1r,tau1i;
	T br,bi,xlsr,xlsi;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;

			tau0r = y.re[index1] + y.re[index2];
			tau0i = y.im[index1] + y.im[index2];

			tau1r = sgn * 0.86602540378 * (y.re[index1] - y.re[index2]);
			tau1i = sgn * 0.86602540378 * (y.im[index1] - y.im[index2]);


			xlsr = y.re[index] - tau0r * 0.5000000000;
			xlsi = y.im[index] - tau0i * 0.5000000000;

			int indexo = k*Ls+j;
			int indexo1 = indexo+Ls*r;
			int indexo2 = indexo1+Ls*r;

			x.re[indexo] = y.re[index] + tau0r;
			x.im[indexo] = y.im[index] + tau0i;

			br = xlsr + tau1i;
			bi = xlsi - tau1r;

			x.re[indexo1] = wlr*br - wli*bi;
			x.im[indexo1] = wlr*bi + wli*br;

			br = xlsr - tau1i;
			bi = xlsi + tau1r;

			x.re[indexo2] = wl2r*br - wl2i*bi;
			x.im[indexo2] = wl2r*bi + wl2i*br;

		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

	}

}

template <typename T>
void inline sh_mixed_radix4_dif(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 4;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T xlsr,xlsi;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;

			tau0r = y.re[index] + y.re[index2];
			tau0i = y.im[index] + y.im[index2];

			tau1r = y.re[index] - y.re[index2];
			tau1i = y.im[index] - y.im[index2];

			tau2r = y.re[index1] + y.re[index3];
			tau2i = y.im[index1] + y.im[index3];

			tau3r = sgn * (y.re[index1] - y.re[index3]);
			tau3i = sgn * (y.im[index1] - y.im[index3]);

			int indexo = k*Ls+j;
			int indexo1 = indexo+Ls*r;
			int indexo2 = indexo1+Ls*r;
			int indexo3 = indexo2+Ls*r;

			x.re[indexo] = tau0r + tau2r;
			x.im[indexo] = tau0i + tau2i;

			xlsr = tau1r + tau3i;
			xlsi = tau1i - tau3r;

			x.re[indexo1] = wlr * xlsr - wli * xlsi;
			x.im[indexo1] = wlr * xlsi + wli * xlsr;

			xlsr = tau0r - tau2r;
			xlsi = tau0i - tau2i;

			x.re[indexo2] = wl2r * xlsr - wl2i * xlsi;
			x.im[indexo2] = wl2r * xlsi + wl2i * xlsr;

			xlsr = tau1r - tau3i;
			xlsi = tau1i + tau3r;

			x.re[indexo3] = wl3r * xlsr - wl3i * xlsi;
			x.im[indexo3] = wl3r * xlsi + wl3i * xlsr;


		}
		T temp = wlr;
		wlr = C * wlr - S * wli;
		wli = S * temp + C * wli;

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;
	}

}

template <typename T>
void inline sh_mixed_radix5_dif(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 5;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T c1 = 0.30901699437;
	T c2 = -0.80901699437;
	T s1 = 0.95105651629;
	T s2 = 0.58778525229;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T br,bi,cr,ci,dr,di,er,ei;

	fft_data<T> y = x;
	int lsr = Ls*r;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;

			tau0r = y.re[index1] + y.re[index4];
			tau0i = y.im[index1] + y.im[index4];

			tau1r = y.re[index2] + y.re[index3];
			tau1i = y.im[index2] + y.im[index3];

			tau2r = y.re[index1] - y.re[index4];
			tau2i = y.im[index1] - y.im[index4];

			tau3r = y.re[index2] - y.re[index3];
			tau3i = y.im[index2] - y.im[index3];

			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;

			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

			br = y.re[index] + tau4r + tau5i;
			bi = y.im[index] + tau4i - tau5r;

			er = y.re[index] + tau4r - tau5i;
			ei = y.im[index] + tau4i + tau5r;

			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;

			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

			cr = y.re[index] + tau4r + tau5i;
			ci = y.im[index] + tau4i - tau5r;

			dr = y.re[index] + tau4r - tau5i;
			di = y.im[index] + tau4i + tau5r;

			int indexo = k*Ls+j;
			int indexo1 = indexo+lsr;
			int indexo2 = indexo1+lsr;
			int indexo3 = indexo2+lsr;
			int indexo4 = indexo3+lsr;

			x.re[indexo]= y.re[index] + tau0r + tau1r;
			x.im[indexo]= y.im[index] + tau0i + tau1i;

			x.re[indexo1] = wlr*br - wli*bi;
			x.im[indexo1] = wlr*bi + wli*br;

			x.re[indexo2] = wl2r*cr - wl2i*ci;
			x.im[indexo2] = wl2r*ci + wl2i*cr;

			x.re[indexo3] = wl3r*dr - wl3i*di;
			x.im[indexo3] = wl3r*di + wl3i*dr;

			x.re[indexo4] = wl4r*er - wl4i*ei;
			x.im[indexo4] = wl4r*ei + wl4i*er;


		}
		T temp = wlr;
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

template <typename T>
void inline sh_mixed_radix2_dif(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 2;
	int r = n / L;
	//int rs = n / Ls;
	T wlr,wli;

	T taur,taui;
	T xlsr,xlsi;

	fft_data<T> y = x;

	for (int j =0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];
		for (int k = 0; k < r; k++) {
			int indexo2 = k*Ls + j ;
			int indexo1 = (k+r) * Ls + j;
			int index1 = k*L+j+Ls;
			int index2 = k*L + j;
			taur = y.re[index1];
			taui = y.im[index1];

			xlsr = y.re[index2] - taur;
			xlsi = y.im[index2] - taui;

			x.re[indexo1] = wlr * xlsr - wli * xlsi;
			x.im[indexo1] = wlr * xlsi + wli * xlsr;

			x.re[indexo2] = y.re[index2] + taur;
			x.im[indexo2] = y.im[index2] + taui;


		}

	}

}

template <typename T>
void inline sh_mixed_radix3_dif(fft_data<T> &x,fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 3;
	int r = n / L;

	T tau0r,tau0i,tau1r,tau1i;
	T br,bi,xlsr,xlsi;
	T wlr,wli,wl2r,wl2i;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];
		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;

			tau0r = y.re[index1] + y.re[index2];
			tau0i = y.im[index1] + y.im[index2];

			tau1r = sgn * 0.86602540378 * (y.re[index1] - y.re[index2]);
			tau1i = sgn * 0.86602540378 * (y.im[index1] - y.im[index2]);


			xlsr = y.re[index] - tau0r * 0.5000000000;
			xlsi = y.im[index] - tau0i * 0.5000000000;

			int indexo = k*Ls+j;
			int indexo1 = indexo+Ls*r;
			int indexo2 = indexo1+Ls*r;

			x.re[indexo] = y.re[index] + tau0r;
			x.im[indexo] = y.im[index] + tau0i;

			br = xlsr + tau1i;
			bi = xlsi - tau1r;

			x.re[indexo1] = wlr*br - wli*bi;
			x.im[indexo1] = wlr*bi + wli*br;

			br = xlsr - tau1i;
			bi = xlsi + tau1r;

			x.re[indexo2] = wl2r*br - wl2i*bi;
			x.im[indexo2] = wl2r*bi + wl2i*br;

		}

	}

}

template <typename T>
void inline sh_mixed_radix4_dif(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 4;
	int r = n / L;

	T wlr,wli,wl2r,wl2i,wl3r,wl3i;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T xlsr,xlsi;

	fft_data<T> y = x;

	for (int j = 0; j < Ls; j++) {
		int ind = r*j;
		wlr = wl.re[ind];
		wli = wl.im[ind];

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;

			tau0r = y.re[index] + y.re[index2];
			tau0i = y.im[index] + y.im[index2];

			tau1r = y.re[index] - y.re[index2];
			tau1i = y.im[index] - y.im[index2];

			tau2r = y.re[index1] + y.re[index3];
			tau2i = y.im[index1] + y.im[index3];

			tau3r = sgn * (y.re[index1] - y.re[index3]);
			tau3i = sgn * (y.im[index1] - y.im[index3]);

			int indexo = k*Ls+j;
			int indexo1 = indexo+Ls*r;
			int indexo2 = indexo1+Ls*r;
			int indexo3 = indexo2+Ls*r;

			x.re[indexo] = tau0r + tau2r;
			x.im[indexo] = tau0i + tau2i;

			xlsr = tau1r + tau3i;
			xlsi = tau1i - tau3r;

			x.re[indexo1] = wlr * xlsr - wli * xlsi;
			x.im[indexo1] = wlr * xlsi + wli * xlsr;

			xlsr = tau0r - tau2r;
			xlsi = tau0i - tau2i;

			x.re[indexo2] = wl2r * xlsr - wl2i * xlsi;
			x.im[indexo2] = wl2r * xlsi + wl2i * xlsr;

			xlsr = tau1r - tau3i;
			xlsi = tau1i + tau3r;

			x.re[indexo3] = wl3r * xlsr - wl3i * xlsi;
			x.im[indexo3] = wl3r * xlsi + wl3i * xlsr;


		}

	}

}

template <typename T>
void inline sh_mixed_radix5_dif(fft_data<T> &x,fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 5;
	int r = n / L;

	T c1 = 0.30901699437;
	T c2 = -0.80901699437;
	T s1 = 0.95105651629;
	T s2 = 0.58778525229;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T br,bi,cr,ci,dr,di,er,ei;

	fft_data<T> y = x;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i;
	int lsr = Ls*r;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

		wl2r = wlr*wlr - wli*wli;
		wl2i = 2.0*wlr*wli;

		wl3r = wl2r*wlr - wli*wl2i;
		wl3i= wl2r*wli + wl2i*wlr;

		wl4r = wl2r*wl2r - wl2i*wl2i;
		wl4i = 2.0*wl2r*wl2i;
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;

			tau0r = y.re[index1] + y.re[index4];
			tau0i = y.im[index1] + y.im[index4];

			tau1r = y.re[index2] + y.re[index3];
			tau1i = y.im[index2] + y.im[index3];

			tau2r = y.re[index1] - y.re[index4];
			tau2i = y.im[index1] - y.im[index4];

			tau3r = y.re[index2] - y.re[index3];
			tau3i = y.im[index2] - y.im[index3];

			tau4r = c1 * tau0r + c2 * tau1r;
			tau4i = c1 * tau0i + c2 * tau1i;

			tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
			tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

			br = y.re[index] + tau4r + tau5i;
			bi = y.im[index] + tau4i - tau5r;

			er = y.re[index] + tau4r - tau5i;
			ei = y.im[index] + tau4i + tau5r;

			tau4r = c2 * tau0r + c1 * tau1r;
			tau4i = c2 * tau0i + c1 * tau1i;

			tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
			tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

			cr = y.re[index] + tau4r + tau5i;
			ci = y.im[index] + tau4i - tau5r;

			dr = y.re[index] + tau4r - tau5i;
			di = y.im[index] + tau4i + tau5r;

			int indexo = k*Ls+j;
			int indexo1 = indexo+lsr;
			int indexo2 = indexo1+lsr;
			int indexo3 = indexo2+lsr;
			int indexo4 = indexo3+lsr;

			x.re[indexo]= y.re[index] + tau0r + tau1r;
			x.im[indexo]= y.im[index] + tau0i + tau1i;

			x.re[indexo1] = wlr*br - wli*bi;
			x.im[indexo1] = wlr*bi + wli*br;

			x.re[indexo2] = wl2r*cr - wl2i*ci;
			x.im[indexo2] = wl2r*ci + wl2i*cr;

			x.re[indexo3] = wl3r*dr - wl3i*di;
			x.im[indexo3] = wl3r*di + wl3i*dr;

			x.re[indexo4] = wl4r*er - wl4i*ei;
			x.im[indexo4] = wl4r*ei + wl4i*er;


		}

	}

}

template <typename T>
void inline radix7_mixed_dit_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 7;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T wl5r = (T) 1.0;
	T wl5i = (T) 0.0;

	T wl6r = (T) 1.0;
	T wl6i = (T) 0.0;

	T c1 = 0.62348980185;
	T c2 = -0.22252093395;
	T c3 = -0.9009688679;
	T s1 = 0.78183148246;
	T s2 = 0.97492791218;
	T s3 = 0.43388373911;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];

			dr = wl3r*x.re[index3] - wl3i*x.im[index3];
			di = wl3r*x.im[index3] + wl3i*x.re[index3];

			er = wl4r*x.re[index4] - wl4i*x.im[index4];
			ei = wl4r*x.im[index4] + wl4i*x.re[index4];

			fr = wl5r*x.re[index5] - wl5i*x.im[index5];
			fi = wl5r*x.im[index5] + wl5i*x.re[index5];

			gr = wl6r*x.re[index6] - wl6i*x.im[index6];
			gi = wl6r*x.im[index6] + wl6i*x.re[index6];

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

			x.re[index] = ar + tau0r + tau1r + tau2r;
			x.im[index] = ai + tau0i + tau1i + tau2i;

			tau6r = ar + c1 * tau0r + c2 * tau1r + c3 * tau2r;
			tau6i = ai + c1 * tau0i + c2 * tau1i + c3 * tau2i;

			tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
			tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

			x.re[index1] = tau6r - tau7i;
			x.im[index1] = tau6i + tau7r;

			x.re[index6] = tau6r + tau7i;
			x.im[index6] = tau6i - tau7r;

			tau6r = ar + c2 * tau0r + c3 * tau1r + c1 * tau2r;
			tau6i = ai + c2 * tau0i + c3 * tau1i + c1 * tau2i;

			tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
			tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

			x.re[index2] = tau6r - tau7i;
			x.im[index2] = tau6i + tau7r;

			x.re[index5] = tau6r + tau7i;
			x.im[index5] = tau6i - tau7r;

			tau6r = ar + c3 * tau0r + c1 * tau1r + c2 * tau2r;
			tau6i = ai + c3 * tau0i + c1 * tau1i + c2 * tau2i;

			tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
			tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

			x.re[index3] = tau6r - tau7i;
			x.im[index3] = tau6i + tau7r;

			x.re[index4] = tau6r + tau7i;
			x.im[index4] = tau6i - tau7r;


		}
		T temp = wlr;
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

template <typename T>
void inline radix7_mixed_dit_inplace(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 7;
	int r = n / L;

	T c1 = 0.62348980185;
	T c2 = -0.22252093395;
	T c3 = -0.9009688679;
	T s1 = 0.78183148246;
	T s2 = 0.97492791218;
	T s3 = 0.43388373911;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i,wl5r,wl5i,wl6r,wl6i;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

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

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];

			dr = wl3r*x.re[index3] - wl3i*x.im[index3];
			di = wl3r*x.im[index3] + wl3i*x.re[index3];

			er = wl4r*x.re[index4] - wl4i*x.im[index4];
			ei = wl4r*x.im[index4] + wl4i*x.re[index4];

			fr = wl5r*x.re[index5] - wl5i*x.im[index5];
			fi = wl5r*x.im[index5] + wl5i*x.re[index5];

			gr = wl6r*x.re[index6] - wl6i*x.im[index6];
			gi = wl6r*x.im[index6] + wl6i*x.re[index6];

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

			x.re[index] = ar + tau0r + tau1r + tau2r;
			x.im[index] = ai + tau0i + tau1i + tau2i;

			tau6r = ar + c1 * tau0r + c2 * tau1r + c3 * tau2r;
			tau6i = ai + c1 * tau0i + c2 * tau1i + c3 * tau2i;

			tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
			tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

			x.re[index1] = tau6r - tau7i;
			x.im[index1] = tau6i + tau7r;

			x.re[index6] = tau6r + tau7i;
			x.im[index6] = tau6i - tau7r;

			tau6r = ar + c2 * tau0r + c3 * tau1r + c1 * tau2r;
			tau6i = ai + c2 * tau0i + c3 * tau1i + c1 * tau2i;

			tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
			tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

			x.re[index2] = tau6r - tau7i;
			x.im[index2] = tau6i + tau7r;

			x.re[index5] = tau6r + tau7i;
			x.im[index5] = tau6i - tau7r;

			tau6r = ar + c3 * tau0r + c1 * tau1r + c2 * tau2r;
			tau6i = ai + c3 * tau0i + c1 * tau1i + c2 * tau2i;

			tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
			tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

			x.re[index3] = tau6r - tau7i;
			x.im[index3] = tau6i + tau7r;

			x.re[index4] = tau6r + tau7i;
			x.im[index4] = tau6i - tau7r;


		}

	}

}

template <typename T>
void inline radix7_mixed_dif_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 7;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T wl5r = (T) 1.0;
	T wl5i = (T) 0.0;

	T wl6r = (T) 1.0;
	T wl6i = (T) 0.0;

	T c1 = 0.62348980185;
	T c2 = -0.22252093395;
	T c3 = -0.9009688679;
	T s1 = 0.78183148246;
	T s2 = 0.97492791218;
	T s3 = 0.43388373911;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;

			tau0r = x.re[index1] + x.re[index6];
			tau0i = x.im[index1] + x.im[index6];

			tau1r = x.re[index2] + x.re[index5];
			tau1i = x.im[index2] + x.im[index5];

			tau2r = x.re[index3] + x.re[index4];
			tau2i = x.im[index3] + x.im[index4];

			tau3r = x.re[index1] - x.re[index6];
			tau3i = x.im[index1] - x.im[index6];

			tau4r = x.re[index2] - x.re[index5];
			tau4i = x.im[index2] - x.im[index5];

			tau5r = x.re[index3] - x.re[index4];
			tau5i = x.im[index3] - x.im[index4];

			tau6r = c1 * tau0r + c2 * tau1r + c3 * tau2r;
			tau6i = c1 * tau0i + c2 * tau1i + c3 * tau2i;

			tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
			tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

			br = x.re[index] + tau6r - tau7i;
			bi = x.im[index] + tau6i + tau7r;

			gr = x.re[index] + tau6r + tau7i;
			gi = x.im[index] + tau6i - tau7r;

			tau6r = c2 * tau0r + c3 * tau1r + c1 * tau2r;
			tau6i = c2 * tau0i + c3 * tau1i + c1 * tau2i;

			tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
			tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

			cr = x.re[index] + tau6r - tau7i;
			ci = x.im[index] + tau6i + tau7r;

			fr = x.re[index] + tau6r + tau7i;
			fi = x.im[index] + tau6i - tau7r;

			tau6r = c3 * tau0r + c1 * tau1r + c2 * tau2r;
			tau6i = c3 * tau0i + c1 * tau1i + c2 * tau2i;

			tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
			tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

			dr = x.re[index] + tau6r - tau7i;
			di = x.im[index] + tau6i + tau7r;

			er = x.re[index] + tau6r + tau7i;
			ei = x.im[index] + tau6i - tau7r;

			x.re[index] += tau0r + tau1r + tau2r;
			x.im[index] += tau0i + tau1i + tau2i;

			x.re[index1] = wlr*br - wli*bi;
			x.im[index1] = wlr*bi + wli*br;

			x.re[index2] = wl2r*cr - wl2i*ci;
			x.im[index2] = wl2r*ci + wl2i*cr;

			x.re[index3] = wl3r*dr - wl3i*di;
			x.im[index3] = wl3r*di + wl3i*dr;

			x.re[index4] = wl4r*er - wl4i*ei;
			x.im[index4] = wl4r*ei + wl4i*er;

			x.re[index5] = wl5r*fr - wl5i*fi;
			x.im[index5] = wl5r*fi + wl5i*fr;

			x.re[index6] = wl6r*gr - wl6i*gi;
			x.im[index6] = wl6r*gi + wl6i*gr;


		}
		T temp = wlr;
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

template <typename T>
void inline radix7_mixed_dif_inplace(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 7;
	int r = n / L;

	T c1 = 0.62348980185;
	T c2 = -0.22252093395;
	T c3 = -0.9009688679;
	T s1 = 0.78183148246;
	T s2 = 0.97492791218;
	T s3 = 0.43388373911;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i,wl5r,wl5i,wl6r,wl6i;

	for (int j = 0; j < Ls; j++) {
		int ind = r*j;
		wlr = wl.re[ind];
		wli = wl.im[ind];

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

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;

			tau0r = x.re[index1] + x.re[index6];
			tau0i = x.im[index1] + x.im[index6];

			tau1r = x.re[index2] + x.re[index5];
			tau1i = x.im[index2] + x.im[index5];

			tau2r = x.re[index3] + x.re[index4];
			tau2i = x.im[index3] + x.im[index4];

			tau3r = x.re[index1] - x.re[index6];
			tau3i = x.im[index1] - x.im[index6];

			tau4r = x.re[index2] - x.re[index5];
			tau4i = x.im[index2] - x.im[index5];

			tau5r = x.re[index3] - x.re[index4];
			tau5i = x.im[index3] - x.im[index4];

			tau6r = c1 * tau0r + c2 * tau1r + c3 * tau2r;
			tau6i = c1 * tau0i + c2 * tau1i + c3 * tau2i;

			tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
			tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

			br = x.re[index] + tau6r - tau7i;
			bi = x.im[index] + tau6i + tau7r;

			gr = x.re[index] + tau6r + tau7i;
			gi = x.im[index] + tau6i - tau7r;

			tau6r = c2 * tau0r + c3 * tau1r + c1 * tau2r;
			tau6i = c2 * tau0i + c3 * tau1i + c1 * tau2i;

			tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
			tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

			cr = x.re[index] + tau6r - tau7i;
			ci = x.im[index] + tau6i + tau7r;

			fr = x.re[index] + tau6r + tau7i;
			fi = x.im[index] + tau6i - tau7r;

			tau6r = c3 * tau0r + c1 * tau1r + c2 * tau2r;
			tau6i = c3 * tau0i + c1 * tau1i + c2 * tau2i;

			tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
			tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

			dr = x.re[index] + tau6r - tau7i;
			di = x.im[index] + tau6i + tau7r;

			er = x.re[index] + tau6r + tau7i;
			ei = x.im[index] + tau6i - tau7r;

			x.re[index] += tau0r + tau1r + tau2r;
			x.im[index] += tau0i + tau1i + tau2i;

			x.re[index1] = wlr*br - wli*bi;
			x.im[index1] = wlr*bi + wli*br;

			x.re[index2] = wl2r*cr - wl2i*ci;
			x.im[index2] = wl2r*ci + wl2i*cr;

			x.re[index3] = wl3r*dr - wl3i*di;
			x.im[index3] = wl3r*di + wl3i*dr;

			x.re[index4] = wl4r*er - wl4i*ei;
			x.im[index4] = wl4r*ei + wl4i*er;

			x.re[index5] = wl5r*fr - wl5i*fi;
			x.im[index5] = wl5r*fi + wl5i*fr;

			x.re[index6] = wl6r*gr - wl6i*gi;
			x.im[index6] = wl6r*gi + wl6i*gr;


		}

	}

}

template <typename T>
void inline sh_mixed_radix7_dit(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 7;
	int r = n / L;
	int rs = n / Ls;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T wl5r = (T) 1.0;
	T wl5i = (T) 0.0;

	T wl6r = (T) 1.0;
	T wl6i = (T) 0.0;

	T c1 = 0.62348980185;
	T c2 = -0.22252093395;
	T c3 = -0.9009688679;
	T s1 = 0.78183148246;
	T s2 = 0.97492791218;
	T s3 = 0.43388373911;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi;

	fft_data<T> y = x;
	int lsr = Ls *r;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = j*rs + k;
			int index1 = index+r;
			int index2 = index1+r;
			int index3 = index2+r;
			int index4 = index3+r;
			int index5 = index4+r;
			int index6 = index5+r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];

			dr = wl3r*y.re[index3] - wl3i*y.im[index3];
			di = wl3r*y.im[index3] + wl3i*y.re[index3];

			er = wl4r*y.re[index4] - wl4i*y.im[index4];
			ei = wl4r*y.im[index4] + wl4i*y.re[index4];

			fr = wl5r*y.re[index5] - wl5i*y.im[index5];
			fi = wl5r*y.im[index5] + wl5i*y.re[index5];

			gr = wl6r*y.re[index6] - wl6i*y.im[index6];
			gi = wl6r*y.im[index6] + wl6i*y.re[index6];

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

			int indexo = j*r + k;
			int indexo1 = indexo + lsr;
			int indexo2 = indexo1 + lsr;
			int indexo3 = indexo2 + lsr;
			int indexo4 = indexo3 + lsr;
			int indexo5 = indexo4 + lsr;
			int indexo6 = indexo5 + lsr;

			x.re[indexo] = ar + tau0r + tau1r + tau2r;
			x.im[indexo] = ai + tau0i + tau1i + tau2i;

			tau6r = ar + c1 * tau0r + c2 * tau1r + c3 * tau2r;
			tau6i = ai + c1 * tau0i + c2 * tau1i + c3 * tau2i;

			tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
			tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

			x.re[indexo1] = tau6r - tau7i;
			x.im[indexo1] = tau6i + tau7r;

			x.re[indexo6] = tau6r + tau7i;
			x.im[indexo6] = tau6i - tau7r;

			tau6r = ar + c2 * tau0r + c3 * tau1r + c1 * tau2r;
			tau6i = ai + c2 * tau0i + c3 * tau1i + c1 * tau2i;

			tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
			tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

			x.re[indexo2] = tau6r - tau7i;
			x.im[indexo2] = tau6i + tau7r;

			x.re[indexo5] = tau6r + tau7i;
			x.im[indexo5] = tau6i - tau7r;

			tau6r = ar + c3 * tau0r + c1 * tau1r + c2 * tau2r;
			tau6i = ai + c3 * tau0i + c1 * tau1i + c2 * tau2i;

			tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
			tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

			x.re[indexo3] = tau6r - tau7i;
			x.im[indexo3] = tau6i + tau7r;

			x.re[indexo4] = tau6r + tau7i;
			x.im[indexo4] = tau6i - tau7r;


		}
		T temp = wlr;
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

template <typename T>
void inline sh_mixed_radix7_dit(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 7;
	int r = n / L;
	int rs = n / Ls;

	T c1 = 0.62348980185;
	T c2 = -0.22252093395;
	T c3 = -0.9009688679;
	T s1 = 0.78183148246;
	T s2 = 0.97492791218;
	T s3 = 0.43388373911;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i,wl5r,wl5i,wl6r,wl6i;

	fft_data<T> y = x;
	int lsr = Ls *r;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

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


		for (int k =0; k < r; k++) {
			int index = j*rs + k;
			int index1 = index+r;
			int index2 = index1+r;
			int index3 = index2+r;
			int index4 = index3+r;
			int index5 = index4+r;
			int index6 = index5+r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];

			dr = wl3r*y.re[index3] - wl3i*y.im[index3];
			di = wl3r*y.im[index3] + wl3i*y.re[index3];

			er = wl4r*y.re[index4] - wl4i*y.im[index4];
			ei = wl4r*y.im[index4] + wl4i*y.re[index4];

			fr = wl5r*y.re[index5] - wl5i*y.im[index5];
			fi = wl5r*y.im[index5] + wl5i*y.re[index5];

			gr = wl6r*y.re[index6] - wl6i*y.im[index6];
			gi = wl6r*y.im[index6] + wl6i*y.re[index6];

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

			int indexo = j*r + k;
			int indexo1 = indexo + lsr;
			int indexo2 = indexo1 + lsr;
			int indexo3 = indexo2 + lsr;
			int indexo4 = indexo3 + lsr;
			int indexo5 = indexo4 + lsr;
			int indexo6 = indexo5 + lsr;

			x.re[indexo] = ar + tau0r + tau1r + tau2r;
			x.im[indexo] = ai + tau0i + tau1i + tau2i;

			tau6r = ar + c1 * tau0r + c2 * tau1r + c3 * tau2r;
			tau6i = ai + c1 * tau0i + c2 * tau1i + c3 * tau2i;

			tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
			tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

			x.re[indexo1] = tau6r - tau7i;
			x.im[indexo1] = tau6i + tau7r;

			x.re[indexo6] = tau6r + tau7i;
			x.im[indexo6] = tau6i - tau7r;

			tau6r = ar + c2 * tau0r + c3 * tau1r + c1 * tau2r;
			tau6i = ai + c2 * tau0i + c3 * tau1i + c1 * tau2i;

			tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
			tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

			x.re[indexo2] = tau6r - tau7i;
			x.im[indexo2] = tau6i + tau7r;

			x.re[indexo5] = tau6r + tau7i;
			x.im[indexo5] = tau6i - tau7r;

			tau6r = ar + c3 * tau0r + c1 * tau1r + c2 * tau2r;
			tau6i = ai + c3 * tau0i + c1 * tau1i + c2 * tau2i;

			tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
			tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

			x.re[indexo3] = tau6r - tau7i;
			x.im[indexo3] = tau6i + tau7r;

			x.re[indexo4] = tau6r + tau7i;
			x.im[indexo4] = tau6i - tau7r;


		}

	}

}

template <typename T>
void inline sh_mixed_radix7_dif(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 7;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T wl5r = (T) 1.0;
	T wl5i = (T) 0.0;

	T wl6r = (T) 1.0;
	T wl6i = (T) 0.0;

	T c1 = 0.62348980185;
	T c2 = -0.22252093395;
	T c3 = -0.9009688679;
	T s1 = 0.78183148246;
	T s2 = 0.97492791218;
	T s3 = 0.43388373911;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi;

	fft_data<T> y = x;
	int lsr = Ls*r;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;

			tau0r = y.re[index1] + y.re[index6];
			tau0i = y.im[index1] + y.im[index6];

			tau1r = y.re[index2] + y.re[index5];
			tau1i = y.im[index2] + y.im[index5];

			tau2r = y.re[index3] + y.re[index4];
			tau2i = y.im[index3] + y.im[index4];

			tau3r = y.re[index1] - y.re[index6];
			tau3i = y.im[index1] - y.im[index6];

			tau4r = y.re[index2] - y.re[index5];
			tau4i = y.im[index2] - y.im[index5];

			tau5r = y.re[index3] - y.re[index4];
			tau5i = y.im[index3] - y.im[index4];

			tau6r = c1 * tau0r + c2 * tau1r + c3 * tau2r;
			tau6i = c1 * tau0i + c2 * tau1i + c3 * tau2i;

			tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
			tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

			br = y.re[index] + tau6r - tau7i;
			bi = y.im[index] + tau6i + tau7r;

			gr = y.re[index] + tau6r + tau7i;
			gi = y.im[index] + tau6i - tau7r;

			tau6r = c2 * tau0r + c3 * tau1r + c1 * tau2r;
			tau6i = c2 * tau0i + c3 * tau1i + c1 * tau2i;

			tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
			tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

			cr = y.re[index] + tau6r - tau7i;
			ci = y.im[index] + tau6i + tau7r;

			fr = y.re[index] + tau6r + tau7i;
			fi = y.im[index] + tau6i - tau7r;

			tau6r = c3 * tau0r + c1 * tau1r + c2 * tau2r;
			tau6i = c3 * tau0i + c1 * tau1i + c2 * tau2i;

			tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
			tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

			dr = y.re[index] + tau6r - tau7i;
			di = y.im[index] + tau6i + tau7r;

			er = y.re[index] + tau6r + tau7i;
			ei = y.im[index] + tau6i - tau7r;

			int indexo = k*Ls+j;
			int indexo1 = indexo+lsr;
			int indexo2 = indexo1+lsr;
			int indexo3 = indexo2+lsr;
			int indexo4 = indexo3+lsr;
			int indexo5 = indexo4+lsr;
			int indexo6 = indexo5+lsr;

			x.re[indexo] = y.re[index] + tau0r + tau1r + tau2r;
			x.im[indexo] = y.im[index] + tau0i + tau1i + tau2i;

			x.re[indexo1] = wlr*br - wli*bi;
			x.im[indexo1] = wlr*bi + wli*br;

			x.re[indexo2] = wl2r*cr - wl2i*ci;
			x.im[indexo2] = wl2r*ci + wl2i*cr;

			x.re[indexo3] = wl3r*dr - wl3i*di;
			x.im[indexo3] = wl3r*di + wl3i*dr;

			x.re[indexo4] = wl4r*er - wl4i*ei;
			x.im[indexo4] = wl4r*ei + wl4i*er;

			x.re[indexo5] = wl5r*fr - wl5i*fi;
			x.im[indexo5] = wl5r*fi + wl5i*fr;

			x.re[indexo6] = wl6r*gr - wl6i*gi;
			x.im[indexo6] = wl6r*gi + wl6i*gr;


		}
		T temp = wlr;
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

template <typename T>
void inline sh_mixed_radix7_dif(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 7;
	int r = n / L;

	T c1 = 0.62348980185;
	T c2 = -0.22252093395;
	T c3 = -0.9009688679;
	T s1 = 0.78183148246;
	T s2 = 0.97492791218;
	T s3 = 0.43388373911;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i,wl5r,wl5i,wl6r,wl6i;

	fft_data<T> y = x;
	int lsr = Ls*r;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

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

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;

			tau0r = y.re[index1] + y.re[index6];
			tau0i = y.im[index1] + y.im[index6];

			tau1r = y.re[index2] + y.re[index5];
			tau1i = y.im[index2] + y.im[index5];

			tau2r = y.re[index3] + y.re[index4];
			tau2i = y.im[index3] + y.im[index4];

			tau3r = y.re[index1] - y.re[index6];
			tau3i = y.im[index1] - y.im[index6];

			tau4r = y.re[index2] - y.re[index5];
			tau4i = y.im[index2] - y.im[index5];

			tau5r = y.re[index3] - y.re[index4];
			tau5i = y.im[index3] - y.im[index4];

			tau6r = c1 * tau0r + c2 * tau1r + c3 * tau2r;
			tau6i = c1 * tau0i + c2 * tau1i + c3 * tau2i;

			tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
			tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

			br = y.re[index] + tau6r - tau7i;
			bi = y.im[index] + tau6i + tau7r;

			gr = y.re[index] + tau6r + tau7i;
			gi = y.im[index] + tau6i - tau7r;

			tau6r = c2 * tau0r + c3 * tau1r + c1 * tau2r;
			tau6i = c2 * tau0i + c3 * tau1i + c1 * tau2i;

			tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
			tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

			cr = y.re[index] + tau6r - tau7i;
			ci = y.im[index] + tau6i + tau7r;

			fr = y.re[index] + tau6r + tau7i;
			fi = y.im[index] + tau6i - tau7r;

			tau6r = c3 * tau0r + c1 * tau1r + c2 * tau2r;
			tau6i = c3 * tau0i + c1 * tau1i + c2 * tau2i;

			tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
			tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

			dr = y.re[index] + tau6r - tau7i;
			di = y.im[index] + tau6i + tau7r;

			er = y.re[index] + tau6r + tau7i;
			ei = y.im[index] + tau6i - tau7r;

			int indexo = k*Ls+j;
			int indexo1 = indexo+lsr;
			int indexo2 = indexo1+lsr;
			int indexo3 = indexo2+lsr;
			int indexo4 = indexo3+lsr;
			int indexo5 = indexo4+lsr;
			int indexo6 = indexo5+lsr;

			x.re[indexo] = y.re[index] + tau0r + tau1r + tau2r;
			x.im[indexo] = y.im[index] + tau0i + tau1i + tau2i;

			x.re[indexo1] = wlr*br - wli*bi;
			x.im[indexo1] = wlr*bi + wli*br;

			x.re[indexo2] = wl2r*cr - wl2i*ci;
			x.im[indexo2] = wl2r*ci + wl2i*cr;

			x.re[indexo3] = wl3r*dr - wl3i*di;
			x.im[indexo3] = wl3r*di + wl3i*dr;

			x.re[indexo4] = wl4r*er - wl4i*ei;
			x.im[indexo4] = wl4r*ei + wl4i*er;

			x.re[indexo5] = wl5r*fr - wl5i*fi;
			x.im[indexo5] = wl5r*fi + wl5i*fr;

			x.re[indexo6] = wl6r*gr - wl6i*gi;
			x.im[indexo6] = wl6r*gi + wl6i*gr;


		}

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

template <typename T>
void inline radix8_mixed_dit_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 8;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T wl5r = (T) 1.0;
	T wl5i = (T) 0.0;

	T wl6r = (T) 1.0;
	T wl6i = (T) 0.0;

	T wl7r = (T) 1.0;
	T wl7i = (T) 0.0;

	T c1 = 0.70710678118654752440084436210485;
	//T c2 = -0.70710678118654752440084436210485;
	T s1 = 0.70710678118654752440084436210485;
	//T s2 = 0.70710678118654752440084436210485;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T tau8r,tau8i,tau9r,tau9i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi,hr,hi;
	T temp1r,temp1i,temp2r,temp2i;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;
			int index7 = index6+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];

			dr = wl3r*x.re[index3] - wl3i*x.im[index3];
			di = wl3r*x.im[index3] + wl3i*x.re[index3];

			er = wl4r*x.re[index4] - wl4i*x.im[index4];
			ei = wl4r*x.im[index4] + wl4i*x.re[index4];

			fr = wl5r*x.re[index5] - wl5i*x.im[index5];
			fi = wl5r*x.im[index5] + wl5i*x.re[index5];

			gr = wl6r*x.re[index6] - wl6i*x.im[index6];
			gi = wl6r*x.im[index6] + wl6i*x.re[index6];

			hr = wl7r*x.re[index7] - wl7i*x.im[index7];
			hi = wl7r*x.im[index7] + wl7i*x.re[index7];

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

			x.re[index] = tau0r + tau1r + tau2r + tau3r;
			x.im[index] = tau0i + tau1i + tau2i + tau3i;

			x.re[index4] = tau0r - tau1r - tau2r + tau3r;
			x.im[index4] = tau0i - tau1i - tau2i + tau3i;

			temp1r = tau1r - tau2r;
			temp1i = tau1i - tau2i;

			temp2r = tau5r + tau6r;
			temp2i = tau5i + tau6i;

			tau8r =  tau4r + c1 * temp1r;
			tau8i =  tau4i + c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r - tau7r);
			tau9i = sgn * ( -s1 * temp2i - tau7i);

			x.re[index1] = tau8r - tau9i;
			x.im[index1] = tau8i + tau9r;

			x.re[index7] = tau8r + tau9i;
			x.im[index7] = tau8i - tau9r;

			tau8r = tau0r - tau3r;
			tau8i = tau0i - tau3i;

			tau9r = sgn * ( -tau5r + tau6r);
			tau9i = sgn * ( -tau5i + tau6i);

			x.re[index2] = tau8r - tau9i;
			x.im[index2] = tau8i + tau9r;

			x.re[index6] = tau8r + tau9i;
			x.im[index6] = tau8i - tau9r;

			tau8r = tau4r - c1 * temp1r;
			tau8i = tau4i - c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r + tau7r);
			tau9i = sgn * ( -s1 * temp2i + tau7i);

			x.re[index3] = tau8r - tau9i;
			x.im[index3] = tau8i + tau9r;

			x.re[index5] = tau8r + tau9i;
			x.im[index5] = tau8i - tau9r;


		}
		T temp = wlr;
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

template <typename T>
void inline radix8_mixed_dit_inplace(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 8;
	int r = n / L;


	T c1 = 0.70710678118654752440084436210485;
	//T c2 = -0.70710678118654752440084436210485;
	T s1 = 0.70710678118654752440084436210485;
	//T s2 = 0.70710678118654752440084436210485;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T tau8r,tau8i,tau9r,tau9i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi,hr,hi;
	T temp1r,temp1i,temp2r,temp2i;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i,wl5r,wl5i,wl6r,wl6i,wl7r,wl7i;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

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

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;
			int index7 = index6+Ls;

			ar = x.re[index];
			ai = x.im[index];

			br = wlr*x.re[index1] - wli*x.im[index1];
			bi = wlr*x.im[index1] + wli*x.re[index1];

			cr = wl2r*x.re[index2] - wl2i*x.im[index2];
			ci = wl2r*x.im[index2] + wl2i*x.re[index2];

			dr = wl3r*x.re[index3] - wl3i*x.im[index3];
			di = wl3r*x.im[index3] + wl3i*x.re[index3];

			er = wl4r*x.re[index4] - wl4i*x.im[index4];
			ei = wl4r*x.im[index4] + wl4i*x.re[index4];

			fr = wl5r*x.re[index5] - wl5i*x.im[index5];
			fi = wl5r*x.im[index5] + wl5i*x.re[index5];

			gr = wl6r*x.re[index6] - wl6i*x.im[index6];
			gi = wl6r*x.im[index6] + wl6i*x.re[index6];

			hr = wl7r*x.re[index7] - wl7i*x.im[index7];
			hi = wl7r*x.im[index7] + wl7i*x.re[index7];

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

			x.re[index] = tau0r + tau1r + tau2r + tau3r;
			x.im[index] = tau0i + tau1i + tau2i + tau3i;

			x.re[index4] = tau0r - tau1r - tau2r + tau3r;
			x.im[index4] = tau0i - tau1i - tau2i + tau3i;

			temp1r = tau1r - tau2r;
			temp1i = tau1i - tau2i;

			temp2r = tau5r + tau6r;
			temp2i = tau5i + tau6i;

			tau8r =  tau4r + c1 * temp1r;
			tau8i =  tau4i + c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r - tau7r);
			tau9i = sgn * ( -s1 * temp2i - tau7i);

			x.re[index1] = tau8r - tau9i;
			x.im[index1] = tau8i + tau9r;

			x.re[index7] = tau8r + tau9i;
			x.im[index7] = tau8i - tau9r;

			tau8r = tau0r - tau3r;
			tau8i = tau0i - tau3i;

			tau9r = sgn * ( -tau5r + tau6r);
			tau9i = sgn * ( -tau5i + tau6i);

			x.re[index2] = tau8r - tau9i;
			x.im[index2] = tau8i + tau9r;

			x.re[index6] = tau8r + tau9i;
			x.im[index6] = tau8i - tau9r;

			tau8r = tau4r - c1 * temp1r;
			tau8i = tau4i - c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r + tau7r);
			tau9i = sgn * ( -s1 * temp2i + tau7i);

			x.re[index3] = tau8r - tau9i;
			x.im[index3] = tau8i + tau9r;

			x.re[index5] = tau8r + tau9i;
			x.im[index5] = tau8i - tau9r;


		}

	}

}

template <typename T>
void inline radix8_mixed_dif_inplace(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 8;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T wl5r = (T) 1.0;
	T wl5i = (T) 0.0;

	T wl6r = (T) 1.0;
	T wl6i = (T) 0.0;

	T wl7r = (T) 1.0;
	T wl7i = (T) 0.0;

	T c1 = 0.70710678118654752440084436210485;
	//T c2 = -0.70710678118654752440084436210485;
	T s1 = 0.70710678118654752440084436210485;
	//T s2 = 0.70710678118654752440084436210485;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T tau8r,tau8i,tau9r,tau9i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi,hr,hi;
	T temp1r,temp1i,temp2r,temp2i;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;
			int index7 = index6+Ls;

			tau0r = x.re[index] + x.re[index4];
			tau0i = x.im[index] + x.im[index4];

			tau1r = x.re[index1] + x.re[index7];
			tau1i = x.im[index1] + x.im[index7];

			tau2r = x.re[index3] + x.re[index5];
			tau2i = x.im[index3] + x.im[index5];

			tau3r = x.re[index2] + x.re[index6];
			tau3i = x.im[index2] + x.im[index6];

			tau4r = x.re[index] - x.re[index4];
			tau4i = x.im[index] - x.im[index4];

			tau5r = x.re[index1] - x.re[index7];
			tau5i = x.im[index1] - x.im[index7];

			tau6r = x.re[index3] - x.re[index5];
			tau6i = x.im[index3] - x.im[index5];

			tau7r = x.re[index2] - x.re[index6];
			tau7i = x.im[index2] - x.im[index6];

			ar = tau0r + tau1r + tau2r + tau3r;
			ai = tau0i + tau1i + tau2i + tau3i;

			er = tau0r - tau1r - tau2r + tau3r;
			ei = tau0i - tau1i - tau2i + tau3i;

			temp1r = tau1r - tau2r;
			temp1i = tau1i - tau2i;

			temp2r = tau5r + tau6r;
			temp2i = tau5i + tau6i;

			tau8r =  tau4r + c1 * temp1r;
			tau8i =  tau4i + c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r - tau7r);
			tau9i = sgn * ( -s1 * temp2i - tau7i);

			br = tau8r - tau9i;
			bi = tau8i + tau9r;

			hr = tau8r + tau9i;
			hi = tau8i - tau9r;

			tau8r = tau0r - tau3r;
			tau8i = tau0i - tau3i;

			tau9r = sgn * ( -tau5r + tau6r);
			tau9i = sgn * ( -tau5i + tau6i);

			cr = tau8r - tau9i;
			ci = tau8i + tau9r;

			gr = tau8r + tau9i;
			gi = tau8i - tau9r;

			tau8r = tau4r - c1 * temp1r;
			tau8i = tau4i - c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r + tau7r);
			tau9i = sgn * ( -s1 * temp2i + tau7i);

			dr = tau8r - tau9i;
			di = tau8i + tau9r;

			fr = tau8r + tau9i;
			fi = tau8i - tau9r;

			x.re[index] = ar;
			x.im[index] = ai;

			x.re[index1] = wlr*br - wli*bi;
			x.im[index1] = wlr*bi + wli*br;

			x.re[index2] = wl2r*cr - wl2i*ci;
			x.im[index2] = wl2r*ci + wl2i*cr;

			x.re[index3] = wl3r*dr - wl3i*di;
			x.im[index3] = wl3r*di + wl3i*dr;

			x.re[index4] = wl4r*er - wl4i*ei;
			x.im[index4] = wl4r*ei + wl4i*er;

			x.re[index5] = wl5r*fr - wl5i*fi;
			x.im[index5] = wl5r*fi + wl5i*fr;

			x.re[index6] = wl6r*gr - wl6i*gi;
			x.im[index6] = wl6r*gi + wl6i*gr;

			x.re[index7] = wl7r*hr - wl7i*hi;
			x.im[index7] = wl7r*hi + wl7i*hr;
		}
		T temp = wlr;
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

template <typename T>
void inline radix8_mixed_dif_inplace(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 8;
	int r = n / L;

	T c1 = 0.70710678118654752440084436210485;
	//T c2 = -0.70710678118654752440084436210485;
	T s1 = 0.70710678118654752440084436210485;
	//T s2 = 0.70710678118654752440084436210485;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T tau8r,tau8i,tau9r,tau9i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi,hr,hi;
	T temp1r,temp1i,temp2r,temp2i;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i,wl5r,wl5i,wl6r,wl6i,wl7r,wl7i;


	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

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

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;
			int index7 = index6+Ls;

			tau0r = x.re[index] + x.re[index4];
			tau0i = x.im[index] + x.im[index4];

			tau1r = x.re[index1] + x.re[index7];
			tau1i = x.im[index1] + x.im[index7];

			tau2r = x.re[index3] + x.re[index5];
			tau2i = x.im[index3] + x.im[index5];

			tau3r = x.re[index2] + x.re[index6];
			tau3i = x.im[index2] + x.im[index6];

			tau4r = x.re[index] - x.re[index4];
			tau4i = x.im[index] - x.im[index4];

			tau5r = x.re[index1] - x.re[index7];
			tau5i = x.im[index1] - x.im[index7];

			tau6r = x.re[index3] - x.re[index5];
			tau6i = x.im[index3] - x.im[index5];

			tau7r = x.re[index2] - x.re[index6];
			tau7i = x.im[index2] - x.im[index6];

			ar = tau0r + tau1r + tau2r + tau3r;
			ai = tau0i + tau1i + tau2i + tau3i;

			er = tau0r - tau1r - tau2r + tau3r;
			ei = tau0i - tau1i - tau2i + tau3i;

			temp1r = tau1r - tau2r;
			temp1i = tau1i - tau2i;

			temp2r = tau5r + tau6r;
			temp2i = tau5i + tau6i;

			tau8r =  tau4r + c1 * temp1r;
			tau8i =  tau4i + c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r - tau7r);
			tau9i = sgn * ( -s1 * temp2i - tau7i);

			br = tau8r - tau9i;
			bi = tau8i + tau9r;

			hr = tau8r + tau9i;
			hi = tau8i - tau9r;

			tau8r = tau0r - tau3r;
			tau8i = tau0i - tau3i;

			tau9r = sgn * ( -tau5r + tau6r);
			tau9i = sgn * ( -tau5i + tau6i);

			cr = tau8r - tau9i;
			ci = tau8i + tau9r;

			gr = tau8r + tau9i;
			gi = tau8i - tau9r;

			tau8r = tau4r - c1 * temp1r;
			tau8i = tau4i - c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r + tau7r);
			tau9i = sgn * ( -s1 * temp2i + tau7i);

			dr = tau8r - tau9i;
			di = tau8i + tau9r;

			fr = tau8r + tau9i;
			fi = tau8i - tau9r;

			x.re[index] = ar;
			x.im[index] = ai;

			x.re[index1] = wlr*br - wli*bi;
			x.im[index1] = wlr*bi + wli*br;

			x.re[index2] = wl2r*cr - wl2i*ci;
			x.im[index2] = wl2r*ci + wl2i*cr;

			x.re[index3] = wl3r*dr - wl3i*di;
			x.im[index3] = wl3r*di + wl3i*dr;

			x.re[index4] = wl4r*er - wl4i*ei;
			x.im[index4] = wl4r*ei + wl4i*er;

			x.re[index5] = wl5r*fr - wl5i*fi;
			x.im[index5] = wl5r*fi + wl5i*fr;

			x.re[index6] = wl6r*gr - wl6i*gi;
			x.im[index6] = wl6r*gi + wl6i*gr;

			x.re[index7] = wl7r*hr - wl7i*hi;
			x.im[index7] = wl7r*hi + wl7i*hr;
		}


	}

}

template <typename T>
void inline sh_mixed_radix8_dit(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 8;
	int r = n / L;
	int rs = n / Ls;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T wl5r = (T) 1.0;
	T wl5i = (T) 0.0;

	T wl6r = (T) 1.0;
	T wl6i = (T) 0.0;

	T wl7r = (T) 1.0;
	T wl7i = (T) 0.0;

	T c1 = 0.70710678118654752440084436210485;
	//T c2 = -0.70710678118654752440084436210485;
	T s1 = 0.70710678118654752440084436210485;
	//T s2 = 0.70710678118654752440084436210485;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T tau8r,tau8i,tau9r,tau9i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi,hr,hi;
	T temp1r,temp1i,temp2r,temp2i;

	fft_data<T> y = x;
	int lsr = Ls * r;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = j*rs + k;
			int index1 = index+r;
			int index2 = index1+r;
			int index3 = index2+r;
			int index4 = index3+r;
			int index5 = index4+r;
			int index6 = index5+r;
			int index7 = index6+r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];

			dr = wl3r*y.re[index3] - wl3i*y.im[index3];
			di = wl3r*y.im[index3] + wl3i*y.re[index3];

			er = wl4r*y.re[index4] - wl4i*y.im[index4];
			ei = wl4r*y.im[index4] + wl4i*y.re[index4];

			fr = wl5r*y.re[index5] - wl5i*y.im[index5];
			fi = wl5r*y.im[index5] + wl5i*y.re[index5];

			gr = wl6r*y.re[index6] - wl6i*y.im[index6];
			gi = wl6r*y.im[index6] + wl6i*y.re[index6];

			hr = wl7r*y.re[index7] - wl7i*y.im[index7];
			hi = wl7r*y.im[index7] + wl7i*y.re[index7];

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

			int indexo = j*r + k;
			int indexo1 = indexo + lsr;
			int indexo2 = indexo1 + lsr;
			int indexo3 = indexo2 + lsr;
			int indexo4 = indexo3 + lsr;
			int indexo5 = indexo4 + lsr;
			int indexo6 = indexo5 + lsr;
			int indexo7 = indexo6 + lsr;

			x.re[indexo] = tau0r + tau1r + tau2r + tau3r;
			x.im[indexo] = tau0i + tau1i + tau2i + tau3i;

			x.re[indexo4] = tau0r - tau1r - tau2r + tau3r;
			x.im[indexo4] = tau0i - tau1i - tau2i + tau3i;

			temp1r = tau1r - tau2r;
			temp1i = tau1i - tau2i;

			temp2r = tau5r + tau6r;
			temp2i = tau5i + tau6i;

			tau8r =  tau4r + c1 * temp1r;
			tau8i =  tau4i + c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r - tau7r);
			tau9i = sgn * ( -s1 * temp2i - tau7i);

			x.re[indexo1] = tau8r - tau9i;
			x.im[indexo1] = tau8i + tau9r;

			x.re[indexo7] = tau8r + tau9i;
			x.im[indexo7] = tau8i - tau9r;

			tau8r = tau0r - tau3r;
			tau8i = tau0i - tau3i;

			tau9r = sgn * ( -tau5r + tau6r);
			tau9i = sgn * ( -tau5i + tau6i);

			x.re[indexo2] = tau8r - tau9i;
			x.im[indexo2] = tau8i + tau9r;

			x.re[indexo6] = tau8r + tau9i;
			x.im[indexo6] = tau8i - tau9r;

			tau8r = tau4r - c1 * temp1r;
			tau8i = tau4i - c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r + tau7r);
			tau9i = sgn * ( -s1 * temp2i + tau7i);

			x.re[indexo3] = tau8r - tau9i;
			x.im[indexo3] = tau8i + tau9r;

			x.re[indexo5] = tau8r + tau9i;
			x.im[indexo5] = tau8i - tau9r;


		}
		T temp = wlr;
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

template <typename T>
void inline sh_mixed_radix8_dit(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 8;
	int r = n / L;
	int rs = n / Ls;

	T c1 = 0.70710678118654752440084436210485;
	//T c2 = -0.70710678118654752440084436210485;
	T s1 = 0.70710678118654752440084436210485;
	//T s2 = 0.70710678118654752440084436210485;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T tau8r,tau8i,tau9r,tau9i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi,hr,hi;
	T temp1r,temp1i,temp2r,temp2i;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i,wl5r,wl5i,wl6r,wl6i,wl7r,wl7i;

	fft_data<T> y = x;
	int lsr = Ls * r;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

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

		for (int k =0; k < r; k++) {
			int index = j*rs + k;
			int index1 = index+r;
			int index2 = index1+r;
			int index3 = index2+r;
			int index4 = index3+r;
			int index5 = index4+r;
			int index6 = index5+r;
			int index7 = index6+r;

			ar = y.re[index];
			ai = y.im[index];

			br = wlr*y.re[index1] - wli*y.im[index1];
			bi = wlr*y.im[index1] + wli*y.re[index1];

			cr = wl2r*y.re[index2] - wl2i*y.im[index2];
			ci = wl2r*y.im[index2] + wl2i*y.re[index2];

			dr = wl3r*y.re[index3] - wl3i*y.im[index3];
			di = wl3r*y.im[index3] + wl3i*y.re[index3];

			er = wl4r*y.re[index4] - wl4i*y.im[index4];
			ei = wl4r*y.im[index4] + wl4i*y.re[index4];

			fr = wl5r*y.re[index5] - wl5i*y.im[index5];
			fi = wl5r*y.im[index5] + wl5i*y.re[index5];

			gr = wl6r*y.re[index6] - wl6i*y.im[index6];
			gi = wl6r*y.im[index6] + wl6i*y.re[index6];

			hr = wl7r*y.re[index7] - wl7i*y.im[index7];
			hi = wl7r*y.im[index7] + wl7i*y.re[index7];

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

			int indexo = j*r + k;
			int indexo1 = indexo + lsr;
			int indexo2 = indexo1 + lsr;
			int indexo3 = indexo2 + lsr;
			int indexo4 = indexo3 + lsr;
			int indexo5 = indexo4 + lsr;
			int indexo6 = indexo5 + lsr;
			int indexo7 = indexo6 + lsr;

			x.re[indexo] = tau0r + tau1r + tau2r + tau3r;
			x.im[indexo] = tau0i + tau1i + tau2i + tau3i;

			x.re[indexo4] = tau0r - tau1r - tau2r + tau3r;
			x.im[indexo4] = tau0i - tau1i - tau2i + tau3i;

			temp1r = tau1r - tau2r;
			temp1i = tau1i - tau2i;

			temp2r = tau5r + tau6r;
			temp2i = tau5i + tau6i;

			tau8r =  tau4r + c1 * temp1r;
			tau8i =  tau4i + c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r - tau7r);
			tau9i = sgn * ( -s1 * temp2i - tau7i);

			x.re[indexo1] = tau8r - tau9i;
			x.im[indexo1] = tau8i + tau9r;

			x.re[indexo7] = tau8r + tau9i;
			x.im[indexo7] = tau8i - tau9r;

			tau8r = tau0r - tau3r;
			tau8i = tau0i - tau3i;

			tau9r = sgn * ( -tau5r + tau6r);
			tau9i = sgn * ( -tau5i + tau6i);

			x.re[indexo2] = tau8r - tau9i;
			x.im[indexo2] = tau8i + tau9r;

			x.re[indexo6] = tau8r + tau9i;
			x.im[indexo6] = tau8i - tau9r;

			tau8r = tau4r - c1 * temp1r;
			tau8i = tau4i - c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r + tau7r);
			tau9i = sgn * ( -s1 * temp2i + tau7i);

			x.re[indexo3] = tau8r - tau9i;
			x.im[indexo3] = tau8i + tau9r;

			x.re[indexo5] = tau8r + tau9i;
			x.im[indexo5] = tau8i - tau9r;


		}

	}

}

template <typename T>
void inline sh_mixed_radix8_dif(fft_data<T> &x, int q, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 8;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	T wlr = (T) 1.0;
	T wli = (T) 0.0;

	T wl2r = (T) 1.0;
	T wl2i = (T) 0.0;

	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;

	T wl4r = (T) 1.0;
	T wl4i = (T) 0.0;

	T wl5r = (T) 1.0;
	T wl5i = (T) 0.0;

	T wl6r = (T) 1.0;
	T wl6i = (T) 0.0;

	T wl7r = (T) 1.0;
	T wl7i = (T) 0.0;

	T c1 = 0.70710678118654752440084436210485;
	//T c2 = -0.70710678118654752440084436210485;
	T s1 = 0.70710678118654752440084436210485;
	//T s2 = 0.70710678118654752440084436210485;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T tau8r,tau8i,tau9r,tau9i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi,hr,hi;
	T temp1r,temp1i,temp2r,temp2i;

	fft_data<T> y = x;
	int lsr = Ls * r;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;
			int index7 = index6+Ls;

			tau0r = y.re[index] + y.re[index4];
			tau0i = y.im[index] + y.im[index4];

			tau1r = y.re[index1] + y.re[index7];
			tau1i = y.im[index1] + y.im[index7];

			tau2r = y.re[index3] + y.re[index5];
			tau2i = y.im[index3] + y.im[index5];

			tau3r = y.re[index2] + y.re[index6];
			tau3i = y.im[index2] + y.im[index6];

			tau4r = y.re[index] - y.re[index4];
			tau4i = y.im[index] - y.im[index4];

			tau5r = y.re[index1] - y.re[index7];
			tau5i = y.im[index1] - y.im[index7];

			tau6r = y.re[index3] - y.re[index5];
			tau6i = y.im[index3] - y.im[index5];

			tau7r = y.re[index2] - y.re[index6];
			tau7i = y.im[index2] - y.im[index6];

			ar = tau0r + tau1r + tau2r + tau3r;
			ai = tau0i + tau1i + tau2i + tau3i;

			er = tau0r - tau1r - tau2r + tau3r;
			ei = tau0i - tau1i - tau2i + tau3i;

			temp1r = tau1r - tau2r;
			temp1i = tau1i - tau2i;

			temp2r = tau5r + tau6r;
			temp2i = tau5i + tau6i;

			tau8r =  tau4r + c1 * temp1r;
			tau8i =  tau4i + c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r - tau7r);
			tau9i = sgn * ( -s1 * temp2i - tau7i);

			br = tau8r - tau9i;
			bi = tau8i + tau9r;

			hr = tau8r + tau9i;
			hi = tau8i - tau9r;

			tau8r = tau0r - tau3r;
			tau8i = tau0i - tau3i;

			tau9r = sgn * ( -tau5r + tau6r);
			tau9i = sgn * ( -tau5i + tau6i);

			cr = tau8r - tau9i;
			ci = tau8i + tau9r;

			gr = tau8r + tau9i;
			gi = tau8i - tau9r;

			tau8r = tau4r - c1 * temp1r;
			tau8i = tau4i - c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r + tau7r);
			tau9i = sgn * ( -s1 * temp2i + tau7i);

			dr = tau8r - tau9i;
			di = tau8i + tau9r;

			fr = tau8r + tau9i;
			fi = tau8i - tau9r;

			int indexo = k*Ls+j;
			int indexo1 = indexo+lsr;
			int indexo2 = indexo1+lsr;
			int indexo3 = indexo2+lsr;
			int indexo4 = indexo3+lsr;
			int indexo5 = indexo4+lsr;
			int indexo6 = indexo5+lsr;
			int indexo7 = indexo6+lsr;

			x.re[indexo] = ar;
			x.im[indexo] = ai;

			x.re[indexo1] = wlr*br - wli*bi;
			x.im[indexo1] = wlr*bi + wli*br;

			x.re[indexo2] = wl2r*cr - wl2i*ci;
			x.im[indexo2] = wl2r*ci + wl2i*cr;

			x.re[indexo3] = wl3r*dr - wl3i*di;
			x.im[indexo3] = wl3r*di + wl3i*dr;

			x.re[indexo4] = wl4r*er - wl4i*ei;
			x.im[indexo4] = wl4r*ei + wl4i*er;

			x.re[indexo5] = wl5r*fr - wl5i*fi;
			x.im[indexo5] = wl5r*fi + wl5i*fr;

			x.re[indexo6] = wl6r*gr - wl6i*gi;
			x.im[indexo6] = wl6r*gi + wl6i*gr;

			x.re[indexo7] = wl7r*hr - wl7i*hi;
			x.im[indexo7] = wl7r*hi + wl7i*hr;
		}
		T temp = wlr;
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

template <typename T>
void inline sh_mixed_radix8_dif(fft_data<T> &x, fft_data<T> &wl, int L, int sgn) {
	int n = x.re.size();
	int Ls = L / 8;
	int r = n / L;


	T c1 = 0.70710678118654752440084436210485;
	//T c2 = -0.70710678118654752440084436210485;
	T s1 = 0.70710678118654752440084436210485;
	//T s2 = 0.70710678118654752440084436210485;

	T tau0r,tau0i,tau1r,tau1i,tau2r,tau2i,tau3r,tau3i;
	T tau4r,tau4i,tau5r,tau5i;
	T tau6r,tau6i,tau7r,tau7i;
	T tau8r,tau8i,tau9r,tau9i;
	T ar,ai,br,bi,cr,ci,dr,di,er,ei;
	T fr,fi,gr,gi,hr,hi;
	T temp1r,temp1i,temp2r,temp2i;
	T wlr,wli,wl2r,wl2i,wl3r,wl3i,wl4r,wl4i,wl5r,wl5i,wl6r,wl6i,wl7r,wl7i;


	fft_data<T> y = x;
	int lsr = Ls * r;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr = wl.re[ind];
		wli = wl.im[ind];

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

		for (int k =0; k < r; k++) {
			int index = k*L+j;
			int index1 = index+Ls;
			int index2 = index1+Ls;
			int index3 = index2+Ls;
			int index4 = index3+Ls;
			int index5 = index4+Ls;
			int index6 = index5+Ls;
			int index7 = index6+Ls;

			tau0r = y.re[index] + y.re[index4];
			tau0i = y.im[index] + y.im[index4];

			tau1r = y.re[index1] + y.re[index7];
			tau1i = y.im[index1] + y.im[index7];

			tau2r = y.re[index3] + y.re[index5];
			tau2i = y.im[index3] + y.im[index5];

			tau3r = y.re[index2] + y.re[index6];
			tau3i = y.im[index2] + y.im[index6];

			tau4r = y.re[index] - y.re[index4];
			tau4i = y.im[index] - y.im[index4];

			tau5r = y.re[index1] - y.re[index7];
			tau5i = y.im[index1] - y.im[index7];

			tau6r = y.re[index3] - y.re[index5];
			tau6i = y.im[index3] - y.im[index5];

			tau7r = y.re[index2] - y.re[index6];
			tau7i = y.im[index2] - y.im[index6];

			ar = tau0r + tau1r + tau2r + tau3r;
			ai = tau0i + tau1i + tau2i + tau3i;

			er = tau0r - tau1r - tau2r + tau3r;
			ei = tau0i - tau1i - tau2i + tau3i;

			temp1r = tau1r - tau2r;
			temp1i = tau1i - tau2i;

			temp2r = tau5r + tau6r;
			temp2i = tau5i + tau6i;

			tau8r =  tau4r + c1 * temp1r;
			tau8i =  tau4i + c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r - tau7r);
			tau9i = sgn * ( -s1 * temp2i - tau7i);

			br = tau8r - tau9i;
			bi = tau8i + tau9r;

			hr = tau8r + tau9i;
			hi = tau8i - tau9r;

			tau8r = tau0r - tau3r;
			tau8i = tau0i - tau3i;

			tau9r = sgn * ( -tau5r + tau6r);
			tau9i = sgn * ( -tau5i + tau6i);

			cr = tau8r - tau9i;
			ci = tau8i + tau9r;

			gr = tau8r + tau9i;
			gi = tau8i - tau9r;

			tau8r = tau4r - c1 * temp1r;
			tau8i = tau4i - c1 * temp1i;

			tau9r = sgn * ( -s1 * temp2r + tau7r);
			tau9i = sgn * ( -s1 * temp2i + tau7i);

			dr = tau8r - tau9i;
			di = tau8i + tau9r;

			fr = tau8r + tau9i;
			fi = tau8i - tau9r;

			int indexo = k*Ls+j;
			int indexo1 = indexo+lsr;
			int indexo2 = indexo1+lsr;
			int indexo3 = indexo2+lsr;
			int indexo4 = indexo3+lsr;
			int indexo5 = indexo4+lsr;
			int indexo6 = indexo5+lsr;
			int indexo7 = indexo6+lsr;

			x.re[indexo] = ar;
			x.im[indexo] = ai;

			x.re[indexo1] = wlr*br - wli*bi;
			x.im[indexo1] = wlr*bi + wli*br;

			x.re[indexo2] = wl2r*cr - wl2i*ci;
			x.im[indexo2] = wl2r*ci + wl2i*cr;

			x.re[indexo3] = wl3r*dr - wl3i*di;
			x.im[indexo3] = wl3r*di + wl3i*dr;

			x.re[indexo4] = wl4r*er - wl4i*ei;
			x.im[indexo4] = wl4r*ei + wl4i*er;

			x.re[indexo5] = wl5r*fr - wl5i*fi;
			x.im[indexo5] = wl5r*fi + wl5i*fr;

			x.re[indexo6] = wl6r*gr - wl6i*gi;
			x.im[indexo6] = wl6r*gi + wl6i*gr;

			x.re[indexo7] = wl7r*hr - wl7i*hi;
			x.im[indexo7] = wl7r*hi + wl7i*hr;
		}

	}

}

#endif /* CODELETS2_H_ */
