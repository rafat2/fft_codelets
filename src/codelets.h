//============================================================================
// Name : FFT
// Author : Rafat Hussain
// Version :
// Copyright : GNU GPL License
// Description : FFT Library Component
//============================================================================
/*
* Copyright (c) 2012 Rafat Hussain
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*/

#ifndef CODELETS_H
#define CODELETS_H

#include "hsfft_alg.h"

typedef vector<double>::iterator dblitr;

typedef vector<float>::iterator fltitr;

template <typename T>
void inline radix2_dit(fft_data<T> &x,fft_data<T> &wl,int q) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
	int Ls = L / 2;
	int r = n / L;
	
	T taur,taui;
	
	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j+Ls;
			taur = wl.re[j]*x.re[index] - wl.im[j]*x.im[index];
			taui = wl.re[j]*x.im[index] + wl.im[j]*x.re[index];
			
			x.re[k*L+j+Ls] = x.re[k*L+j] - taur;
			x.im[k*L+j+Ls] = x.im[k*L+j] - taui;
			
			x.re[k*L+j] = x.re[k*L+j] + taur;
			x.im[k*L+j] = x.im[k*L+j] + taui;
			
			
		}
	}
}

template <typename T>
void inline radix2_dit_us(fft_data<T> &x,fft_data<T> &wl,int q) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
	int Ls = L / 2;
	int r = n / L;
	T taur,taui;
	
	for (int k = 0; k < r; k++) {
		for (int j = 0; j < Ls; j++) {
			taur = wl.re[Ls - 1 + j] *x.re[k * L + j + Ls] - wl.im[Ls - 1 + j] *x.im[k * L + j + Ls];
			taui = wl.re[Ls - 1 + j] *x.im[k * L + j + Ls] + wl.im[Ls - 1 + j] *x.re[k * L + j + Ls];
			
			x.re[k*L+j+Ls] = x.re[k*L+j] - taur;
			x.im[k*L+j+Ls] = x.im[k*L+j] - taui;
			
			x.re[k*L+j] = x.re[k*L+j] + taur;
			x.im[k*L+j] = x.im[k*L+j] + taui;
			
		}
	}
	
}

template <typename T>
void inline radix2_dit_inplace(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
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
void inline sh_radix2_dit_us(fft_data<T> &x,fft_data<T> &wl,int q) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
	int Ls = L / 2;
	int r = n / L;
	T taur,taui;
	
	fft_data<T> y = x;
	
	for (int k = 0; k < r; k++) {
		for (int j = 0; j < Ls; j++) {
			int index = (k+r) * Ls + j;
			taur = wl.re[Ls - 1 + j] *y.re[index] - wl.im[Ls - 1 + j] *y.im[index];
			taui = wl.re[Ls - 1 + j] *y.im[index] + wl.im[Ls - 1 + j] *y.re[index];
			
			x.re[k*L+j+Ls] = y.re[k*Ls+j] - taur;
			x.im[k*L+j+Ls] = y.im[k*Ls+j] - taui;
			
			x.re[k*L+j] = y.re[k*Ls+j] + taur;
			x.im[k*L+j] = y.im[k*Ls+j] + taui;
			
		}
	}
	
}

template <typename T>
void inline sh_radix2_dit(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
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
void inline radix2_dif_inplace(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
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
void inline sh_radix2_dif(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
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
void inline radix2_dif_us(fft_data<T> &x,fft_data<T> &wl,int q) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
	int Ls = L / 2;
	int r = n / L;
	T taur,taui;
	T xlsr,xlsi;
	
	for (int k = 0; k < r; k++) {
		for (int j = 0; j < Ls; j++) {
			int index1 = k*L+Ls+j;
			int index2 = k*L+j;
			taur = x.re[index1];
			taui = x.im[index1];
			
			xlsr = x.re[index2] - taur;
			xlsi = x.im[index2] - taui;
			
			x.re[index1] = wl.re[Ls - 1 + j] *xlsr - wl.im[Ls - 1 + j] *xlsi;
			x.im[index1] = wl.re[Ls - 1 + j] *xlsi + wl.im[Ls - 1 + j] *xlsr;
			
			
			x.re[index2] = x.re[index2] + taur;
			x.im[index2] = x.im[index2] + taui;
			
		}
	}
	
}


template <typename T>
void inline radix4_dit_inplace(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(4.0, (double)q);
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
void inline radix4_dif_inplace(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(4.0, (double)q);
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
void inline sh_radix4_dit(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(4.0, (double)q);
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
void inline sh_radix4_dif(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(4.0, (double)q);
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
void inline split_radix(fft_data<T> &x,vector<int> &beta, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
	int m = L / 2;
	int r = n / L;
	int p = L/4;
	
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;
	
	T S = (T) sin(theta);
	T C = (T) cos(theta);
	
	T PI6 = (T) 3.0*6.28318530717958647692528676655900577;
	T theta3 = (T) -1.0 * sgn * PI6/L;
	
	T S3 = (T) sin(theta3);
	T C3 = (T) cos(theta3);
	
	T wlr = (T) 1.0;
	T wli = (T) 0.0;
	
	//T wl2r = (T) 1.0;
	//T wl2i = (T) 0.0;
	
	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;
	
	T tau1r,tau1i,tau2r,tau2i;
	T ur,ui,vr,vi;
	
	if (sgn == 1) {
		for (int j = 0; j < p; j++) {
			for (int k =0; k < r; k++) {
				if (beta[k] == 1) {
					if (j == 0) {
						int index = k*L;
						int index1 = index+m;
						int index2 = index1+p;
						int index3 = index+p;
						
						tau1r = x.re[index1];
						tau1i = x.im[index1];
						
						tau2r = x.re[index2];
						tau2i = x.im[index2];
						
						ur = tau1r + tau2r;
						ui = tau1i + tau2i;
						
						vr = tau1r - tau2r;
						vi = tau1i - tau2i;
						
						x.re[index2] = x.re[index3] - vi;
						x.im[index2] = x.im[index3] + vr;
						
						x.re[index1] = x.re[index] - ur;
						x.im[index1] = x.im[index] - ui;
						
						x.re[index3] = x.re[index3] + vi;
						x.im[index3] = x.im[index3] - vr;
						
						x.re[index] = x.re[index] + ur;
						x.im[index] = x.im[index] + ui;
					
					} else {
						int index = k*L+j;
						int index1 = index+m;
						int index2 = index1+p;
						int index3 = index+p;
						tau1r = wlr * x.re[index1] - wli * x.im[index1];
						tau1i = wlr * x.im[index1] + wli * x.re[index1];
						
						tau2r = wl3r * x.re[index2] - wl3i * x.im[index2];
						tau2i = wl3r * x.im[index2] + wl3i * x.re[index2];
						
						ur = tau1r + tau2r;
						ui = tau1i + tau2i;
						
						vr = tau1r - tau2r;
						vi = tau1i - tau2i;
						
						x.re[index2] = x.re[index3] - vi;
						x.im[index2] = x.im[index3] + vr;
						
						x.re[index1] = x.re[index] - ur;
						x.im[index1] = x.im[index] - ui;
						
						x.re[index3] = x.re[index3] + vi;
						x.im[index3] = x.im[index3] - vr;
						
						x.re[index] = x.re[index] + ur;
						x.im[index] = x.im[index] + ui;
						
					}
					
				}
			}
			T temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			temp = wl3r;
			wl3r = C3 * wl3r - S3 * wl3i;
			wl3i = S3 * temp + C3 * wl3i;
			
			//wl2r = wlr*wlr - wli*wli;
			//wl2i = 2.0*wlr*wli;
			
			//wl3r = wl2r*wlr - wli*wl2i;
			//wl3i= wl2r*wli + wl2i*wlr;
		}
	} else if (sgn == -1) {
		for (int j = 0; j < p; j++) {
			for (int k =0; k < r; k++) {
				if (beta[k] == 1) {
					if (j == 0) {
						int index = k*L;
						int index1 = index+m;
						int index2 = index1+p;
						int index3 = index+p;
						tau1r = x.re[index1];
						tau1i = x.im[index1];
						
						tau2r = x.re[index2];
						tau2i = x.im[index2];
						
						ur = tau1r + tau2r;
						ui = tau1i + tau2i;
						
						vr = tau2r - tau1r;
						vi = tau2i - tau1i;
						
						x.re[index2] = x.re[index3] - vi;
						x.im[index2] = x.im[index3] + vr;
						
						x.re[index1] = x.re[index] - ur;
						x.im[index1] = x.im[index] - ui;
						
						x.re[index3] = x.re[index3] + vi;
						x.im[index3] = x.im[index3] - vr;
						
						x.re[index] = x.re[index] + ur;
						x.im[index] = x.im[index] + ui;
					
					} else {
						int index = k*L+j;
						int index1 = index+m;
						int index2 = index1+p;
						int index3 = index+p;
						tau1r = wlr * x.re[index1] - wli * x.im[index1];
						tau1i = wlr * x.im[index1] + wli * x.re[index1];
						
						tau2r = wl3r * x.re[index2] - wl3i * x.im[index2];
						tau2i = wl3r * x.im[index2] + wl3i * x.re[index2];
						
						ur = tau1r + tau2r;
						ui = tau1i + tau2i;
						
						vr = tau2r - tau1r;
						vi = tau2i - tau1i;
						
						x.re[index2] = x.re[index3] - vi;
						x.im[index2] = x.im[index3] + vr;
						
						x.re[index1] = x.re[index] - ur;
						x.im[index1] = x.im[index] - ui;
						
						x.re[index3] = x.re[index3] + vi;
						x.im[index3] = x.im[index3] - vr;
						
						x.re[index] = x.re[index] + ur;
						x.im[index] = x.im[index] + ui;
						
					}
					
				}
			}
			T temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			temp = wl3r;
			wl3r = C3 * wl3r - S3 * wl3i;
			wl3i = S3 * temp + C3 * wl3i;
			
			//wl2r = wlr*wlr - wli*wli;
			//wl2i = 2.0*wlr*wli;
			
			//wl3r = wl2r*wlr - wli*wl2i;
			//wl3i= wl2r*wli + wl2i*wlr;
		}
	}
	
}


void inline split_radix_rec(vector<double>::iterator r, vector<double>::iterator i, int sgn,int N) {
	if (N == 1) {
		return;
	} else if (N == 2) {

		double taur =  *(r+1);
		double taui =  *(i+1);
			
		*(r+1) = *r - taur; 
		*(i+1) = *i - taui; 
			
		*r = *r + taur; 
		*i = *i + taui; 

	} else if (N == 4) {
	
		double taur =  *(r+1);
		double taui =  *(i+1);
			
		*(r+1) = *r - taur; 
		*(i+1) = *i - taui; 
			
		*r = *r + taur; 
		*i = *i + taui; 
		
		double ur = *(r+2) + *(r+3);
		double ui = *(i+2) + *(i+3);
						
		double vr = sgn* (*(r+2) - *(r+3));
		double vi = sgn* (*(i+2) - *(i+3));
						
		*(r+3) = *(r+1) - vi;
		*(i+3) = *(i+1) + vr;
						
		*(r+2) = *r - ur;
		*(i+2) = *i - ui;
						
		*(r+1) = *(r+1) + vi;
		*(i+1) = *(i+1) - vr;
						
		*r = *r + ur;
		*i = *i + ui;
		
		//N = 2;
	
	} else if (N == 8) {
		//r 4
		double taur =  *(r+1);
		double taui =  *(i+1);
			
		*(r+1) = *r - taur; 
		*(i+1) = *i - taui; 
			
		*r = *r + taur; 
		*i = *i + taui; 
		
		double ur = *(r+2) + *(r+3);
		double ui = *(i+2) + *(i+3);
						
		double vr = sgn* (*(r+2) - *(r+3));
		double vi = sgn* (*(i+2) - *(i+3));
						
		*(r+3) = *(r+1) - vi;
		*(i+3) = *(i+1) + vr;
						
		*(r+2) = *r - ur;
		*(i+2) = *i - ui;
						
		*(r+1) = *(r+1) + vi;
		*(i+1) = *(i+1) - vr;
						
		*r = *r + ur;
		*i = *i + ui;
		//r+4 2
		taur =  *(r+5);
		taui =  *(i+5);
			
		*(r+5) = *(r+4) - taur; 
		*(i+5) = *(i+4) - taui; 
			
		*(r+4) = *(r+4) + taur; 
		*(i+4) = *(i+4) + taui;
//r+6 2
		taur =  *(r+7);
		taui =  *(i+7);
			
		*(r+7) = *(r+6) - taur; 
		*(i+7) = *(i+6) - taui; 
			
		*(r+6) = *(r+6) + taur; 
		*(i+6) = *(i+6) + taui;

		//index1 = 4, index2 = 6, index3 = 2

		double tau1r = *(r+4);
		double tau1i = *(i+4);

		double tau2r = *(r+6);
		double tau2i = *(i+6);

		ur = *(r+4) + *(r+6);
		ui = *(i+4) + *(i+6);

		vr = sgn * (*(r+4) - *(r+6));
		vi = sgn * (*(i+4) - *(i+6));

		*(r+6) = *(r+2) - vi;
		*(i+6) = *(i+2) + vr;
						
		*(r+4) = *r - ur;
		*(i+4) = *i - ui;
						
		*(r+2) = *(r+2) + vi;
		*(i+2) = *(i+2) - vr;
						
		*r = *r + ur;
		*i = *i + ui;


		double wlr = 0.70710678118654752440084436210485;
		double wli = sgn * -0.70710678118654752440084436210485;

		double wl3r = -0.70710678118654752440084436210485;
		double wl3i = sgn * -0.70710678118654752440084436210485;

		//index1 = 5, index2 = 7, index3 = 3;

		tau1r = *(r+5)*wlr - *(i+5)*wli;
		tau1i = *(i+5)*wlr + *(r+5)*wli;

		tau2r = *(r+7)*wl3r - *(i+7)*wl3i;
		tau2i = *(i+7)*wl3r + *(r+7)*wl3i;

		ur = tau1r + tau2r;
		ui = tau1i + tau2i;
						
		vr = sgn* (tau1r - tau2r);
		vi = sgn* (tau1i - tau2i);

		*(r+7) = *(r+3) - vi;
		*(i+7) = *(i+3) + vr;
						
		*(r+5) = *(r+1) - ur;
		*(i+5) = *(i+1) - ui;
						
		*(r+3) = *(r+3) + vi;
		*(i+3) = *(i+3) - vr;
						
		*(r+1) = *(r+1) + ur;
		*(i+1) = *(i+1) + ui;
	
	} else {
		int m = N/2;
		int p = N/4;
		
		split_radix_rec(r,i,sgn,m);
		split_radix_rec(r+m,i+m,sgn,p);
		split_radix_rec(r+m+p,i+m+p,sgn,p);
		
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
			//wlr =  cos(theta*j);
		    //wli =  sin(theta*j);
			
			//wl3r =  cos(theta3*j);
		    //wl3i =  sin(theta3*j);
			
			int index1 = j+m;
			int index2 = index1+p;
			int index3 = j+p;
			
			cout << m << " " << *(r+j) << endl;

			tau1r = *(r+index1)*wlr - *(i+index1)*wli;
			tau1i = *(i+index1)*wlr + *(r+index1)*wli;
						
			tau2r = *(r+index2)*wl3r - *(i+index2)*wl3i;;
			tau2i = *(i+index2)*wl3r + *(r+index2)*wl3i;
						
			ur = tau1r + tau2r;
			ui = tau1i + tau2i;
						
			vr = sgn* (tau1r - tau2r);
			vi = sgn* (tau1i - tau2i);
						
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
		
	}
	
	
}


template <typename T>
void inline split_radix_dit(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
	int m = L / 2;
	int r = n / L;
	int p = L/4;
	
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;
	
	T S = (T) sin(theta);
	T C = (T) cos(theta);
	
	T PI6 = (T) 3.0*6.28318530717958647692528676655900577;
	T theta3 = (T) -1.0 * sgn * PI6/L;
	
	T S3 = (T) sin(theta3);
	T C3 = (T) cos(theta3);
	
	T wlr = (T) 1.0;
	T wli = (T) 0.0;
	
	//T wl2r = (T) 1.0;
	//T wl2i = (T) 0.0;
	
	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;
	
	T tau1r,tau1i,tau2r,tau2i;
	T ur,ui,vr,vi;
	
	
	if (sgn == 1) {
		for (int j = 0; j < p; j++) {
			int c1 = j;
			int step = 2*L;
			while (c1 < n-1) {
				for (int k =c1; k < n-1; k+=step) {
						int index = k;
						int index1 = index+m;
						int index2 = index1+p;
						int index3 = index+p;
						tau1r = wlr * x.re[index1] - wli * x.im[index1];
						tau1i = wlr * x.im[index1] + wli * x.re[index1];
						
						tau2r = wl3r * x.re[index2] - wl3i * x.im[index2];
						tau2i = wl3r * x.im[index2] + wl3i * x.re[index2];
						
						ur = tau1r + tau2r;
						ui = tau1i + tau2i;
						
						vr = tau1r - tau2r;
						vi = tau1i - tau2i;
						
						x.re[index2] = x.re[index3] - vi;
						x.im[index2] = x.im[index3] + vr;
						
						x.re[index1] = x.re[index] - ur;
						x.im[index1] = x.im[index] - ui;
						
						x.re[index3] = x.re[index3] + vi;
						x.im[index3] = x.im[index3] - vr;
						
						x.re[index] = x.re[index] + ur;
						x.im[index] = x.im[index] + ui;
						
				
					
				}
				c1 = 2*step - L + j;
				step = 4 * step;
			}
			T temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			temp = wl3r;
			wl3r = C3 * wl3r - S3 * wl3i;
			wl3i = S3 * temp + C3 * wl3i;
			
			//wl2r = wlr*wlr - wli*wli;
			//wl2i = 2.0*wlr*wli;
			
			//wl3r = wl2r*wlr - wli*wl2i;
			//wl3i= wl2r*wli + wl2i*wlr;
		}
	} else if (sgn == -1) {
		for (int j = 0; j < p; j++) {
			int c1 = j;
			int step = 2*L;
			while (c1 < n-1) {
				for (int k =c1; k < n-1; k+=step) {
						int index = k;
						int index1 = index+m;
						int index2 = index1+p;
						int index3 = index+p;
						tau1r = wlr * x.re[index1] - wli * x.im[index1];
						tau1i = wlr * x.im[index1] + wli * x.re[index1];
						
						tau2r = wl3r * x.re[index2] - wl3i * x.im[index2];
						tau2i = wl3r * x.im[index2] + wl3i * x.re[index2];
						
						ur = tau1r + tau2r;
						ui = tau1i + tau2i;
						
						vr = tau2r - tau1r;
						vi = tau2i - tau1i;
						
						x.re[index2] = x.re[index3] - vi;
						x.im[index2] = x.im[index3] + vr;
						
						x.re[index1] = x.re[index] - ur;
						x.im[index1] = x.im[index] - ui;
						
						x.re[index3] = x.re[index3] + vi;
						x.im[index3] = x.im[index3] - vr;
						
						x.re[index] = x.re[index] + ur;
						x.im[index] = x.im[index] + ui;
						
					
				}
				c1 = 2*step - L + j;
				step = 4 * step;
			}
			T temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			temp = wl3r;
			wl3r = C3 * wl3r - S3 * wl3i;
			wl3i = S3 * temp + C3 * wl3i;
			
		}
	}
	
}

template <typename T>
void inline split_radix_dif(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
	int m = L / 2;
	int p = L/4;
	
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;
	
	T S = (T) sin(theta);
	T C = (T) cos(theta);
	
	T PI6 = (T) 3.0*6.28318530717958647692528676655900577;
	T theta3 = (T) -1.0 * sgn * PI6/L;
	
	T S3 = (T) sin(theta3);
	T C3 = (T) cos(theta3);
	
	T wlr = (T) 1.0;
	T wli = (T) 0.0;
	
	
	T wl3r = (T) 1.0;
	T wl3i = (T) 0.0;
	
	T tau1r,tau1i,tau2r,tau2i;
	T ur,ui,vr,vi;
	
	
	if (sgn == 1) {
		for (int j = 0; j < p; j++) {
			int c1 = j;
			int step = 2*L;
			while (c1 < n-1) {
				for (int k =c1; k < n; k+=step) {
						int index = k;
						int index1 = index+m;
						int index2 = index1+p;
						int index3 = index+p;
						
						ur = x.re[index] - x.re[index1];
						ui = x.im[index] - x.im[index1];
						
						vr = x.re[index3] - x.re[index2];
						vi = x.im[index3] - x.im[index2];
						
						tau1r = ur + vi;
						tau1i = ui - vr;
						
						tau2r = ur - vi;
						tau2i = ui + vr;
						
						
						x.re[index3] = x.re[index3] + x.re[index2];
						x.im[index3] = x.im[index3] + x.im[index2];
						
						x.re[index] = x.re[index] + x.re[index1];
						x.im[index] = x.im[index] + x.im[index1];
						
						x.re[index1] = wlr * tau1r - wli * tau1i;
						x.im[index1] = wlr * tau1i + wli * tau1r;
						
						x.re[index2] = wl3r * tau2r - wl3i * tau2i;
						x.im[index2] = wl3r * tau2i + wl3i * tau2r;
						
					
				}
				c1 = 2*step - L + j;
				step = 4 * step;
			}
			T temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			temp = wl3r;
			wl3r = C3 * wl3r - S3 * wl3i;
			wl3i = S3 * temp + C3 * wl3i;
			
		}
	} else if (sgn == -1) {
		for (int j = 0; j < p; j++) {
			int c1 = j;
			int step = 2*L;
			while (c1 < n-1) {
				for (int k =c1; k < n; k+=step) {
						int index = k;
						int index1 = index+m;
						int index2 = index1+p;
						int index3 = index+p;
						ur = x.re[index] - x.re[index1];
						ui = x.im[index] - x.im[index1];
						
						vr = x.re[index3] - x.re[index2];
						vi = x.im[index3] - x.im[index2];
						
						tau1r = ur - vi;
						tau1i = ui + vr;
						
						tau2r = ur + vi;
						tau2i = ui - vr;
						
						
						x.re[index3] = x.re[index3] + x.re[index2];
						x.im[index3] = x.im[index3] + x.im[index2];
						
						x.re[index] = x.re[index] + x.re[index1];
						x.im[index] = x.im[index] + x.im[index1];
						
						x.re[index1] = wlr * tau1r - wli * tau1i;
						x.im[index1] = wlr * tau1i + wli * tau1r;
						
						x.re[index2] = wl3r * tau2r - wl3i * tau2i;
						x.im[index2] = wl3r * tau2i + wl3i * tau2r;

				}
				c1 = 2*step - L + j;
				step = 4 * step;
			}
			T temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;
			
			temp = wl3r;
			wl3r = C3 * wl3r - S3 * wl3i;
			wl3i = S3 * temp + C3 * wl3i;
			
		}
	}
	
}

template <typename T>
void inline bluestein_hl(fft_data<T> &hl, fft_data<T> &hlt, int len, int M) {
	T PI = (T) 3.1415926535897932384626433832795;
	T theta = PI / len;
	int l2 = 0;
	int len2 = 2 * len;
	T angle;
	hl.re.resize(M,(T) 0.0);
	hl.im.resize(M,(T) 0.0);
	//fft_data<T> hlt;
	hlt.re.resize(len,(T) 0.0);
	hlt.im.resize(len,(T) 0.0);
	
	for (int i = 0 ; i < len; ++i) {
		angle = theta * l2;
		hl.re[i]=hlt.re[i] = cos(angle);
		hl.im[i]=hlt.im[i] = sin(angle);
		l2+=2*i+1;
		while (l2 > len2) {
			l2-=len2;
		}
		
	}
	
	for (int i = M - len + 1; i < M; i++) {
		hl.re[i] = hlt.re[M-i];
		hl.im[i] = hlt.im[M-i];
	}
	
}

template <typename T>
void inline radix2_wfta(fft_data<T> &x, int stride,int N, int sgn) {

	T R1;
	int scale = 2;
	
	if (stride == N) 
		scale = 1;
	
	for (int i = 0; i < N; ++i) {
		int v = i*scale;
		int w = v+stride;
		R1 = x.re[v];
		x.re[v] = R1 + x.re[w];
		x.re[w] = R1 - x.re[w];
		
		R1 = x.im[v];
		x.im[v] = R1 + x.im[w];
		x.im[w] = R1 - x.im[w];
	}
	
}

template <typename T>
void inline radix3_wfta(fft_data<T> &x, int stride,int N, int sgn) {
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/3;
	
	T cs = cos(theta) - 1;
	T ss = sin(theta);

	T R1,R2,S1,S2;
	
	int scale = 3;
	
	if (stride == N) 
		scale = 1;	
		
	for (int i = 0; i < N; ++i) {
		int u = i * scale;
		int v = u+stride;
		int w = v+stride;
		
		R2 = ss * (x.re[v] - x.re[w]);
		R1 = x.re[v] + x.re[w];
		
		x.re[u]+=R1;
		R1 = x.re[u] + R1*cs;
		
		S2 = ss * (x.im[v] - x.im[w]);
		S1 = x.im[v] + x.im[w];
		
		x.im[u]+=S1;
		S1 = x.im[u] + S1*cs;
		
		// Use sgn below
		
		x.re[v] = R1 - S2;
		x.re[w] = R1 + S2;
		
		x.im[v] = S1 + R2;
		x.im[w] = S1 - R2;
		
		
	}
	
}

template <typename T>
void inline radix4_wfta(fft_data<T> &x, int stride,int N, int sgn) {

	T R1,R2,S1,S2;
	
	int scale = 4;
	
	if (stride == N) 
		scale = 1;	
		
	for (int i = 0; i < N; ++i) {
		int u = i * scale;
		int v = u+stride;
		int w = v+stride;
		int y = w+stride;
		
		R1 = x.re[u] + x.re[w];
		S1 = x.re[u] - x.re[w];
		R2 = x.re[v] + x.re[y];
		
		x.re[u] = R1 + R2;
		x.re[w] = R1 - R2;
		
		R1 = x.im[u] + x.im[w];
		S2 = x.im[u] - x.im[w];
		R2 = x.im[v] + x.im[y];
		
		x.im[u] = R1 + R2;
		x.im[w] = R1 - R2;
		
		R1 = x.re[v] - x.re[y];
		R2 = x.im[v] - x.im[y];

		
		// Use sgn below
		
		x.re[v] = S1 + R2;
		x.re[y] = S1 - R2;
		
		x.im[v] = S2 - R1;
		x.im[y] = S2 + R1;
		
		
	}
	
}

template <typename T>
void inline radix5_wfta(fft_data<T> &x, int stride,int N, int sgn) {
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) 1.0 * sgn * PI2/5;
	
	T b1 = 0.5 * (cos(theta) + cos(2*theta)) - 1;
	T b2 = 0.5 * (cos(theta) - cos(2*theta));
	T b3 = sin(theta);
	T b4 = sin(theta) + sin(2*theta);
	T b5 = -sin(theta) + sin(2*theta);

	T T0,T1,T2,T3,T4,T1i,T2i,T3i,T4i,A1,A2,A3,A4,A5,A1i,A2i,A3i,A4i,A5i;
	
	int scale = 5;
	
	if (stride == N) 
		scale = 1;	
		
	for (int i = 0; i < N; ++i) {
		int i0 = i * scale;
		int i1 = i0+stride;
		int i2 = i1+stride;
		int i3 = i2+stride;
		int i4 = i3+stride;
		
		T0 = x.re[i1] + x.re[i4];
		T1 = x.re[i2] + x.re[i3];
		
		A4 = x.re[i3] - x.re[i2];
		A5 = x.re[i1] - x.re[i4];
		A1 = T0 + T1;
		A2 = T0 - T1;
		A3 = A4 + A5;
		x.re[i0] = x.re[i0] + A1;
		
		T0 = x.im[i1] + x.im[i4];
		T1 = x.im[i2] + x.im[i3];
		
		A4i = x.im[i3] - x.im[i2];
		A5i = x.im[i1] - x.im[i4];
		A1i = T0 + T1;
		A2i = T0 - T1;
		A3i = A4 + A5;
		x.im[i0] = x.im[i0] + A1i;
		
		A1 = A1 * b1;
		A2 = A2 * b2;
		A3 = A3 * b3;
		A4 = A4 * b4;
		A5 = A5 * b5;
		
		T0 = x.re[i0] + A1;
		T1 = A3 - A4;
		T2 = A3 + A5;
		T3 = T0 + A2;
		T4 = T0 - A2;
		
		A1i = A1i * b1;
		A2i = A2i * b2;
		A3i = A3i * b3;
		A4i = A4i * b4;
		A5i = A5i * b5;
		
		T0 = x.im[i0] + A1i;
		T1i = A3i - A4i;
		T2i = A3i + A5i;
		T3i = T0 + A2i;
		T4i = T0 - A2i;
		
		
		// Use sgn below
		
		x.re[i1] = T3 + T1i;
		x.re[i2] = T4 + T2i;
		x.re[i4] = T3 - T1i;
		x.re[i3] = T4 - T2i;
		
		x.im[i1] = T3i - T1;
		x.im[i2] = T4i - T2;
		x.im[i4] = T3i + T1;
		x.im[i3] = T4i + T2;
		
	}
	
}

template <typename T>
void inline radix2_dit_stride(vector<double>::iterator xr,vector<double>::iterator xi,fft_data<T> &wl,int q, int stride, int n) {
	//int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
	int Ls = L / 2;
	int r = n / L;
	T taur,taui;
	
	
	for (int k = 0; k < r; k++) {
		
		for (int j = 0; j < Ls; j++) {
			int index1 = (k * L + j + Ls) * stride;
			int index2 = (k * L + j) * stride;
			//cout << index2 << " : "<< index1 << endl;
			taur = *(xr+index1)*wl.re[Ls - 1 + j] - *(xi+index1)*wl.im[Ls - 1 + j] ;
			taui = *(xi+index1)*wl.re[Ls - 1 + j] + *(xr+index1)*wl.im[Ls - 1 + j] ;
			
			*(xr+index1)= *(xr+index2) - taur;
			*(xi+index1)= *(xi+index2) - taui;
			
			*(xr+index2) = *(xr+index2) + taur;
			*(xi+index2) = *(xi+index2)  + taui;
			
		}
	}
	
}

template <typename T>
void inline radix2_dit_stride(vector<float>::iterator xr,vector<float>::iterator xi,fft_data<T> &wl,int q, int stride, int n) {
	//int n = x.re.size();
	int L = (int) pow(2.0, (double)q);
	int Ls = L / 2;
	int r = n / L;
	T taur,taui;
	
	
	for (int k = 0; k < r; k++) {
		
		for (int j = 0; j < Ls; j++) {
			int index1 = (k * L + j + Ls) * stride;
			int index2 = (k * L + j) * stride;
			//cout << index2 << " : "<< index1 << endl;
			taur = *(xr+index1)*wl.re[Ls - 1 + j] - *(xi+index1)*wl.im[Ls - 1 + j] ;
			taui = *(xi+index1)*wl.re[Ls - 1 + j] + *(xr+index1)*wl.im[Ls - 1 + j] ;
			
			*(xr+index1)= *(xr+index2) - taur;
			*(xi+index1)= *(xi+index2) - taui;
			
			*(xr+index2) = *(xr+index2) + taur;
			*(xi+index2) = *(xi+index2)  + taui;
			
		}
	}
	
}


template <typename T>
void inline radix3_dit_inplace(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(3.0, (double)q);
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
void inline radix3_dit_inplace(fft_data<T> &x,fft_data<T> &wl, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(3.0, (double)q);
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
void inline radix3_dif_inplace(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(3.0, (double)q);
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
void inline radix3_dif_inplace(fft_data<T> &x,fft_data<T> &wl, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(3.0, (double)q);
	int Ls = L / 3;
	int r = n / L;
	
	
	T tau0r,tau0i,tau1r,tau1i;
	T ar,ai,br,bi,cr,ci,xlsr,xlsi;
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
void inline sh_radix3_dit(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(3.0, (double)q);
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
void inline sh_radix3_dit(fft_data<T> &x,fft_data<T> &wl, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(3.0, (double)q);
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
void inline sh_radix3_dif(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(3.0, (double)q);
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
	T ar,ai,br,bi,cr,ci,xlsr,xlsi;
	
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
void inline sh_radix3_dif(fft_data<T> &x,fft_data<T> &wl, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(3.0, (double)q);
	int Ls = L / 3;
	int r = n / L;
	
	T tau0r,tau0i,tau1r,tau1i;
	T ar,ai,br,bi,cr,ci,xlsr,xlsi;
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
void inline radix5_dit_inplace(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(5.0, (double)q);
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
void inline radix5_dit_inplace(fft_data<T> &x,fft_data<T> &wl, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(5.0, (double)q);
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
void inline radix5_dif_inplace(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(5.0, (double)q);
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
void inline radix5_dif_inplace(fft_data<T> &x, fft_data<T> &wl, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(5.0, (double)q);
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
void inline sh_radix5_dit(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(5.0, (double)q);
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
void inline sh_radix5_dit(fft_data<T> &x,fft_data<T> &wl, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(5.0, (double)q);
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
void inline sh_radix5_dif(fft_data<T> &x, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(5.0, (double)q);
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
void inline sh_radix5_dif(fft_data<T> &x,fft_data<T> &wl, int q, int sgn) {
	int n = x.re.size();
	int L = (int) pow(5.0, (double)q);
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
void inline radix2_dit_rec(fft_data<T> &y, fft_data<T> &data, int sgn,int N) {
	if (N == 1) {
		return;
	} else {
		int m = N/2;
		fft_data<T> a0,a1,y0,y1;
		evenodd(data,a0,a1,m);
		radix2_dit_rec(y0,a0,sgn,m);
		radix2_dit_rec(y1,a1,sgn,m);

		T PI2 = 6.28318530717958647692528676655900577;
		T theta =  -1.0 * sgn * PI2/N;

		T S =  sin(theta);
		T C =  cos(theta);

		T wlr =  1.0;
		T wli =  0.0;

		//T wl2r = (T) 1.0;
		//T wl2i = (T) 0.0;


		T tau1r,tau1i;

		for (int j = 0; j < m; j++) {

			//int index1 = j+m;

			tau1r = y1.re[j]*wlr - y1.im[j]*wli;
			tau1i = y1.im[j]*wlr + y1.re[j]*wli;

			y.re[j] = y0.re[j] + tau1r;
			y.im[j] = y0.im[j] + tau1i;

			tau1r = y1.re[j+m]*wlr - y1.im[j+m]*wli;
			tau1i = y1.im[j+m]*wlr + y1.re[j+m]*wli;

			y.re[j+m] = y0.re[j+m] - tau1r;
			y.im[j+m] = y0.im[j+m] - tau1i;

			T temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;

		}

	}


}

void inline radix2_dit_rec(dblitr opr, dblitr opi,fft_data<double> &data, int sgn,int N) {
	
	if (N == 1) {
		*opr = data.re[0];
		*opi = data.im[0];
		//return;
	} else {
		int m = N / 2;
		fft_data<double> a0,a1;
		evenodd(data,a0,a1,m);
		radix2_dit_rec(opr,opi,a0,sgn,m);
		radix2_dit_rec(opr+m,opi+m,a1,sgn,m);

		double PI2 = 6.28318530717958647692528676655900577;
		double theta =  -1.0 * sgn * PI2/N;

		double S =  sin(theta);
		double C =  cos(theta);

		double wlr =  1.0;
		double wli =  0.0;

		double tau1r,tau1i,tau2r,tau2i;

		for (int k = 0; k < m; ++k) {
			tau1r = *(opr+k);
			tau1i = *(opi+k);

			tau2r = *(opr+k+m)*wlr - *(opi+k+m)*wli;
			tau2i = *(opi+k+m)*wlr + *(opr+k+m)*wli;

			*(opr+k) = tau1r + tau2r;
			*(opi+k) = tau1i + tau2i;

			*(opr+k+m) = tau1r - tau2r;
			*(opi+k+m) = tau1i - tau2i;			

			double temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;			

		}

	}
}

void inline radix2_dit_rec(dblitr opr, dblitr opi,dblitr ipr, dblitr ipi, int sgn,int N,int l) {

	if (N == 1) {
		*opr = *ipr;
		*opi = *ipi;
		//return;
	} else {
		int m = N / 2;

		radix2_dit_rec(opr,opi,ipr,ipi,sgn,m,l*2);
		radix2_dit_rec(opr+m,opi+m,ipr+l,ipi+l,sgn,m,l*2);

		double PI2 = 6.28318530717958647692528676655900577;
		double theta =  -1.0 * sgn * PI2/N;

		double S =  sin(theta);
		double C =  cos(theta);

		double wlr =  1.0;
		double wli =  0.0;

		double tau1r,tau1i,tau2r,tau2i;

		for (int k = 0; k < m; ++k) {
			tau1r = *(opr+k);
			tau1i = *(opi+k);

			tau2r = *(opr+k+m)*wlr - *(opi+k+m)*wli;
			tau2i = *(opi+k+m)*wlr + *(opr+k+m)*wli;

			*(opr+k) = tau1r + tau2r;
			*(opi+k) = tau1i + tau2i;

			*(opr+k+m) = tau1r - tau2r;
			*(opi+k+m) = tau1i - tau2i;

			double temp = wlr;
			wlr = C * wlr - S * wli;
			wli = S * temp + C * wli;

		}

	}
}

#endif 
