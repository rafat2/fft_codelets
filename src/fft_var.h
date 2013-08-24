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

#ifndef FFT_VAR_H
#define FFT_VAR__H

#include "codelets4.h"

template <typename T>
void inline fftct_radix2_dit(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int L;
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	bitreverse2(data);
	
	for (int i = 1; i < t+1; i++) {
		L = (int) pow(2.0, (double)i);
		fft_data<T> twi;
		twiddle_rcbs(twi,L);
		if (sgn == 1) {
			transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
		}
		radix2_dit(data,twi,i);
		
	}
	
	
	
}

template <typename T>
void inline fftct_radix2_dit_us(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int L;
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	bitreverse2(data);
	fft_data<T> twi;
	longvector(twi,N);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	for (int i = 1; i < t+1; i++) {		
		radix2_dit_us(data,twi,i);
	}
	
	
	
}

template <typename T>
void inline fftct_radix2_dit_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	bitreverse2(data);
	
	for (int i=1; i < t+1; i++) {
    radix2_dit_inplace(data,i,sgn);
    
    }
	
}

template <typename T>
void inline fftsh_radix2_dit_us(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	fft_data<T> twi;
	longvector(twi,N);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	for (int i = 1; i < t+1; i++) {		
		sh_radix2_dit_us(data,twi,i);
	}
	
}

template <typename T>
void inline fftsh_radix2_dit(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	
	for (int i=1; i < t+1; i++) {
    sh_radix2_dit(data,i,sgn);
    
    }
	
}

template <typename T>
void inline fftct_radix2_dif_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	
	for (int i=t; i > 0; i--) {
    radix2_dif_inplace(data,i,sgn);
    
    }
	
	bitreverse2(data);
	
}

template <typename T>
void inline fftsh_radix2_dif(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	
	for (int i=t; i > 0; i--)  {
    sh_radix2_dif(data,i,sgn);
    
    }
	
}

template <typename T>
void inline fftct_radix2_dif_us(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	
	fft_data<T> twi;
	longvector(twi,N);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	for (int i=t; i > 0; i--) {
    radix2_dif_us(data,twi,i);
    
    }
	
	bitreverse2(data);
	
}


template <typename T>
void inline fftct_radix4_dit_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(4.0,ceil(log10(static_cast<double>(N))/log10(4.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(4.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(4);
	}
	
	indrev4(data);
	//bitreverse2(data);
	
	for (int i=1; i < t+1; i++) {
    radix4_dit_inplace(data,i,sgn);
    
    }
	
}

template <typename T>
void inline fftct_radix4_dif_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(4.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(4.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(4);
	}
	
	//indrev(data,index);
	
	for (int i=t; i > 0; i--) {
    radix4_dif_inplace(data,i,sgn);
    
    }
	
	indrev4(data);
	
}

template <typename T>
void inline fftsh_radix4_dit(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(4.0,ceil(log10(static_cast<double>(N))/log10(4.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(4.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(4);
	}
	
	//indrev(data,index);
	
	for (int i=1; i < t+1; i++) {
    sh_radix4_dit(data,i,sgn);
    
    }
	
}

template <typename T>
void inline fftsh_radix4_dif(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(4.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(4.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(4);
	}
	
	//indrev(data,index);
	
	for (int i=t; i > 0; i--) {
    sh_radix4_dif(data,i,sgn);
    
    }
	
}

template <typename T>
void inline fft_split_radix(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int L;
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	indrev2(data);
	
	vector<int> bvec;
	bvector(bvec,t);
	T taur,taui;
	
	for (int k = 0; k < N/2; k++) {
		if (bvec[k] == 1) {
			int index = 2*k;
			int index1 = index+1;
			taur = data.re[index1];
			taui = data.im[index1];
			
			data.re[index1] = data.re[index] - taur; 
			data.im[index1] = data.im[index] - taui; 
			
			data.re[index] = data.re[index] + taur; 
			data.im[index] = data.im[index] + taui; 
		}
	}
	
	for (int i = 2; i < t+1; i++) {
		split_radix(data,bvec,i,sgn);
	}
	
}


void inline fft_split_radix_rec(fft_data<double> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((double) 0.0);
			data.im.push_back((double) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((double) 0.0);
			data.im.push_back((double) 0.0);
		}
		N = K;
	}
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	indrev2(data);
	
	
	split_radix_rec(data.re.begin(),data.im.begin(),sgn,N);
}


template <typename T>
void inline fft_split_radix_dit(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int L;
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	indrev2(data);
	
	T taur,taui;
	
	int c1 = 0;
	int r = 4;
	
	while (c1 < N - 1) {
		for  (int i = c1 ; i < N; i+=r) {
			int index = i;
			int index1 = i+1;
			taur = data.re[index1];
			taui = data.im[index1];
			
			data.re[index1] = data.re[index] - taur; 
			data.im[index1] = data.im[index] - taui; 
			
			data.re[index] = data.re[index] + taur; 
			data.im[index] = data.im[index] + taui; 
		}
		c1 = 2*r - 2;
		r = r*4;
	}
	
	for (int i = 2; i < t+1; i++) {
		split_radix_dit(data,i,sgn);
	}
	
}

template <typename T>
void inline fft_split_radix_dif(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	if (N != len) {
		for (unsigned int al = 0; al < N - len; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
	}
	
	unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
	
	if (N < K) {
		for (unsigned int al = 0; al < K - N; al++) {
			data.re.push_back((T) 0.0);
			data.im.push_back((T) 0.0);
		}
		N = K;
	}
	
	int L;
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(2);
	}
	
	//indrev(data,index);
	
	T taur,taui;
	
	int c1 = 0;
	int r = 4;
	
	for (int i = t; i > 1; --i) {
		split_radix_dif(data,i,sgn);
	}
	
	while (c1 < N - 1) {
		for  (int i = c1 ; i < N; i+=r) {
			int index = i;
			int index1 = i+1;
			taur = data.re[index1];
			taui = data.im[index1];
			
			data.re[index1] = data.re[index] - taur; 
			data.im[index1] = data.im[index] - taui; 
			
			data.re[index] = data.re[index] + taur; 
			data.im[index] = data.im[index] + taui; 
		}
		c1 = 2*r - 2;
		r = r*4;
	}
	
	indrev2(data);
	
}

template <typename T>
void inline fft_bluestein(fft_data<T> &data,int sgn) {
	unsigned int len = data.re.size(); 
	unsigned int K,M;
	
	K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(len))/log10(2.0)));
	
	if (K < 2 * len - 2) {
		M = K * 2;
	} else {
		M = K;
	}
	
	
	int t = (int) ceil(log10(static_cast<double>(M))/log10(2.0));
	
	// Initialize the Toeplitz matrix scalar coefficients h_0 to h_(len-1)
	fft_data<T> hk,hlt,yn,twi;
	yn.re.resize(M,(T) 0.0);
	yn.im.resize(M,(T) 0.0);
	
	bluestein_hl(hk,hlt,len,M);
	
	//fft_data<T> hk = hl;
	
	T scale = (T) 1.0/M;
	
	// Apply FFT to hl
	transform(hk.im.begin(), hk.im.end(), hk.im.begin(),bind1st(multiplies<T>(),(T) scale));
	transform(hk.re.begin(), hk.re.end(), hk.re.begin(),bind1st(multiplies<T>(),(T) scale));
	
	indrev2(hk);
	
	longvector(twi,M);
	
	transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	
	for (int i = 1; i < t+1; i++) {		
		radix2_dit_us(hk,twi,i);
	}
	
	if (sgn == 1) {
		for (unsigned int i = 0; i < len; i++) {
			yn.re[i] = data.re[i] * hlt.re[i] + data.im[i] * hlt.im[i];
			yn.im[i] = -data.re[i] * hlt.im[i] + data.im[i] * hlt.re[i];
		}
	} else {
		for (unsigned int i = 0; i < len; i++) {
			yn.re[i] = data.re[i] * hlt.re[i] - data.im[i] * hlt.im[i];
			yn.im[i] = data.re[i] * hlt.im[i] + data.im[i] * hlt.re[i];
		}
	}
	
	indrev2(yn);
	
	for (int i = 1; i < t+1; i++) {		
		radix2_dit_us(yn,twi,i);
	}
	
	
	if (sgn == 1) {
		for (unsigned int i = 0; i < M; i++) {
			T temp = (T) yn.re[i] * hk.re[i] - yn.im[i] * hk.im[i];
			yn.im[i] = yn.re[i] * hk.im[i] + yn.im[i] * hk.re[i];
			yn.re[i] = temp;
		} 
	} else {
		for (unsigned int i = 0; i < M; i++) {
			T temp = (T) yn.re[i] * hk.re[i] + yn.im[i] * hk.im[i];
			yn.im[i] = -yn.re[i] * hk.im[i] + yn.im[i] * hk.re[i];
			yn.re[i] = temp;
		} 
	
	}
	
	
	//IFFT
	
	transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	
	indrev2(yn);
	for (int i = 1; i < t+1; i++) {		
		radix2_dit_us(yn,twi,i);
	}
	
	//transform(yn.im.begin(), yn.im.end(), yn.im.begin(),bind1st(multiplies<T>(),(T) scale));
	//transform(yn.re.begin(), yn.re.end(), yn.re.begin(),bind1st(multiplies<T>(),(T) scale));

	if (sgn == 1) {
		for (unsigned int i = 0; i < len; i++) {
			data.re[i] = yn.re[i] * hlt.re[i] + yn.im[i] * hlt.im[i];
			data.im[i] = -yn.re[i] * hlt.im[i] + yn.im[i] * hlt.re[i];
		}
	} else {
		for (unsigned int i = 0; i < len; i++) {
			data.re[i] = yn.re[i] * hlt.re[i] - yn.im[i] * hlt.im[i];
			data.im[i] = yn.re[i] * hlt.im[i] + yn.im[i] * hlt.re[i];
		}
		
	}
	
	
}


template <typename T>
void inline fftct_radix2_dit_stride(fft_data<T> &data,int sgn, unsigned int stride, unsigned int N) {
	unsigned int len = data.re.size(); 
	unsigned int iter = len/N;
	int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vector<int> index;
	
	//indrev(data,index);
	//indrev2(data);
	fft_data<T> twi;
	longvector(twi,N);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	
	if (stride == 1) {
		for (int i = 0; i < iter; ++i) {
			indrev2_stride(data.re.begin()+i*N,data.im.begin()+i*N,1,N);
			for (int j = 1; j < t+1; ++j) {		
				radix2_dit_stride(data.re.begin()+i*N,data.im.begin()+i*N,twi,j,1,N);
			}
		}
	} else {
		for (int i = 0; i < iter; ++i) {
			indrev2_stride(data.re.begin()+i,data.im.begin()+i,stride,N);
			for (int j = 1; j < t+1; ++j) {		
				radix2_dit_stride(data.re.begin()+i,data.im.begin()+i,twi,j,stride,N);
			}
		}
	}
	
}

template <typename T>
void inline fft_prime(fft_data<T> &data,int sgn, int n1, int n2) {
	//int N = n1*n2;
	// IFFT needs work.
	fft_data<T> y;
	rmapT(data,y,n1,n2);
	
	switch(n1) {
		case 2 : radix2_wfta(y,1,n2,sgn);break;
		case 3 : radix3_wfta(y,1,n2,sgn);break;
		case 4 : radix4_wfta(y,1,n2,sgn);break;
		case 5 : radix5_wfta(y,1,n2,sgn);break;
		case 8 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 16 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 32 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 64 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 128 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 256 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 512 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 1024 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 2048 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 4096 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 8192 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 16384 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 32768 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
		case 65536 : fftct_radix2_dit_stride(y,sgn,1,n1);break;
	}
	
	switch(n2) {
		case 2 : radix2_wfta(y,n1,n1,sgn);break;
		case 3 : radix3_wfta(y,n1,n1,sgn);break;
		case 4 : radix4_wfta(y,n1,n1,sgn);break;
		case 5 : radix5_wfta(y,n1,n1,sgn);break;
		case 8 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 16 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 32 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 64 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 128 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 256 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 512 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 1024 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 2048 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 4096 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 8192 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 16384 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 32768 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
		case 65536 : fftct_radix2_dit_stride(y,sgn,n1,n2);break;
	}
	
	cmap(y,data,n1,n2);
	
}

template <typename T>
void inline fftct_radix3_dit_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(3.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(3);
	}
	
	indrev(data,index);
	//bitreverse2(data);
	fft_data<T> twi;
	
	twiddle(twi,N,3);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	for (int i=1; i < t+1; i++) {
    radix3_dit_inplace(data,twi,i,sgn);
    
    }
	
}

template <typename T>
void inline fftct_radix3_dif_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(3.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(3);
	}
	
	//indrev(data,index);
	fft_data<T> twi;
	
	twiddle(twi,N,3);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	for (int i=t; i > 0; i--) {
    radix3_dif_inplace(data,twi,i,sgn);
    
    }
	
	indrev(data,index);
	
}

template <typename T>
void inline fftsh_radix3_dit(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(3.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(3);
	}
	
	fft_data<T> twi;
	
	twiddle(twi,N,3);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	//indrev(data,index);
	//bitreverse2(data);
	
	for (int i=1; i < t+1; i++) {
    sh_radix3_dit(data,twi,i,sgn);
    
    }
	
}

template <typename T>
void inline fftsh_radix3_dif(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(3.0));

	//indrev(data,index);
	fft_data<T> twi;
	
	twiddle(twi,N,3);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	for (int i=t; i > 0; i--) {
    sh_radix3_dif(data,twi,i,sgn);
    
    }
	
}


template <typename T>
void inline fftct_radix5_dit_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(5.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(5);
	}
	
	indrev(data,index);
	//bitreverse2(data);
	fft_data<T> twi;
	
	twiddle(twi,N,5);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	for (int i=1; i < t+1; i++) {
    radix5_dit_inplace(data,twi,i,sgn);
    
    }
	
}

template <typename T>
void inline fftct_radix5_dif_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(5.0));
	vector<int> index;
	
	for (int i = 0; i < t; i++) {
		index.push_back(5);
	}
	
	//indrev(data,index);
	fft_data<T> twi;
	
	twiddle(twi,N,5);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	

	for (int i=t; i > 0; i--) {
    radix5_dif_inplace(data,twi,i,sgn);
    
    }
	
	indrev(data,index);
	
}

template <typename T>
void inline fftsh_radix5_dit(fft_data<T> &data,int sgn, unsigned int N) {
	unsigned int len = data.re.size(); 
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(5.0));
	fft_data<T> twi;
	
	twiddle(twi,N,5);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	//indrev(data,index);
	//bitreverse2(data);
	
	for (int i=1; i < t+1; i++) {
    sh_radix5_dit(data,twi,i,sgn);
    
    }
	
}

template <typename T>
void inline fftsh_radix5_dif(fft_data<T> &data,int sgn, unsigned int N) {
	//unsigned int len = data.re.size(); 
	
	int t = (int) ceil(log10(static_cast<double>(N))/log10(5.0));

	//indrev(data,index);
	fft_data<T> twi;
	
	twiddle(twi,N,5);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}
	
	for (int i=t; i > 0; i--) {
    sh_radix5_dif(data,twi,i,sgn);
    
    }
	
}

template <typename T>
void inline fftct_mixed_dit_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	/*
	 * Only works if N can be factored into integers 2,3,4 and 5
	 * For other lengths use M = findnext(N) to find the closest length
	 *  that can be factored into 2,3,4,5 and resize the input.
	 *  See the next few commented lines below.
	 */
	// For Example
	//unsigned int M = findnext(N);
	//if (M > N) {
	//	data.re.resize(M,(T) 0.0);
	//	data.im.resize(M,(T) 0.0);
	//}


	vector<int> index;
	factors(N,index);
	int t = index.size();
	indrevasym(data,index);
    int L = 1;
    fft_data<T> twi;

    twiddle(twi,N,2);
    if (sgn == 1) {
    	transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
    }

	for (int i=1; i < t+1; i++) {
		int p = index[i-1];
		L = L *p;
		switch(p) {
			case 2 : radix2_mixed_dit_inplace(data,twi,L,sgn);break;
			case 3 : radix3_mixed_dit_inplace(data,twi,L,sgn);break;
			case 4 : radix4_mixed_dit_inplace(data,twi,L,sgn);break;
			case 5 : radix5_mixed_dit_inplace(data,twi,L,sgn);break;
			case 7 : radix7_mixed_dit_inplace(data,twi,L,sgn);break;
			case 8 : radix8_mixed_dit_inplace(data,twi,L,sgn);break;
			default : radixN_mixed_dit_inplace(data,twi,L,p,sgn);break;

		}

    }

}

template <typename T>
void inline fftct_mixed_dif_inplace(fft_data<T> &data,int sgn, unsigned int N) {
	/*
	 * Only works if N can be factored into integers 2,3,4 and 5
	 * For other lengths use M = findnext(N) to find the closest length
	 *  that can be factored into 2,3,4,5 and resize the input.
	 *  See the next few commented lines below.
	 */
	// For Example
	//unsigned int M = findnext(N);
	//if (M > N) {
	//	data.re.resize(M,(T) 0.0);
	//	data.im.resize(M,(T) 0.0);
	//}


	vector<int> index;
	factors(N,index);
	int t = index.size();
    int L = N;
    fft_data<T> twi;

    twiddle(twi,N,2);
    if (sgn == 1) {
        transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
    }

	for (int i=t; i > 0; --i) {
		int p = index[t-i];
		switch(p) {
			case 2 : radix2_mixed_dif_inplace(data,twi,L,sgn);break;
			case 3 : radix3_mixed_dif_inplace(data,twi,L,sgn);break;
			case 4 : radix4_mixed_dif_inplace(data,twi,L,sgn);break;
			case 5 : radix5_mixed_dif_inplace(data,twi,L,sgn);break;
			case 7 : radix7_mixed_dif_inplace(data,twi,L,sgn);break;
			case 8 : radix8_mixed_dif_inplace(data,twi,L,sgn);break;
			default : radixN_mixed_dif_inplace(data,twi,L,p,sgn);break;

		}
		L = L / p;

    }

	indrevasym(data,index);


}


template <typename T>
void inline fftsh_mixed_dit(fft_data<T> &data,int sgn, unsigned int N) {
	/*
	 * Only works if N can be factored into integers 2,3,4 and 5
	 * For other lengths use M = findnext(N) to find the closest length
	 *  that can be factored into 2,3,4,5 and resize the input.
	 *  See the next few commented lines below.
	 */
	// For Example
	//unsigned int M = findnext(N);
	//if (M > N) {
	//	data.re.resize(M,(T) 0.0);
	//	data.im.resize(M,(T) 0.0);
	//}


	vector<int> index;
	factors(N,index);
	int t = index.size();
    int L = 1;
    fft_data<T> twi;

    twiddle(twi,N,2);
    if (sgn == 1) {
        transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
    }

	for (int i=1; i < t+1; i++) {
		int p = index[i-1];
		L = L *p;
		switch(p) {
			case 2 : sh_mixed_radix2_dit(data,twi,L,sgn);break;
			case 3 : sh_mixed_radix3_dit(data,twi,L,sgn);break;
			case 4 : sh_mixed_radix4_dit(data,twi,L,sgn);break;
			case 5 : sh_mixed_radix5_dit(data,twi,L,sgn);break;
			case 7 : sh_mixed_radix7_dit(data,twi,L,sgn);break;
			case 8 : sh_mixed_radix8_dit(data,twi,L,sgn);break;
			default : sh_mixed_radixN_dit(data,twi,L,p,sgn);break;

		}

    }

}

template <typename T>
void inline fftsh_mixed_dif(fft_data<T> &data,int sgn, unsigned int N) {
	/*
	 * Only works if N can be factored into integers 2,3,4 and 5
	 * For other lengths use M = findnext(N) to find the closest length
	 *  that can be factored into 2,3,4,5 and resize the input.
	 *  See the next few commented lines below.
	 */
	// For Example
	//unsigned int M = findnext(N);
	//if (M > N) {
	//	data.re.resize(M,(T) 0.0);
	//	data.im.resize(M,(T) 0.0);
	//}


	vector<int> index;
	factors(N,index);
	int t = index.size();
    int L = N;
    fft_data<T> twi;

    twiddle(twi,N,2);
    if (sgn == 1) {
        transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
    }

	for (int i=t; i > 0; --i) {
		int p = index[t-i];
		switch(p) {
			case 2 : sh_mixed_radix2_dif(data,twi,L,sgn);break;
			case 3 : sh_mixed_radix3_dif(data,twi,L,sgn);break;
			case 4 : sh_mixed_radix4_dif(data,twi,L,sgn);break;
			case 5 : sh_mixed_radix5_dif(data,twi,L,sgn);break;
			case 7 : sh_mixed_radix7_dif(data,twi,L,sgn);break;
			case 8 : sh_mixed_radix8_dif(data,twi,L,sgn);break;
			default : sh_mixed_radixN_dif(data,twi,L,p,sgn);break;

		}
		L = L / p;

    }


}

template <typename T>
void inline fftct_radix7_dit_inplace(fft_data<T> &data,int sgn, unsigned int N) {

	int t = (int) ceil(log10(static_cast<double>(N))/log10(7.0));
	vector<int> index;

	for (int i = 0; i < t; i++) {
		index.push_back(7);
	}
	int L = 1;

	indrev(data,index);
	//bitreverse2(data);
	fft_data<T> twi;

	twiddle(twi,N,7);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}


	for (int i=1; i < t+1; i++) {
		int p = index[i-1];
		L = L *p;
		radix7_mixed_dit_inplace(data,twi,L,sgn);

    }

}

template <typename T>
void inline fftct_radix7_dif_inplace(fft_data<T> &data,int sgn, unsigned int N) {

	int t = (int) ceil(log10(static_cast<double>(N))/log10(7.0));
	vector<int> index;

	for (int i = 0; i < t; i++) {
		index.push_back(7);
	}

	//indrev(data,index);

	int L = N;
	fft_data<T> twi;

	twiddle(twi,N,7);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}

	for (int i=t; i > 0; i--) {
		int p = index[t-i];
		radix7_mixed_dif_inplace(data,twi,L,sgn);
		L = L / p;

    }

	indrev(data,index);

}

template <typename T>
void inline fftsh_radix7_dit(fft_data<T> &data,int sgn, unsigned int N) {

	int t = (int) ceil(log10(static_cast<double>(N))/log10(7.0));
	vector<int> index;

	for (int i = 0; i < t; i++) {
		index.push_back(7);
	}
	int L = 1;
	fft_data<T> twi;

	twiddle(twi,N,7);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}

	for (int i=1; i < t+1; i++) {
		int p = index[i-1];
		L = L *p;
		sh_mixed_radix7_dit(data,twi,L,sgn);

    }

}

template <typename T>
void inline fftsh_radix7_dif(fft_data<T> &data,int sgn, unsigned int N) {

	int t = (int) ceil(log10(static_cast<double>(N))/log10(7.0));
	vector<int> index;

	for (int i = 0; i < t; i++) {
		index.push_back(7);
	}

	//indrev(data,index);

	int L = N;
	fft_data<T> twi;

	twiddle(twi,N,7);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}


	for (int i=t; i > 0; i--) {
		int p = index[t-i];
		sh_mixed_radix7_dif(data,twi,L,sgn);
		L = L / p;

    }

}

template <typename T>
void inline fftct_radix8_dit_inplace(fft_data<T> &data,int sgn, unsigned int N) {

	vector<int> index;

	factorize(N,8,index);
	int t;
	t = (int) index.size();

	int L = 1;

	indrev8(data);
	//bitreverse2(data);
	fft_data<T> twi;

	twiddle(twi,N,8);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}


	for (int i=1; i < t+1; i++) {
		int p = index[i-1];
		L = L *p;
		radix8_mixed_dit_inplace(data,twi,L,sgn);

    }

}

template <typename T>
void inline fftct_radix8_dif_inplace(fft_data<T> &data,int sgn, unsigned int N) {

	vector<int> index;

	factorize(N,8,index);
	int t;
	t = (int) index.size();

	int L = N;

	//bitreverse2(data);
	fft_data<T> twi;

	twiddle(twi,N,8);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}


	for (int i=t; i > 0; i--) {
		int p = index[t-i];
		radix8_mixed_dif_inplace(data,twi,L,sgn);
		L = L / p;

    }

	indrev8(data);

}

template <typename T>
void inline fftsh_radix8_dit(fft_data<T> &data,int sgn, unsigned int N) {

	vector<int> index;

	factorize(N,8,index);
	int t;
	t = (int) index.size();

	int L = 1;

	//indrev8(data);
	//bitreverse2(data);
	fft_data<T> twi;

	twiddle(twi,N,8);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}

	for (int i=1; i < t+1; i++) {
		int p = index[i-1];
		L = L *p;
		sh_mixed_radix8_dit(data,twi,L,sgn);

    }

}

template <typename T>
void inline fftsh_radix8_dif(fft_data<T> &data,int sgn, unsigned int N) {

	vector<int> index;

	factorize(N,8,index);
	int t;
	t = (int) index.size();

	int L = N;

	//bitreverse2(data);
	fft_data<T> twi;

	twiddle(twi,N,8);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}

	for (int i=t; i > 0; i--) {
		int p = index[t-i];
		sh_mixed_radix8_dif(data,twi,L,sgn);
		L = L / p;

    }


}

template <typename T>
void inline fftct_radixN_dit_inplace(fft_data<T> &data,int sgn, unsigned int N, unsigned int radix) {

	vector<int> index;

	factorize(N,radix,index);
	int t;
	t = (int) index.size();
	int L = 1;

	indrev(data,index);
	//bitreverse2(data);
	fft_data<T> twi;

	twiddle(twi,N,radix);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}


	for (int i=1; i < t+1; i++) {
		int p = index[i-1];
		L = L *p;
		radixN_mixed_dit_inplace(data,twi,L,radix,sgn);

    }

}


template <typename T>
void inline fftct_radixN_dif_inplace(fft_data<T> &data,int sgn, unsigned int N, unsigned int radix) {

	vector<int> index;

	factorize(N,radix,index);
	int t;
	t = (int) index.size();
	int L = N;

	//bitreverse2(data);
	fft_data<T> twi;

	twiddle(twi,N,radix);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}



	for (int i=t; i > 0; i--) {
		int p = index[i-1];
		radixN_mixed_dif_inplace(data,twi,L,radix,sgn);
		L = L/p;

    }
	indrev(data,index);

}

template <typename T>
void inline fftsh_radixN_dit(fft_data<T> &data,int sgn, unsigned int N, unsigned int radix) {

	vector<int> index;

	factorize(N,radix,index);
	int t;
	t = (int) index.size();
	int L = 1;

	fft_data<T> twi;

	twiddle(twi,N,radix);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}


	for (int i=1; i < t+1; i++) {
		int p = index[i-1];
		L = L *p;
		sh_mixed_radixN_dit(data,i,L,radix,sgn);

    }

}

template <typename T>
void inline fftsh_radixN_dif(fft_data<T> &data,int sgn, unsigned int N, unsigned int radix) {

	vector<int> index;

	factorize(N,radix,index);
	int t;
	t = (int) index.size();
	int L = N;

	//bitreverse2(data);
	fft_data<T> twi;

	twiddle(twi,N,radix);
	if (sgn == 1) {
		transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));
	}




	for (int i=t; i > 0; i--) {
		int p = index[i-1];
		sh_mixed_radixN_dif(data,i,L,radix,sgn);
		L = L/p;

    }

}


template <typename T>
void inline fft_bluestein2(fft_data<T> &data,int sgn) {
	unsigned int len = data.re.size();
	unsigned int K,M;

	K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(len))/log10(2.0)));

	if (K < 2 * len - 2) {
		M = K * 2;
	} else {
		M = K;
	}


	int t = (int) ceil(log10(static_cast<double>(M))/log10(2.0));

	// Initialize the Toeplitz matrix scalar coefficients h_0 to h_(len-1)
	fft_data<T> hk,hlt,yn,twi;
	yn.re.resize(M,(T) 0.0);
	yn.im.resize(M,(T) 0.0);

	bluestein_hl(hk,hlt,len,M);

	//fft_data<T> hk = hl;

	T scale = (T) 1.0/M;

	// Apply FFT to hl

	for (unsigned int ii = 0; ii < hk.im.size(); ++ii) {
		hk.im[ii] = hk.im[ii] * scale;
		hk.re[ii] = hk.re[ii] * scale;
	}


	indrev2(hk);

	nlongvector(twi,M);

	//transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));

	for (int i = 1; i < t+1; i++) {
		radix2_dit_us(hk,twi,i);
	}

	if (sgn == 1) {
		for (unsigned int i = 0; i < len; i++) {
			yn.re[i] = data.re[i] * hlt.re[i] + data.im[i] * hlt.im[i];
			yn.im[i] = -data.re[i] * hlt.im[i] + data.im[i] * hlt.re[i];
		}
	} else {
		for (unsigned int i = 0; i < len; i++) {
			yn.re[i] = data.re[i] * hlt.re[i] - data.im[i] * hlt.im[i];
			yn.im[i] = data.re[i] * hlt.im[i] + data.im[i] * hlt.re[i];
		}
	}

	indrev2(yn);

	for (int i = 1; i < t+1; i++) {
		radix2_dit_us(yn,twi,i);
	}


	if (sgn == 1) {
		for (unsigned int i = 0; i < M; i++) {
			T temp = (T) yn.re[i] * hk.re[i] - yn.im[i] * hk.im[i];
			yn.im[i] = yn.re[i] * hk.im[i] + yn.im[i] * hk.re[i];
			yn.re[i] = temp;
		}
	} else {
		for (unsigned int i = 0; i < M; i++) {
			T temp = (T) yn.re[i] * hk.re[i] + yn.im[i] * hk.im[i];
			yn.im[i] = -yn.re[i] * hk.im[i] + yn.im[i] * hk.re[i];
			yn.re[i] = temp;
		}

	}


	//IFFT

	for (unsigned int ii = 0; ii < twi.im.size(); ++ii) {
		twi.im[ii] = -twi.im[ii];
	}
	//transform(twi.im.begin(), twi.im.end(), twi.im.begin(),bind1st(multiplies<T>(),(T) -1.0));

	indrev2(yn);
	for (int i = 1; i < t+1; i++) {
		radix2_dit_us(yn,twi,i);
	}

	//transform(yn.im.begin(), yn.im.end(), yn.im.begin(),bind1st(multiplies<T>(),(T) scale));
	//transform(yn.re.begin(), yn.re.end(), yn.re.begin(),bind1st(multiplies<T>(),(T) scale));

	if (sgn == 1) {
		for (unsigned int i = 0; i < len; i++) {
			data.re[i] = yn.re[i] * hlt.re[i] + yn.im[i] * hlt.im[i];
			data.im[i] = -yn.re[i] * hlt.im[i] + yn.im[i] * hlt.re[i];
		}
	} else {
		for (unsigned int i = 0; i < len; i++) {
			data.re[i] = yn.re[i] * hlt.re[i] - yn.im[i] * hlt.im[i];
			data.im[i] = yn.re[i] * hlt.im[i] + yn.im[i] * hlt.re[i];
		}

	}


}

template <typename T>
void inline fftct_radix2_dit_rec(fft_data<T> &data,fft_data<T> &y,int sgn, unsigned int N) {


	y = data;

	radix2_dit_rec(y.re.begin(),y.im.begin(),data.re.begin(),data.im.begin(),sgn,N,1);

	//data = y;



}

template <typename T>
void inline fftct_radix3_dit_rec(fft_data<T> &data,fft_data<T> &y,int sgn, unsigned int N) {


	y = data;

	radix3_dit_rec(y.re.begin(),y.im.begin(),data.re.begin(),data.im.begin(),sgn,N,1);

	//data = y;

}

template <typename T>
void inline fftct_radix4_dit_rec(fft_data<T> &data,fft_data<T> &y,int sgn, unsigned int N) {


	y = data;

	radix4_dit_rec(y.re.begin(),y.im.begin(),data.re.begin(),data.im.begin(),sgn,N,1);

	//data = y;

}

template <typename T>
void inline fftct_radix5_dit_rec(fft_data<T> &data,fft_data<T> &y,int sgn, unsigned int N) {


	y = data;

	radix5_dit_rec(y.re.begin(),y.im.begin(),data.re.begin(),data.im.begin(),sgn,N,1);

	//data = y;

}

template <typename T>
void inline fftct_radix7_dit_rec(fft_data<T> &data,fft_data<T> &y,int sgn, unsigned int N) {


	y = data;

	radix7_dit_rec(y.re.begin(),y.im.begin(),data.re.begin(),data.im.begin(),sgn,N,1);

	//data = y;

}

template <typename T>
void inline fftct_radix8_dit_rec(fft_data<T> &data,fft_data<T> &y,int sgn, unsigned int N) {


	y = data;

	radix8_dit_rec(y.re.begin(),y.im.begin(),data.re.begin(),data.im.begin(),sgn,N,1);

	//data = y;

}



#endif
