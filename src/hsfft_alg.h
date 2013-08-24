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

#ifndef HSFFT_ALG_H
#define HSFFT_ALG_H
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

//#define PI 3.141592653589793238462643383279502884197169399375105820974944;

template <typename T>
class fft_data{
public:	
	vector<T> re;
	vector<T> im;

};

template <typename T>
void inline subvec_scaling(fft_data<T> &vec,int N) {
	int q = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	T PI = (T) 3.141592653589793238462643383279502884197169399375105820974944;
	//cout << q << endl;
	vec.re.resize(N/2,(T) 0.0);
	vec.im.resize(N/2,(T) 0.0);
	vec.re[0] = (T) 1.0;
	vec.im[0] = (T) 0.0;
	T theta,wr,wi;
	int th1 = 2;
	int L = N;
	
	
	for (int j = 1; j < q; ++j) {
		theta = (T) th1 * PI / L;
		wr = cos(theta);
		wi = sin(theta);
		int k = 0;
		for (int i = th1 / 2;i < th1; ++i ) {
			vec.re[i] = wr * vec.re[k] - wi * vec.im[k];
			vec.im[i] = wr * vec.im[k] + wi * vec.re[k];
			++k;
		}
		
		th1*=2;
	}
	
}

template <typename T>
void inline rec_bisection(fft_data<T> &vec,int N) {
	int q = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) PI2/N;
	vec.re.resize(N,(T) 0.0);
	vec.im.resize(N,(T) 0.0);
	vec.re[0] = (T) 1.0;
	vec.im[0] = (T) 0.0;
	
	int p = 1;
	T h;
	
	for (int k=0; k < q; ++k) {
		vec.re[p] = cos(p * theta);
		vec.im[p] = sin(p * theta);
		p = p * 2;
	}
	int i;
	
	for (int j = 1; j < q-1; ++j) {
		p = (int) pow(2.0,(double) q-j-2);
		h = 1.0 / (2.0 * vec.re[p]);
		int lam = (int) pow(2.0,(double) j) - 1;
		for (int k = 0; k < lam; ++k) {
			i = p * (3 + 2*k);
			vec.re[i] = h * (vec.re[i - p] + vec.re[i + p]);
			vec.im[i] = h * (vec.im[i - p] + vec.im[i + p]);
		}
		
	}
	vec.re.resize(N/2);
	vec.im.resize(N/2);
	
}

template <typename T>
void inline twiddle_rcbs(fft_data<T> &vec,int N) {
	int q = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	vec.re.resize(N/2,(T) 0.0);
	vec.im.resize(N/2,(T) 0.0);
	vec.re[0] = (T) 1.0;
	vec.im[0] = (T) 0.0;
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta2 = (T) PI2/N;
	
	if (N < 16) {
		for (int K = 1; K < N/2; K++) {
			vec.re[K] = (T) cos(theta2 * K);
			vec.im[K] = (T) sin(theta2 * K);
		}
	} else {
		int t = q - 2;
		int p = 1;
		T h;
		
		for (int k=0; k < t; ++k) {
			vec.re[p] = cos(p * theta2);
			vec.im[p] = sin(p * theta2);
			p = p * 2;
		}
		int i;
		
		for (int j = 1; j < t-1; ++j) {
			p = (int) pow(2.0,(double) t-j-2);
			h = 1.0 / (2.0 * vec.re[p]);
			int lam = (int) pow(2.0,(double) j) - 1;
			for (int k = 0; k < lam; ++k) {
				i = p * (3 + 2*k);
				vec.re[i] = h * (vec.re[i - p] + vec.re[i + p]);
				vec.im[i] = h * (vec.im[i - p] + vec.im[i + p]);
			}
			
		}
		int N8 = N/8;
		vec.re[N8] = (T) sqrt(2.0) / 2.0;
		vec.im[N8] = vec.re[N/8];
		
		for (int k = 1; k <= N8 - 1; k++) {
			vec.re[N8 + k] = vec.im[N8 - k];
			vec.im[N8 + k] = vec.re[N8 - k];
		}
		int N4 = N/4;
		vec.im[N4] = (T) 1.0;
		for (int k=1; k <= N4 - 1; k++) {
			vec.re[N4 + k] = -1.0 * vec.re[N4 - k];
			vec.im[N4 + k] = vec.im[N4 - k];
		}
	}

	
}



template <typename T>
void inline twiddle_svsc(fft_data<T> &vec,int N) {
	int q = (int) ceil(log10(static_cast<double>(N))/log10(2.0));
	T PI = (T) 3.141592653589793238462643383279502884197169399375105820974944;
	T theta,wr,wi;
	vec.re.resize(N/2,(T) 0.0);
	vec.im.resize(N/2,(T) 0.0);
	vec.re[0] = (T) 1.0;
	vec.im[0] = (T) 0.0;
	int th1 = 2;
	int L = N;
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta2 = (T) PI2/N;
	
	if (N < 16) {
		for (int K = 1; K < N/2; K++) {
			vec.re[K] = (T) cos(theta2 * K);
			vec.im[K] = (T) sin(theta2 * K);
		}
	} else {
		for (int j = 1; j < q-2; ++j) {
			theta = (T) th1 * PI / L;
			wr = cos(theta);
			wi = sin(theta);
			int k = 0;
			for (int i = th1 / 2;i < th1; ++i ) {
				vec.re[i] = wr * vec.re[k] - wi * vec.im[k];
				vec.im[i] = wr * vec.im[k] + wi * vec.re[k];
				++k;
			}
			
			th1*=2;
		}
		int N8 = N/8;
		vec.re[N8] = (T) sqrt(2.0) / 2.0;
		vec.im[N8] = vec.re[N/8];
		
		for (int k = 1; k <= N8 - 1; k++) {
			vec.re[N8 + k] = vec.im[N8 - k];
			vec.im[N8 + k] = vec.re[N8 - k];
		}
		int N4 = N/4;
		vec.im[N4] = (T) 1.0;
		for (int k=1; k <= N4 - 1; k++) {
			vec.re[N4 + k] = -1.0 * vec.re[N4 - k];
			vec.im[N4 + k] = vec.im[N4 - k];
		}
	}

	
}


template <typename T>
void inline twiddle(fft_data<T> &vec,int N){
	// Calculates twiddle factors for radix-2
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) PI2/N;
	T S = (T) sin(theta);
	T C = (T) cos(theta);
	vec.re.resize(N/2,(T) 0.0);
	vec.im.resize(N/2,(T) 0.0);
	vec.re[0] = (T) 1.0;
	
	if (N < 16) {
		for (int K = 1; K < N/2; K++) {
			vec.re[K] = (T) cos(theta * K);
			vec.im[K] = (T) sin(theta * K);
		}
	} else {
		for (int k = 0; k <= (N/8) - 2; k++) {
			vec.re[k+1] = (T) C * vec.re[k] - S * vec.im[k];
			vec.im[k+1] = (T) S * vec.re[k] + C * vec.im[k];
		}
		int N8 = N/8;
		vec.re[N8] = (T) sqrt(2.0) / 2.0;
		vec.im[N8] = vec.re[N/8];
		
		for (int k = 1; k <= N8 - 1; k++) {
			vec.re[N8 + k] = vec.im[N8 - k];
			vec.im[N8 + k] = vec.re[N8 - k];
		}
		int N4 = N/4;
		vec.im[N4] = (T) 1.0;
		for (int k=1; k <= N4 - 1; k++) {
			vec.re[N4 + k] = -1.0 * vec.re[N4 - k];
			vec.im[N4 + k] = vec.im[N4 - k];
		}
	}

	
		
}

template <typename T>
void inline twiddle(fft_data<T> &vec,int N,int radix){
	// Calculates twiddle factors for radix-2
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) PI2/N;
	vec.re.resize(N/radix,(T) 0.0);
	vec.im.resize(N/radix,(T) 0.0);
	vec.re[0] = (T) 1.0;
	
	for (int K = 1; K < N/radix; K++) {
		vec.re[K] = (T) cos(theta * K);
		vec.im[K] = (T) sin(theta * K);
	}

	
}

template <typename T>
void inline ntwiddle(fft_data<T> &vec,int N,int radix){
	// Calculates twiddle factors for radix-2
	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) PI2/N;
	vec.re.resize(N/radix,(T) 0.0);
	vec.im.resize(N/radix,(T) 0.0);
	vec.re[0] = (T) 1.0;

	for (int K = 1; K < N/radix; K++) {
		vec.re[K] = (T) cos(theta * K);
		vec.im[K] = (T) -sin(theta * K);
	}


}

template <typename T>
void indrev(fft_data<T> &sig, vector<int> &index) {

  int t = index.size();
  int n = sig.re.size();
  int j,m,s,alpha;
  T tempr, tempi;
  for(int k = 0; k < n; k++){
    m = k;
    j = 0;
    for (int q = t; q > 0; q--){
      s = m/index[q-1];
      alpha = m - s*index[q-1];
      j = j * index[q-1] + alpha;
      m = s;
    }
    if (j > k){
      tempr = sig.re[j];
	  tempi = sig.im[j];
      sig.re[j] = sig.re[k];
	  sig.im[j] = sig.im[k];
      sig.re[k] = tempr;
	  sig.im[k] = tempi;
    }
  }
}

template <typename T>
void indrevasym(fft_data<T> &sig, vector<int> &index) {

  int t = index.size();
  int n = sig.re.size();
  int j,m,s,alpha;
  fft_data<T> sig2 = sig;
  for(int k = 0; k < n; k++){
    m = k;
    j = 0;
    for (int q = t; q > 0; q--){
      s = m/index[q-1];
      alpha = m - s*index[q-1];
      j = j * index[q-1] + alpha;
      m = s;
    }
	sig.re[j] = sig2.re[k];
	sig.im[j] = sig2.im[k];
  }
}


template <typename T>
void indrev4(fft_data<T> &sig) {

  int n = sig.re.size();
  int j = 0;
  int n2 = n>>2;
  T tempr, tempi;
  for(int i = 0; i < n-1; i++){
    if (j > i){
      tempr = sig.re[j];
	  tempi = sig.im[j];
      sig.re[j] = sig.re[i];
	  sig.im[j] = sig.im[i];
      sig.re[i] = tempr;
	  sig.im[i] = tempi;
    }
	int k = n2;
	while (j >= k*3) {
		j-=k*3;
		k>>=2;
	}
	j+=k;
  }
}

template <typename T>
void indrev2(fft_data<T> &sig) {

  int n = sig.re.size();
  int j = 0;
  int n2 = n>>1;
  T tempr, tempi;
  for(int i = 0; i < n-1; i++){
    if (j > i){
      tempr = sig.re[j];
	  tempi = sig.im[j];
      sig.re[j] = sig.re[i];
	  sig.im[j] = sig.im[i];
      sig.re[i] = tempr;
	  sig.im[i] = tempi;
    }
	int k = n2;
	while (j >= k) {
		j-=k;
		k>>=1;
	}
	j+=k;
  }
}

template <typename T>
void indrev8(fft_data<T> &sig) {

  int n = sig.re.size();
  int j = 0;
  int n2 = n>>3;
  T tempr, tempi;
  for(int i = 0; i < n-1; i++){
    if (j > i){
      tempr = sig.re[j];
	  tempi = sig.im[j];
      sig.re[j] = sig.re[i];
	  sig.im[j] = sig.im[i];
      sig.re[i] = tempr;
	  sig.im[i] = tempi;
    }
	int k = n2;
	while (j >= k*7) {
		j-=k*7;
		k>>=3;
	}
	j+=k;
  }
}

void inline indrev2_stride(vector<double>::iterator sigr, vector<double>::iterator sigi, int stride, int n) {

  //int n = sig.re.size();
  int j = 0;
  int jj,ii;
  int n2 = n>>1;
  double tempr, tempi;
  for(int i = 0; i < n-1; i++){
    if (j > i){
	  jj = j * stride;
	  ii = i * stride;
      tempr = *(sigr+jj);
	  tempi = *(sigi+jj);
      *(sigr+jj) = *(sigr+ii);
	  *(sigi+jj) = *(sigi+ii);
      *(sigr+ii) = tempr;
	  *(sigi+ii) = tempi;
    }
	int k = n2;
	while (j >= k) {
		j-=k;
		k>>=1;
	}
	j+=k;
  }
}

void inline indrev2_stride(vector<float>::iterator sigr, vector<float>::iterator sigi, int stride, int n) {

  //int n = sig.re.size();
  int j = 0;
  int jj,ii;
  int n2 = n>>1;
  float tempr, tempi;
  for(int i = 0; i < n-1; i++){
    if (j > i){
	  jj = j * stride;
	  ii = i * stride;
      tempr = *(sigr+jj);
	  tempi = *(sigi+jj);
      *(sigr+jj) = *(sigr+ii);
	  *(sigi+jj) = *(sigi+ii);
      *(sigr+ii) = tempr;
	  *(sigi+ii) = tempi;
    }
	int k = n2;
	while (j >= k) {
		j-=k;
		k>>=1;
	}
	j+=k;
  }
}


template <typename T>
void inline bitreverse2(fft_data<T> &sig) {
  unsigned int len = sig.re.size();
  unsigned int N = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(len))/log10(2.0)));
  unsigned int rev = 0;
// Processing Input Data
  for (unsigned int iter = 0; iter < N; ++iter)
    {
      if (rev > iter)
	{
// Replacing current values with reversed values

	  T tempr = sig.re[rev];
	  T tempi = sig.im[rev];
	  sig.re[rev] = sig.re[iter];
	  sig.im[rev] = sig.im[iter];
	  sig.re[iter] = tempr;
	  sig.im[iter] = tempi;

	}
// Using filter "filt" such that the value of reverse changes with each iteration
      unsigned int filt = N;
      while (rev & (filt >>= 1)) {
	rev &= ~filt;
      }
      rev |= filt;
    }
  

}


template <typename T>
void  inline downsamp(fft_data<T> &sig, int M, fft_data<T> &sig_d){
	int len = sig.re.size();
	T tempr,tempi;
	double len_n = ceil( (double) len / (double) M);
	for (int i = 0; i < (int) len_n; i++) {
		tempr = sig.re[i*M];
		tempi = sig.im[i*M];
		
		sig_d.re.push_back(tempr);
		sig_d.im.push_back(tempi);
	}
}

template <typename T>
void  inline evenodd(fft_data<T> &sig, fft_data<T> &sig_e, fft_data<T> &sig_o, int len){
	//int len = sig.re.size();

	sig_e.re.reserve(len);
	sig_e.im.reserve(len);
	sig_o.re.reserve(len);
	sig_o.im.reserve(len);
	T tempr,tempi;
	//double len_n = ceil( (double) len / (double) M);
	for (int i = 0; i < len; i++) {
		tempr = sig.re[2*i];
		tempi = sig.im[2*i];

		sig_e.re.push_back(tempr);
		sig_e.im.push_back(tempi);

		tempr = sig.re[2*i + 1];
		tempi = sig.im[2*i + 1];

		sig_o.re.push_back(tempr);
		sig_o.im.push_back(tempi);
	}
}

template <typename T>
void inline split(vector<T> &sig, vector<T> &even, vector<T> &odd)  {

	for (int i=0; i < (int) sig.size(); i++) {
		if (i%2 == 0) {
			even.push_back(sig[i]);
		} else {
			odd.push_back(sig[i]);
		}
	}
	
	
}

template <typename T>
void inline longvector(fft_data<T> &longvec,int N) {
	fft_data<T> twi;
	twiddle(twi,N);
	int N2 = N/2;
	longvec.re.push_back((T) 1.0);
	longvec.im.push_back((T) 0.0);
	int tx = (int) ceil(log10(static_cast<double>(N2))/log10(2.0));
	for (int i = 1; i < tx; i++) {
		int po =(int) pow(2.0, (double) tx - i);
		fft_data<T> sig2;
		downsamp(twi,po,sig2);
		longvec.re.insert(longvec.re.end(),sig2.re.begin(),sig2.re.end());
		longvec.im.insert(longvec.im.end(),sig2.im.begin(),sig2.im.end());
	}
	longvec.re.insert(longvec.re.end(),twi.re.begin(),twi.re.end());
	longvec.im.insert(longvec.im.end(),twi.im.begin(),twi.im.end());
}

template <typename T>
void inline nlongvector(fft_data<T> &longvec,int N) {
	fft_data<T> twi;
	ntwiddle(twi,N,2);
	int N2 = N/2;
	longvec.re.push_back((T) 1.0);
	longvec.im.push_back((T) 0.0);
	int tx = (int) ceil(log10(static_cast<double>(N2))/log10(2.0));
	for (int i = 1; i < tx; i++) {
		int po =(int) pow(2.0, (double) tx - i);
		fft_data<T> sig2;
		downsamp(twi,po,sig2);
		longvec.re.insert(longvec.re.end(),sig2.re.begin(),sig2.re.end());
		longvec.im.insert(longvec.im.end(),sig2.im.begin(),sig2.im.end());
	}
	longvec.re.insert(longvec.re.end(),twi.re.begin(),twi.re.end());
	longvec.im.insert(longvec.im.end(),twi.im.begin(),twi.im.end());
}

template <typename T>
void inline vecmult(vector<T> &sig, double x) {
	
	for (int i=0; i < (int) sig.size(); i++ ) {
		sig[i]= (T) (x*sig[i]);
	}
	
}

template <typename T>
void inline merge(vector<T> &sig, vector<T> &even, vector<T> &odd) {
	
	int N = even.size() + odd.size();
	
	for (int i=0; i < N; i++) {
		if (i%2 == 0) {
			sig.push_back(even[i/2]);
		} else {
			sig.push_back(odd[i/2]);
		}
	}
}

template <typename T> 
void inline transpose(vector<T> &sig, int rows, int cols,vector<T> &col) {
	//Replace
	for (int i=0; i < cols; i++) {
		for (int j=0; j < rows; j++) {
			col.push_back(sig[i+j*cols]);
		}
	}
	
}

void inline bvector(vector<int> &beta, int t) {
	
	//beta.resize(N/2,0);
	int N = (int) pow(2.0,t);
	beta.reserve(N/2);
	beta.push_back(1);
	for (int q = 0; q < t-1; q++) {
		beta.insert(beta.end(),beta.begin(),beta.end());
		int s = (int) pow(2.0,(double)q+1) - 1;
		if (beta[s] == 0) {
			beta[s] = 1;
		} else {
			beta[s] = 0;
		}
	}
	
}

void inline idempotent(int N1, int N2, int &e1, int &e2) {
	int N = N1*N2;
	for (int i = 1; i < N ; i++) {
		if (i % N2 == 0 && i % N1 == 1)
			e1 = i;	
		if (i % N1 == 0 && i % N2 == 1)
			e2 = i;
	}
	
}

void inline primes(int M, vector<int> &arr) {
	int num = 2;
	int mult;

	while (M > 1) {
		mult = num*6;
		int m1 = mult-1;
		int m2 = mult+1;
		while (M%m1 == 0 ) {
			arr.push_back(m1);
			M = M / m1;
		}
		while (M%m2 == 0 ) {
			arr.push_back(m2);
			M = M / m2;
		}
		num+=1;

	}


}

void inline factors(int M, vector<int> &arr) {
	int N = M;
	arr.resize(0);
	while (N%8 == 0){
			N = N/8;
			arr.push_back(8);
		}
	while (N%7 == 0){
			N = N/7;
			arr.push_back(7);
		}
	while (N%4 == 0){
		N = N/4;
		arr.push_back(4);
	}
	while (N%3 == 0){
			N = N/3;
			arr.push_back(3);
		}
	while (N%5 == 0){
			N = N/5;
			arr.push_back(5);
		}
	while (N%2 == 0){
			N = N/2;
			arr.push_back(2);
		}
	if (N > 8)
		primes(N,arr);

}


int inline factorf(int M) {
	int N = M;
	while (N%8 == 0){
			N = N/8;
		}
	while (N%7 == 0){
			N = N/7;
		}
	while (N%4 == 0){
		N = N/4;
	}
	while (N%3 == 0){
			N = N/3;
		}
	while (N%5 == 0){
			N = N/5;
		}
	while (N%2 == 0){
			N = N/2;
		}

	return N;
}


int inline findnext(int M) {
	int N = M;

	while (factorf(N) != 1) {
		++N;
	}

	return N;

}

void inline factorize(int N,int M, vector<int> &arr) {
	while (N % M == 0) {
		arr.push_back(M);
		N = N / M;
	}

}

template <typename T>
void inline cmapT(fft_data<T> &data,fft_data<T> &y,int n1, int n2) {
	int e1,e2;
	int N = n1*n2;
	idempotent(n1,n2,e1,e2);
	y = data;
	vector<int> index(n1,0);
	int temp;
	int temp2 = 0;
	for (int i = 0; i < n1; ++i) {
		while (temp2 >= N)
			temp2-= N;
		index[i] = temp2;
		temp2+=e1;
		
		//index[i] = (e1*i) % N;
	}
	for (int k = 0; k < n2; ++k) {
		int i = 0;
		for (int j = k*n1; j < (k+1)*n1; ++j) {
			y.re[j] = data.re[index[i]];
			y.im[j] = data.im[index[i]];
			++i;
		}
		temp = index[n1 - 1];
		
		for (int j = n1-1; j > 0; --j) {
			index[j] = index[j-1] + 1;
		}
		index[0] = temp+1;
	}	
	
}

template <typename T>
void inline cmap(fft_data<T> &data,fft_data<T> &y,int n1, int n2) {
	int e1,e2;
	int N = n1*n2;
	idempotent(n1,n2,e1,e2);
	y = data;
	vector<int> index(n1,0);
	//cout << data.re.size() << " MAX" << y.re.size() << endl;
	int temp;
	int temp2 = 0;
	for (int i = 0; i < n1; ++i) {
		while (temp2 >= N)
			temp2-= N;
		index[i] = temp2;
		temp2+=e1;
		
		//index[i] = (e1*i) % N;
	}
	for (int k = 0; k < n2; ++k) {
		int i = 0;
		for (int j = k*n1; j < (k+1)*n1; ++j) {
			//cout << j << " : " << index[i] << endl;
			y.re[index[i]] = data.re[j];
			y.im[index[i]] = data.im[j];
			++i;
		}
		temp = index[n1 - 1];
		
		for (int j = n1-1; j > 0; --j) {
			index[j] = index[j-1] + 1;
		}
		index[0] = temp+1;
	}	
	
}

template <typename T>
void inline rmapT(fft_data<T> &data,fft_data<T> &y,int n1, int n2) {
	int e1,e2;
	int N = n1*n2;
	idempotent(n1,n2,e1,e2);
	y = data;
	vector<int> index(n1,0);
	int temp = 0;
	
	for (int i = 0; i < n1; ++i) {
		index[i] = temp;
		temp+=n2;
		
	}
	
	for (int k = 0; k < n2; ++k) {
		int i = 0;
		for (int j = k*n1; j < (k+1)*n1; ++j) {
			y.re[j] = data.re[index[i]];
			y.im[j] = data.im[index[i]];
			++i;
		}
		
		//int temp2 = 0;
		for (int l = 0; l < n1; ++l) {
			index[l] = (index[l] + n1) % N;
		}
		
	}		
	
}	

template <typename T>
void inline rmap(fft_data<T> &data,fft_data<T> &y,int n1, int n2) {
	int e1,e2;
	int N = n1*n2;
	idempotent(n1,n2,e1,e2);
	y = data;
	vector<int> index(n1,0);
	int temp = 0;
	
	for (int i = 0; i < n1; ++i) {
		index[i] = temp;
		temp+=n2;
		
	}
	
	for (int k = 0; k < n2; ++k) {
		int i = 0;
		for (int j = k*n1; j < (k+1)*n1; ++j) {
			y.re[index[i]] = data.re[j];
			y.im[index[i]] = data.im[j];
			++i;
		}
		
		
		for (int l = 0; l < n1; ++l) {
			index[l] = (index[l] + n1) % N;
		}
		
	}		
	
}	

#endif
