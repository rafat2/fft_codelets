/*
 * codelets3.h
 *
 *  Created on: Feb 21, 2013
 *      Author: USER
 */

#ifndef CODELETS3_H_
#define CODELETS3_H_

#include "codelets2.h"

template <typename T>
void inline radixN_mixed_dit_inplace(fft_data<T> &x, int q, int L, int N, int sgn) {
	int n = x.re.size();
	int Ls = L / N;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	int  M = (N-1) / 2;

	vector<T> wlr(N-1, (T) 1.0);
	vector<T> wli(N-1, (T) 0.0);
	vector<T> c1(N-1);
	vector<T> s1(N-1);
	vector<T> taur(N-1);
	vector<T> taui(N-1);

	for (int i = 1; i < M+1;++i) {
		c1[i-1] = cos(i*PI2/N);
		s1[i-1] = sin(i*PI2/N);
	}

	for (int i = 0; i < M;++i) {
		s1[i+M] = -s1[M-1-i];
		c1[i+M] =  c1[M-1-i];
	}



	T temp1r,temp1i,temp2r,temp2i;
	vector<T> yr(N);
	vector<T> yi(N);

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;

			yr[0] = x.re[index];
			yi[0] = x.im[index];

			for (int i = 1; i < N; ++i) {
				int index1 = index + Ls *i;
				yr[i] = wlr[i-1]*x.re[index1] - wli[i-1]*x.im[index1];
				yi[i] = wlr[i-1]*x.im[index1] + wli[i-1]*x.re[index1];
			}

			for (int i = 0; i < M; ++i) {
				taur[i] = yr[i+1] + yr[N-1-i];
				taui[i+M] = yi[i+1] - yi[N-1-i];
				taui[i] = yi[i+1] + yi[N-1-i];
				taur[i+M] = yr[i+1] - yr[N-1-i];
			}

			temp1r = yr[0];
			temp1i = yi[0];

			for (int i = 0; i < M; ++i) {
				temp1r+= taur[i];
				temp1i+= taui[i];
			}

			x.re[index] = temp1r;
			x.im[index] = temp1i;

			for (int u = 0; u < M; u++) {
				temp1r = yr[0];
				temp1i = yi[0];
				temp2r = (T) 0.0;
				temp2i = (T) 0.0;
				for (int v = 0; v < M; v++) {
					//int ind2 = (u+v)%M;
					int t = (u+1)*(v+1);
					while(t >= N)
						t-=N;
					int tt = t-1;

					temp1r+= c1[tt]*taur[v];
					temp1i+= c1[tt]*taui[v];
					temp2r-= s1[tt]*taur[v+M];
					temp2i-= s1[tt]*taui[v+M];
				}
				temp2r = sgn * temp2r;
				temp2i = sgn * temp2i;


				x.re[index + (u+1)*Ls] = temp1r - temp2i;
				x.im[index + (u+1)*Ls] = temp1i + temp2r;

				x.re[index + (N-u-1)*Ls] = temp1r + temp2i;
				x.im[index + (N-u-1)*Ls] = temp1i - temp2r;
			}

		}
		T temp = wlr[0];
		wlr[0] = C * wlr[0] - S * wli[0];
		wli[0] = S * temp + C * wli[0];

		wlr[1] = wlr[0]*wlr[0] - wli[0]*wli[0];
		wli[1] = 2.0*wlr[0]*wli[0];

		for (int i = 0; i < M-1; ++i) {
			int u = 2*(i+1);
			int v = u+1;
			int w = (u+1)/2;
			wlr[u] = wlr[w]*wlr[w-1] - wli[w-1]*wli[w];
			wli[u] = wlr[w]*wli[w-1] + wli[w]*wlr[w-1];

			wlr[v] = wlr[w]*wlr[w] - wli[w]*wli[w];
			wli[v] = 2.0*wlr[w]*wli[w];
		}

	}

}

template <typename T>
void inline radixN_mixed_dit_inplace(fft_data<T> &x, fft_data<T> &wl, int L, int N, int sgn) {
	int n = x.re.size();
	int Ls = L / N;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;

	int  M = (N-1) / 2;

	vector<T> wlr(N-1, (T) 1.0);
	vector<T> wli(N-1, (T) 0.0);
	vector<T> c1(N-1);
	vector<T> s1(N-1);
	vector<T> taur(N-1);
	vector<T> taui(N-1);

	for (int i = 1; i < M+1;++i) {
		c1[i-1] = cos(i*PI2/N);
		s1[i-1] = sin(i*PI2/N);
	}

	for (int i = 0; i < M;++i) {
		s1[i+M] = -s1[M-1-i];
		c1[i+M] =  c1[M-1-i];
	}



	T temp1r,temp1i,temp2r,temp2i;
	vector<T> yr(N);
	vector<T> yi(N);

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr[0] = wl.re[ind];
		wli[0] = wl.im[ind];
		wlr[1] = wlr[0]*wlr[0] - wli[0]*wli[0];
		wli[1] = 2.0*wlr[0]*wli[0];

		for (int i = 0; i < M-1; ++i) {
			int u = 2*(i+1);
			int v = u+1;
			int w = (u+1)/2;
			wlr[u] = wlr[w]*wlr[w-1] - wli[w-1]*wli[w];
			wli[u] = wlr[w]*wli[w-1] + wli[w]*wlr[w-1];

			wlr[v] = wlr[w]*wlr[w] - wli[w]*wli[w];
			wli[v] = 2.0*wlr[w]*wli[w];
		}
		for (int k =0; k < r; k++) {
			int index = k*L+j;

			yr[0] = x.re[index];
			yi[0] = x.im[index];

			for (int i = 1; i < N; ++i) {
				int index1 = index + Ls *i;
				yr[i] = wlr[i-1]*x.re[index1] - wli[i-1]*x.im[index1];
				yi[i] = wlr[i-1]*x.im[index1] + wli[i-1]*x.re[index1];
			}

			for (int i = 0; i < M; ++i) {
				taur[i] = yr[i+1] + yr[N-1-i];
				taui[i+M] = yi[i+1] - yi[N-1-i];
				taui[i] = yi[i+1] + yi[N-1-i];
				taur[i+M] = yr[i+1] - yr[N-1-i];
			}

			temp1r = yr[0];
			temp1i = yi[0];

			for (int i = 0; i < M; ++i) {
				temp1r+= taur[i];
				temp1i+= taui[i];
			}

			x.re[index] = temp1r;
			x.im[index] = temp1i;

			for (int u = 0; u < M; u++) {
				temp1r = yr[0];
				temp1i = yi[0];
				temp2r = (T) 0.0;
				temp2i = (T) 0.0;
				for (int v = 0; v < M; v++) {
					//int ind2 = (u+v)%M;
					int t = (u+1)*(v+1);
					while(t >= N)
						t-=N;
					int tt = t-1;

					temp1r+= c1[tt]*taur[v];
					temp1i+= c1[tt]*taui[v];
					temp2r-= s1[tt]*taur[v+M];
					temp2i-= s1[tt]*taui[v+M];
				}
				temp2r = sgn * temp2r;
				temp2i = sgn * temp2i;


				x.re[index + (u+1)*Ls] = temp1r - temp2i;
				x.im[index + (u+1)*Ls] = temp1i + temp2r;

				x.re[index + (N-u-1)*Ls] = temp1r + temp2i;
				x.im[index + (N-u-1)*Ls] = temp1i - temp2r;
			}

		}

	}

}

template <typename T>
void inline radixN_mixed_dif_inplace(fft_data<T> &x, int q, int L, int N, int sgn) {
	int n = x.re.size();
	int Ls = L / N;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	int  M = (N-1) / 2;

	vector<T> wlr(N-1, (T) 1.0);
	vector<T> wli(N-1, (T) 0.0);
	vector<T> c1(N-1);
	vector<T> s1(N-1);
	vector<T> taur(N-1);
	vector<T> taui(N-1);

	for (int i = 1; i < M+1;++i) {
		c1[i-1] = cos(i*PI2/N);
		s1[i-1] = sin(i*PI2/N);
	}

	for (int i = 0; i < M;++i) {
		s1[i+M] = -s1[M-1-i];
		c1[i+M] =  c1[M-1-i];
	}



	T temp1r,temp1i,temp2r,temp2i;
	vector<T> yr(N);
	vector<T> yi(N);

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;

			for (int i = 0; i < M; ++i) {
				taur[i] = x.re[index+(i+1)*Ls] + x.re[index+(N-1-i)*Ls];
				taui[i+M] = x.im[index+(i+1)*Ls] - x.im[index+(N-1-i)*Ls];
				taui[i] = x.im[index+(i+1)*Ls] + x.im[index+(N-1-i)*Ls];
				taur[i+M] = x.re[index+(i+1)*Ls] - x.re[index+(N-1-i)*Ls];
			}

			yr[0] = x.re[index];
			yi[0] = x.im[index];

			for (int u = 0; u < M; u++) {
				temp1r = yr[0];
				temp1i = yi[0];
				temp2r = (T) 0.0;
				temp2i = (T) 0.0;
				for (int v = 0; v < M; v++) {
					//int ind2 = (u+v)%M;
					int t = (u+1)*(v+1);
					while(t >= N)
						t-=N;
					int tt = t-1;

					temp1r+= c1[tt]*taur[v];
					temp1i+= c1[tt]*taui[v];
					temp2r-= s1[tt]*taur[v+M];
					temp2i-= s1[tt]*taui[v+M];
				}
				temp2r = sgn * temp2r;
				temp2i = sgn * temp2i;


				yr[u+1] = temp1r - temp2i;
				yi[u+1] = temp1i + temp2r;

				yr[N-u-1] = temp1r + temp2i;
				yi[N-u-1] = temp1i - temp2r;
			}

			for (int i = 1; i < N; ++i) {
				int index1 = index + Ls *i;
				x.re[index1] = wlr[i-1]*yr[i] - wli[i-1]*yi[i];
				x.im[index1] = wlr[i-1]*yi[i] + wli[i-1]*yr[i];
			}



			temp1r = yr[0];
			temp1i = yi[0];

			for (int i = 0; i < M; ++i) {
				temp1r+= taur[i];
				temp1i+= taui[i];
			}

			x.re[index] = temp1r;
			x.im[index] = temp1i;



		}
		T temp = wlr[0];
		wlr[0] = C * wlr[0] - S * wli[0];
		wli[0] = S * temp + C * wli[0];

		wlr[1] = wlr[0]*wlr[0] - wli[0]*wli[0];
		wli[1] = 2.0*wlr[0]*wli[0];

		for (int i = 0; i < M-1; ++i) {
			int u = 2*(i+1);
			int v = u+1;
			int w = (u+1)/2;
			wlr[u] = wlr[w]*wlr[w-1] - wli[w-1]*wli[w];
			wli[u] = wlr[w]*wli[w-1] + wli[w]*wlr[w-1];

			wlr[v] = wlr[w]*wlr[w] - wli[w]*wli[w];
			wli[v] = 2.0*wlr[w]*wli[w];
		}

	}

}

template <typename T>
void inline radixN_mixed_dif_inplace(fft_data<T> &x, fft_data<T> &wl, int L, int N, int sgn) {
	int n = x.re.size();
	int Ls = L / N;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;

	int  M = (N-1) / 2;

	vector<T> wlr(N-1, (T) 1.0);
	vector<T> wli(N-1, (T) 0.0);
	vector<T> c1(N-1);
	vector<T> s1(N-1);
	vector<T> taur(N-1);
	vector<T> taui(N-1);

	for (int i = 1; i < M+1;++i) {
		c1[i-1] = cos(i*PI2/N);
		s1[i-1] = sin(i*PI2/N);
	}

	for (int i = 0; i < M;++i) {
		s1[i+M] = -s1[M-1-i];
		c1[i+M] =  c1[M-1-i];
	}



	T temp1r,temp1i,temp2r,temp2i;
	vector<T> yr(N);
	vector<T> yi(N);

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr[0] = wl.re[ind];
		wli[0] = wl.im[ind];

		wlr[1] = wlr[0]*wlr[0] - wli[0]*wli[0];
		wli[1] = 2.0*wlr[0]*wli[0];

		for (int i = 0; i < M-1; ++i) {
			int u = 2*(i+1);
			int v = u+1;
			int w = (u+1)/2;
			wlr[u] = wlr[w]*wlr[w-1] - wli[w-1]*wli[w];
			wli[u] = wlr[w]*wli[w-1] + wli[w]*wlr[w-1];

			wlr[v] = wlr[w]*wlr[w] - wli[w]*wli[w];
			wli[v] = 2.0*wlr[w]*wli[w];
		}


		for (int k =0; k < r; k++) {
			int index = k*L+j;

			for (int i = 0; i < M; ++i) {
				taur[i] = x.re[index+(i+1)*Ls] + x.re[index+(N-1-i)*Ls];
				taui[i+M] = x.im[index+(i+1)*Ls] - x.im[index+(N-1-i)*Ls];
				taui[i] = x.im[index+(i+1)*Ls] + x.im[index+(N-1-i)*Ls];
				taur[i+M] = x.re[index+(i+1)*Ls] - x.re[index+(N-1-i)*Ls];
			}

			yr[0] = x.re[index];
			yi[0] = x.im[index];

			for (int u = 0; u < M; u++) {
				temp1r = yr[0];
				temp1i = yi[0];
				temp2r = (T) 0.0;
				temp2i = (T) 0.0;
				for (int v = 0; v < M; v++) {
					//int ind2 = (u+v)%M;
					int t = (u+1)*(v+1);
					while(t >= N)
						t-=N;
					int tt = t-1;

					temp1r+= c1[tt]*taur[v];
					temp1i+= c1[tt]*taui[v];
					temp2r-= s1[tt]*taur[v+M];
					temp2i-= s1[tt]*taui[v+M];
				}
				temp2r = sgn * temp2r;
				temp2i = sgn * temp2i;


				yr[u+1] = temp1r - temp2i;
				yi[u+1] = temp1i + temp2r;

				yr[N-u-1] = temp1r + temp2i;
				yi[N-u-1] = temp1i - temp2r;
			}

			for (int i = 1; i < N; ++i) {
				int index1 = index + Ls *i;
				x.re[index1] = wlr[i-1]*yr[i] - wli[i-1]*yi[i];
				x.im[index1] = wlr[i-1]*yi[i] + wli[i-1]*yr[i];
			}



			temp1r = yr[0];
			temp1i = yi[0];

			for (int i = 0; i < M; ++i) {
				temp1r+= taur[i];
				temp1i+= taui[i];
			}

			x.re[index] = temp1r;
			x.im[index] = temp1i;

		}


	}

}

template <typename T>
void inline sh_mixed_radixN_dit(fft_data<T> &x, int q, int L, int N, int sgn) {
	int n = x.re.size();
	int Ls = L / N;
	int r = n / L;
	int rs = n /Ls;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	int  M = (N-1) / 2;

	vector<T> wlr(N-1, (T) 1.0);
	vector<T> wli(N-1, (T) 0.0);
	vector<T> c1(N-1);
	vector<T> s1(N-1);
	vector<T> taur(N-1);
	vector<T> taui(N-1);

	fft_data<T> y = x;
	int lsr = Ls *r;

	for (int i = 1; i < M+1;++i) {
		c1[i-1] = cos(i*PI2/N);
		s1[i-1] = sin(i*PI2/N);
	}

	for (int i = 0; i < M;++i) {
		s1[i+M] = -s1[M-1-i];
		c1[i+M] =  c1[M-1-i];
	}



	T temp1r,temp1i,temp2r,temp2i;
	vector<T> yr(N);
	vector<T> yi(N);

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = j*rs+k;

			yr[0] = y.re[index];
			yi[0] = y.im[index];

			for (int i = 1; i < N; ++i) {
				int index1 = index + r *i;
				yr[i] = wlr[i-1]*y.re[index1] - wli[i-1]*y.im[index1];
				yi[i] = wlr[i-1]*y.im[index1] + wli[i-1]*y.re[index1];
			}

			for (int i = 0; i < M; ++i) {
				taur[i] = yr[i+1] + yr[N-1-i];
				taui[i+M] = yi[i+1] - yi[N-1-i];
				taui[i] = yi[i+1] + yi[N-1-i];
				taur[i+M] = yr[i+1] - yr[N-1-i];
			}

			temp1r = yr[0];
			temp1i = yi[0];

			for (int i = 0; i < M; ++i) {
				temp1r+= taur[i];
				temp1i+= taui[i];
			}

			int indexo = j*r+k;

			x.re[indexo] = temp1r;
			x.im[indexo] = temp1i;


			for (int u = 0; u < M; u++) {
				temp1r = yr[0];
				temp1i = yi[0];
				temp2r = (T) 0.0;
				temp2i = (T) 0.0;
				for (int v = 0; v < M; v++) {
					//int ind2 = (u+v)%M;
					int t = (u+1)*(v+1);
					while(t >= N)
						t-=N;
					int tt = t-1;

					temp1r+= c1[tt]*taur[v];
					temp1i+= c1[tt]*taui[v];
					temp2r-= s1[tt]*taur[v+M];
					temp2i-= s1[tt]*taui[v+M];
				}
				temp2r = sgn * temp2r;
				temp2i = sgn * temp2i;


				x.re[indexo + (u+1)*lsr] = temp1r - temp2i;
				x.im[indexo + (u+1)*lsr] = temp1i + temp2r;

				x.re[indexo + (N-u-1)*lsr] = temp1r + temp2i;
				x.im[indexo + (N-u-1)*lsr] = temp1i - temp2r;
			}

		}
		T temp = wlr[0];
		wlr[0] = C * wlr[0] - S * wli[0];
		wli[0] = S * temp + C * wli[0];

		wlr[1] = wlr[0]*wlr[0] - wli[0]*wli[0];
		wli[1] = 2.0*wlr[0]*wli[0];

		for (int i = 0; i < M-1; ++i) {
			int u = 2*(i+1);
			int v = u+1;
			int w = (u+1)/2;
			wlr[u] = wlr[w]*wlr[w-1] - wli[w-1]*wli[w];
			wli[u] = wlr[w]*wli[w-1] + wli[w]*wlr[w-1];

			wlr[v] = wlr[w]*wlr[w] - wli[w]*wli[w];
			wli[v] = 2.0*wlr[w]*wli[w];
		}

	}

}

template <typename T>
void inline sh_mixed_radixN_dit(fft_data<T> &x, fft_data<T> &wl, int L, int N, int sgn) {
	int n = x.re.size();
	int Ls = L / N;
	int r = n / L;
	int rs = n /Ls;

	T PI2 = (T) 6.28318530717958647692528676655900577;

	int  M = (N-1) / 2;

	vector<T> wlr(N-1, (T) 1.0);
	vector<T> wli(N-1, (T) 0.0);
	vector<T> c1(N-1);
	vector<T> s1(N-1);
	vector<T> taur(N-1);
	vector<T> taui(N-1);

	fft_data<T> y = x;
	int lsr = Ls *r;

	for (int i = 1; i < M+1;++i) {
		c1[i-1] = cos(i*PI2/N);
		s1[i-1] = sin(i*PI2/N);
	}

	for (int i = 0; i < M;++i) {
		s1[i+M] = -s1[M-1-i];
		c1[i+M] =  c1[M-1-i];
	}



	T temp1r,temp1i,temp2r,temp2i;
	vector<T> yr(N);
	vector<T> yi(N);

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr[0] = wl.re[ind];
		wli[0] = wl.im[ind];

		wlr[1] = wlr[0]*wlr[0] - wli[0]*wli[0];
		wli[1] = 2.0*wlr[0]*wli[0];

		for (int i = 0; i < M-1; ++i) {
			int u = 2*(i+1);
			int v = u+1;
			int w = (u+1)/2;
			wlr[u] = wlr[w]*wlr[w-1] - wli[w-1]*wli[w];
			wli[u] = wlr[w]*wli[w-1] + wli[w]*wlr[w-1];

			wlr[v] = wlr[w]*wlr[w] - wli[w]*wli[w];
			wli[v] = 2.0*wlr[w]*wli[w];
		}

		for (int k =0; k < r; k++) {
			int index = j*rs+k;

			yr[0] = y.re[index];
			yi[0] = y.im[index];

			for (int i = 1; i < N; ++i) {
				int index1 = index + r *i;
				yr[i] = wlr[i-1]*y.re[index1] - wli[i-1]*y.im[index1];
				yi[i] = wlr[i-1]*y.im[index1] + wli[i-1]*y.re[index1];
			}

			for (int i = 0; i < M; ++i) {
				taur[i] = yr[i+1] + yr[N-1-i];
				taui[i+M] = yi[i+1] - yi[N-1-i];
				taui[i] = yi[i+1] + yi[N-1-i];
				taur[i+M] = yr[i+1] - yr[N-1-i];
			}

			temp1r = yr[0];
			temp1i = yi[0];

			for (int i = 0; i < M; ++i) {
				temp1r+= taur[i];
				temp1i+= taui[i];
			}

			int indexo = j*r+k;

			x.re[indexo] = temp1r;
			x.im[indexo] = temp1i;


			for (int u = 0; u < M; u++) {
				temp1r = yr[0];
				temp1i = yi[0];
				temp2r = (T) 0.0;
				temp2i = (T) 0.0;
				for (int v = 0; v < M; v++) {
					//int ind2 = (u+v)%M;
					int t = (u+1)*(v+1);
					while(t >= N)
						t-=N;
					int tt = t-1;

					temp1r+= c1[tt]*taur[v];
					temp1i+= c1[tt]*taui[v];
					temp2r-= s1[tt]*taur[v+M];
					temp2i-= s1[tt]*taui[v+M];
				}
				temp2r = sgn * temp2r;
				temp2i = sgn * temp2i;


				x.re[indexo + (u+1)*lsr] = temp1r - temp2i;
				x.im[indexo + (u+1)*lsr] = temp1i + temp2r;

				x.re[indexo + (N-u-1)*lsr] = temp1r + temp2i;
				x.im[indexo + (N-u-1)*lsr] = temp1i - temp2r;
			}

		}

	}

}

template <typename T>
void inline sh_mixed_radixN_dif(fft_data<T> &x, int q, int L, int N, int sgn) {
	int n = x.re.size();
	int Ls = L / N;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;
	T theta = (T) -1.0 * sgn * PI2/L;

	T S = (T) sin(theta);
	T C = (T) cos(theta);

	int  M = (N-1) / 2;

	vector<T> wlr(N-1, (T) 1.0);
	vector<T> wli(N-1, (T) 0.0);
	vector<T> c1(N-1);
	vector<T> s1(N-1);
	vector<T> taur(N-1);
	vector<T> taui(N-1);

	for (int i = 1; i < M+1;++i) {
		c1[i-1] = cos(i*PI2/N);
		s1[i-1] = sin(i*PI2/N);
	}

	for (int i = 0; i < M;++i) {
		s1[i+M] = -s1[M-1-i];
		c1[i+M] =  c1[M-1-i];
	}



	T temp1r,temp1i,temp2r,temp2i;
	vector<T> yr(N);
	vector<T> yi(N);

	fft_data<T> y = x;
	int lsr = Ls*r;

	for (int j = 0; j < Ls; j++) {
		for (int k =0; k < r; k++) {
			int index = k*L+j;

			for (int i = 0; i < M; ++i) {
				taur[i] = y.re[index+(i+1)*Ls] + y.re[index+(N-1-i)*Ls];
				taui[i+M] = y.im[index+(i+1)*Ls] - y.im[index+(N-1-i)*Ls];
				taui[i] = y.im[index+(i+1)*Ls] + y.im[index+(N-1-i)*Ls];
				taur[i+M] = y.re[index+(i+1)*Ls] - y.re[index+(N-1-i)*Ls];
			}

			yr[0] = y.re[index];
			yi[0] = y.im[index];

			for (int u = 0; u < M; u++) {
				temp1r = yr[0];
				temp1i = yi[0];
				temp2r = (T) 0.0;
				temp2i = (T) 0.0;
				for (int v = 0; v < M; v++) {
					//int ind2 = (u+v)%M;
					int t = (u+1)*(v+1);
					while(t >= N)
						t-=N;
					int tt = t-1;

					temp1r+= c1[tt]*taur[v];
					temp1i+= c1[tt]*taui[v];
					temp2r-= s1[tt]*taur[v+M];
					temp2i-= s1[tt]*taui[v+M];
				}
				temp2r = sgn * temp2r;
				temp2i = sgn * temp2i;


				yr[u+1] = temp1r - temp2i;
				yi[u+1] = temp1i + temp2r;

				yr[N-u-1] = temp1r + temp2i;
				yi[N-u-1] = temp1i - temp2r;
			}

			int indexo = k*Ls+j;

			for (int i = 1; i < N; ++i) {
				int index1 = indexo + lsr *i;
				x.re[index1] = wlr[i-1]*yr[i] - wli[i-1]*yi[i];
				x.im[index1] = wlr[i-1]*yi[i] + wli[i-1]*yr[i];
			}



			temp1r = yr[0];
			temp1i = yi[0];

			for (int i = 0; i < M; ++i) {
				temp1r+= taur[i];
				temp1i+= taui[i];
			}

			x.re[indexo] = temp1r;
			x.im[indexo] = temp1i;



		}
		T temp = wlr[0];
		wlr[0] = C * wlr[0] - S * wli[0];
		wli[0] = S * temp + C * wli[0];

		wlr[1] = wlr[0]*wlr[0] - wli[0]*wli[0];
		wli[1] = 2.0*wlr[0]*wli[0];

		for (int i = 0; i < M-1; ++i) {
			int u = 2*(i+1);
			int v = u+1;
			int w = (u+1)/2;
			wlr[u] = wlr[w]*wlr[w-1] - wli[w-1]*wli[w];
			wli[u] = wlr[w]*wli[w-1] + wli[w]*wlr[w-1];

			wlr[v] = wlr[w]*wlr[w] - wli[w]*wli[w];
			wli[v] = 2.0*wlr[w]*wli[w];
		}

	}

}


template <typename T>
void inline sh_mixed_radixN_dif(fft_data<T> &x, fft_data<T> &wl, int L, int N, int sgn) {
	int n = x.re.size();
	int Ls = L / N;
	int r = n / L;

	T PI2 = (T) 6.28318530717958647692528676655900577;

	int  M = (N-1) / 2;

	vector<T> wlr(N-1, (T) 1.0);
	vector<T> wli(N-1, (T) 0.0);
	vector<T> c1(N-1);
	vector<T> s1(N-1);
	vector<T> taur(N-1);
	vector<T> taui(N-1);

	for (int i = 1; i < M+1;++i) {
		c1[i-1] = cos(i*PI2/N);
		s1[i-1] = sin(i*PI2/N);
	}

	for (int i = 0; i < M;++i) {
		s1[i+M] = -s1[M-1-i];
		c1[i+M] =  c1[M-1-i];
	}



	T temp1r,temp1i,temp2r,temp2i;
	vector<T> yr(N);
	vector<T> yi(N);

	fft_data<T> y = x;
	int lsr = Ls*r;

	for (int j = 0; j < Ls; j++) {
		int ind = j*r;
		wlr[0] = wl.re[ind];
		wli[0] = wl.im[ind];
		wlr[1] = wlr[0]*wlr[0] - wli[0]*wli[0];
		wli[1] = 2.0*wlr[0]*wli[0];

		for (int i = 0; i < M-1; ++i) {
			int u = 2*(i+1);
			int v = u+1;
			int w = (u+1)/2;
			wlr[u] = wlr[w]*wlr[w-1] - wli[w-1]*wli[w];
			wli[u] = wlr[w]*wli[w-1] + wli[w]*wlr[w-1];

			wlr[v] = wlr[w]*wlr[w] - wli[w]*wli[w];
			wli[v] = 2.0*wlr[w]*wli[w];
		}

		for (int k =0; k < r; k++) {
			int index = k*L+j;

			for (int i = 0; i < M; ++i) {
				taur[i] = y.re[index+(i+1)*Ls] + y.re[index+(N-1-i)*Ls];
				taui[i+M] = y.im[index+(i+1)*Ls] - y.im[index+(N-1-i)*Ls];
				taui[i] = y.im[index+(i+1)*Ls] + y.im[index+(N-1-i)*Ls];
				taur[i+M] = y.re[index+(i+1)*Ls] - y.re[index+(N-1-i)*Ls];
			}

			yr[0] = y.re[index];
			yi[0] = y.im[index];

			for (int u = 0; u < M; u++) {
				temp1r = yr[0];
				temp1i = yi[0];
				temp2r = (T) 0.0;
				temp2i = (T) 0.0;
				for (int v = 0; v < M; v++) {
					//int ind2 = (u+v)%M;
					int t = (u+1)*(v+1);
					while(t >= N)
						t-=N;
					int tt = t-1;

					temp1r+= c1[tt]*taur[v];
					temp1i+= c1[tt]*taui[v];
					temp2r-= s1[tt]*taur[v+M];
					temp2i-= s1[tt]*taui[v+M];
				}
				temp2r = sgn * temp2r;
				temp2i = sgn * temp2i;


				yr[u+1] = temp1r - temp2i;
				yi[u+1] = temp1i + temp2r;

				yr[N-u-1] = temp1r + temp2i;
				yi[N-u-1] = temp1i - temp2r;
			}

			int indexo = k*Ls+j;

			for (int i = 1; i < N; ++i) {
				int index1 = indexo + lsr *i;
				x.re[index1] = wlr[i-1]*yr[i] - wli[i-1]*yi[i];
				x.im[index1] = wlr[i-1]*yi[i] + wli[i-1]*yr[i];
			}



			temp1r = yr[0];
			temp1i = yi[0];

			for (int i = 0; i < M; ++i) {
				temp1r+= taur[i];
				temp1i+= taui[i];
			}

			x.re[indexo] = temp1r;
			x.im[indexo] = temp1i;

		}

	}

}

template <typename T>
void inline ss_mixed_radix2_dif(fft_data<T> &x,int q, int L, int sgn) {
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


#endif /* CODELETS3_H_ */
