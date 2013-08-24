#include <stdio.h>
#include <complex>
#include "src/fft_var.h"
#include "src/tools/cycle.h"


int main(int argc, char **argv)
{
	//int radix = 13;
	int N = 2048;
	//vector<complex<double> > sig1;
	fft_data<double> sig2,sig1,sig3;
	
	for (int i =0; i < N; i++){
	//sig1.push_back(complex<double>((double)1.0, 0.0));
		sig2.re.push_back((double) i);
		sig2.im.push_back((double) i);
		sig1.re.push_back((double) i);
		sig1.im.push_back((double) i);
		//cout << real(sig1[i]) << " " << imag(sig1[i]) << endl;
	}

    //sig1 = sig2;
	double tdl, tct;
	
	ticks t0 = getticks();
    //int N1 = 4, N2=5;
	//fftsh_radixN_dif(sig1,1,N,radix);
	//longvector(sig3,N);


	fft_split_radix_rec(sig1,1,N);
	ticks t1 = getticks();
	ticks t2 = getticks();
	//fft_split_radix_rec(sig1,-1,N);
	//fftct_radix2_dit_rec(sig3,sig1,-1,N);
	//fftct_radix2_dit_rec(sig1,-1,N);
	fft_bluestein(sig2,1);
	//fftct_radixN_dif_inplace(sig1,-1,N,radix);
	//fftct_radixN_dit_inplace(sig1,-1,N,radix);
	//fftsh_radix8_dif(sig1,-1,N);
	//fftct_radix8_dif_inplace(sig1,-1,N);
	//fftsh_radix7_dif(sig1,-1,N);
	//twiddle_rcbs(sig1,N);
	ticks t3 = getticks();

	tdl = elapsed(t1,t0);
	
	tct = elapsed(t3,t2);
	

	for (int i =0; i < N; i++){
		cout << sig1.re[i] - sig2.re[i] << " " << sig1.im[i] - sig2.im[i]  << endl;
	}

        //IFFT
/*
        fft_bluestein(sig2,-1);
        //fft(sig1,-1,N);
        cout << "IFFT - signal" << endl;
        for (int i =0; i < N; i++){
                cout << (sig2.re[i]/ N) << " " << (sig2.im[i]/N) << endl;
        }
*/        
    cout << "tdl : " << tdl << "  " << "tct : " << tct << endl;
    
    //cout << sig3.re.size() << endl;
     
	return 0;

}
