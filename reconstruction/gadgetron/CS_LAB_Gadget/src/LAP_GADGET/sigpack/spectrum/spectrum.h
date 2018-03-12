// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_SPECTRUM_H
#define SP_SPECTRUM_H
namespace sp
{
    ///
    /// @defgroup spectrum Spectrum
    /// \brief Spectrum functions.
    /// @{

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Windowed spectrum calculation.
    ///
    /// The spectrum is calculated using the fast fourier transform of the windowed input data vector
    /// @returns A complex spectrum vector
    /// @param x Input vector
    /// @param W Window function vector. NB! Must be same size as input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::cx_vec spectrum(const arma::Col<T1>& x, const arma::vec& W)
    {
        arma::cx_vec Pxx(x.size());
        double wc = sum(W);     // Window correction factor
        Pxx = fft(x % W)/wc;    // FFT calc
        return Pxx;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Power spectrum density calculation using windowed data.
    /// @returns A real valued PSD vector
    /// @param x Input vector
    /// @param W Window function vector. NB! Must be same size as input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::vec psd(const arma::Col<T1>& x, const arma::vec& W)
    {
        arma::cx_vec X(x.size());
        arma::vec Pxx(x.size());
        X = spectrum(x,W);          // FFT calc
        Pxx = real(X % conj(X));    // Calc power spectra
        return Pxx;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Power spectrum density calculation using Hamming windowed data.
    /// @returns A real valued PSD vector
    /// @param x Input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::vec psd(const arma::Col<T1>& x)
    {
        arma::vec W;
        W = hanning(x.size());
        return psd(x,W);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Spectrogram calculation using Hamming windowed data.
    ///
    /// See spectrogram at [Wikipedia](https://en.wikipedia.org/wiki/Spectrogram)
    /// @returns A complex spectrogram matrix
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::cx_mat specgram_cx(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::cx_mat Pw;

        //Def params
        arma::uword N = x.size();
        arma::uword D = Nfft-Noverl;
        arma::uword m = 0;
        if(N > Nfft)
        {
            arma::Col<T1> xk(Nfft);
            arma::vec W(Nfft);

            W = hanning(Nfft);
            arma::uword U = floor((N-Noverl)/double(D));
            Pw.set_size(Nfft,U);
            Pw.zeros();

            // Avg loop
            for(arma::uword k=0; k<N-Nfft; k+=D)
            {
                xk = x.rows(k,k+Nfft-1);       // Pick out chunk
                Pw.col(m++) = spectrum(xk,W);  // Calculate spectrum
            }
        }
        else
        {
            arma::vec W(N);
            W = hanning(N);
            Pw.set_size(N,1);
            Pw = spectrum(x,W);
        }
        return Pw;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Power spectrogram calculation.
    ///
    /// See spectrogram at [Wikipedia](https://en.wikipedia.org/wiki/Spectrogram)
    /// @returns A power spectrogram matrix
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::mat specgram(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::cx_mat Pw;
        arma::mat Sg;
        Pw = specgram_cx(x,Nfft,Noverl);
        Sg = real(Pw % conj(Pw));              // Calculate power spectrum
        return Sg;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Phase spectrogram calculation.
    ///
    /// See spectrogram at [Wikipedia](https://en.wikipedia.org/wiki/Spectrogram)
    /// @returns A phase spectrogram matrix
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::mat specgram_ph(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::cx_mat Pw;
        arma::mat Sg;
        Pw = specgram_cx(x,Nfft,Noverl);
        Sg = angle(Pw);                        // Calculate phase spectrum
        return Sg;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Phase spectrum calculation using Welsh's method.
    ///
    /// See Welsh's method at [Wikipedia](https://en.wikipedia.org/wiki/Welch%27s_method)
    /// @returns A phase spectrum vector
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::vec pwelch_ph(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::mat Ph;
        Ph  = specgram_ph(x,Nfft,Noverl);
        return arma::mean(Ph,1);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Power spectrum calculation using Welsh's method.
    ///
    /// _abs(pwelch(x,Nfft,Noverl))_ is equivalent to Matlab's: _pwelch(x,Nfft,Noverl,'twosided','power')_ <br>
    /// See Welsh's method at [Wikipedia](https://en.wikipedia.org/wiki/Welch%27s_method)
    /// @returns A power spectrum vector
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::vec pwelch(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::mat Pxx;
        Pxx = specgram(x,Nfft,Noverl);
        return arma::mean(Pxx,1);
    }
    /// @}
} // end namepace
#endif
