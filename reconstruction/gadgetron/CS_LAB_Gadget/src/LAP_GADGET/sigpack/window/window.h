// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_WINDOW_H
#define SP_WINDOW_H
namespace sp
{
    ///
    /// @defgroup window Window
    /// \brief Window functions.
    ///
    /// See window functions at [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    /// @{

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Generic fifth order symmetric cos window.
    ///
    /// \f$ w_i = a_0-a_1\ cos(2\pi i /(N-1))+a_2\ cos(4\pi i /(N-1))-a_3\ cos(6\pi i /(N-1))+a_4\ cos(8\pi i /(N-1))\f$
    /// @returns The cosinus window based on the <b>a</b> vector
    /// @param N Number of window taps
    /// @param a A vector of cosinus coefficients
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec cos_win( const arma::uword N, const arma::vec& a )
    {
        arma::vec h(N);
        for(arma::uword i=0; i<N; i++)
        {
            h[i] = a[0] - a[1]*std::cos(1.0*PI_2*i/(N-1)) + a[2]*std::cos(2.0*PI_2*i/(N-1)) \
                   - a[3]*std::cos(3.0*PI_2*i/(N-1)) + a[4]*std::cos(4.0*PI_2*i/(N-1));
        }
        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Hamming window.
    ///
    /// \f$ w_i = 0.54-0.46\ cos(2\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec hamming( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.54;
        a[1] = 0.46;
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Hann window.
    ///
    /// \f$ w_i = 0.5-0.5\ cos(2\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec hann( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.5;
        a[1] = 0.5;
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Blackman window.
    ///
    /// \f$ w_i = 0.42-0.5\ cos(2\pi i /(N-1))+0.08\ cos(4\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec blackman( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.42; // 7938/18608.0
        a[1] = 0.5;  // 9240/18608.0
        a[2] = 0.08; // 1430/18608.0
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Blackman-Harris window.
    /// Symmetric BH4 window
    ///
    /// \f$ w_i = 0.359-0.488\ cos(2\pi i /(N-1))+0.141\ cos(4\pi i /(N-1))-0.011\ cos(6\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec blackmanharris( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.35875;
        a[1] = 0.48829;
        a[2] = 0.14128;
        a[3] = 0.01168;
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Flattop window.
    ///
    /// \f$ w_i = 0.216-0.417\ cos(2\pi i /(N-1))+0.277\ cos(4\pi i /(N-1))-0.084\ cos(6\pi i /(N-1))+0.007\ cos(8\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec flattopwin( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.21557895;
        a[1] = 0.41663158;
        a[2] = 0.277263158;
        a[3] = 0.083578947;
        a[4] = 0.006947368;
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Hanning window.
    ///
    /// \f$ w_i = 0.5-0.5\ cos(2\pi (i+1) /(N+1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec hanning( const arma::uword N )
    {
        arma::vec h(N);
        for(arma::uword i=0; i<N; i++)
        {
            h[i] = 0.5-0.5*std::cos(PI_2*(i+1)/(N+1));
        }
        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Kaiser window.
    ///
    /// See Kaiser window at [Wikipedia](https://en.wikipedia.org/wiki/Window_function#Kaiser_window)
    /// @param N Nr of taps
    /// @param beta Beta factor
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec kaiser( const arma::uword N, double beta )
    {
        arma::vec h(N);
        double bb = besseli0(beta);
        for( arma::uword i=0; i<N; i++)
        {
            h[i] = besseli0(beta*sqrt(4.0*i*(N-1-i))/(N-1))/bb;
        }
        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Triangle window.
    ///
    /// See Triangle window at [Wikipedia](https://en.wikipedia.org/wiki/Window_function#Triangular_window)
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec triang( const arma::uword N )
    {
        arma::vec h(N);
        if(N%2)    // Odd
        {
            for(arma::uword i=0; i<(N-1)/2; i++)
            {
                h[i]     = 2.0*(i+1)/(N+1);
                h[N-i-1] = h[i];
            }
            h[(N-1)/2] = 1.0;
        }
        else      // Even
        {
            for(arma::uword i=0; i<N/2; i++)
            {
                h[i]     = (2.0*i+1)/N;
                h[N-i-1] = h[i];
            }
        }
        return h;
    }
    /// @}
} // end namespace
#endif
