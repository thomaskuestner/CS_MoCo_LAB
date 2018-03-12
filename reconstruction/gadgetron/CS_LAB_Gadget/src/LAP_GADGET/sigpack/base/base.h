// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_BASE_H
#define SP_BASE_H
namespace sp
{
    ///
    /// @defgroup math Math
    /// \brief Math functions.
    /// @{

    const double PI   = 3.14159265358979323846;     ///< ... _or use arma::datum::pi_
    const double PI_2 = 6.28318530717958647692;

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A sinc, sin(x)/x, function.
    /// @param x The angle in radians
    ////////////////////////////////////////////////////////////////////////////////////////////
    double sinc( double x )
    {
        if(x==0.0)
            return 1.0;
        else
            return std::sin(PI*x)/(PI*x);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A sinc, sin(x)/x, function.
    /// @param x The angle in radians
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec sinc(const arma::vec& x)
    {
        arma::vec out;
        out.copy_size(x);
        for (unsigned int n = 0; n < out.size(); n++)
        {
            out(n) = sinc(x(n));
        }
        return out;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Modified first kind bessel function order zero.
    ///
    /// See bessel functions on [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function)
    /// @param x
    ////////////////////////////////////////////////////////////////////////////////////////////
    double besseli0( double x )
    {
        double y=1.0,s=1.0,x2=x*x;
        int n = 1;
        while (s > y*1.0e-9)
        {
            s *= x2/4.0/(n*n);
            y += s;
            n++;
        }
        return y;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Calculates angle in radians for complex input.
    /// @param x Complex input value
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    double angle( const std::complex<T>& x )
    {
        return std::arg(x);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Calculates angle in radians for complex input.
    /// @param x Complex input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec angle( const arma::cx_vec& x )
    {
        arma::vec P;
        P.copy_size(x);
        for(unsigned int r=0; r<x.n_rows; r++)
            P(r) = std::arg(x(r));
        return P;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Calculates angle in radians for complex input.
    /// @param x Complex input matrix
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::mat angle( const arma::cx_mat& x )
    {
        arma::mat P;
        P.copy_size(x);
        for(unsigned int r=0; r<x.n_rows; r++)
            for(unsigned int c=0; c<x.n_cols; c++)
                P(r,c) = std::arg(x(r,c));
        return P;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Unwraps the angle vector x, accumulates phase.
    /// @param x Complex input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec unwrap( const arma::vec& x )
    {
        arma::vec P;
        double pacc = 0, pdiff = 0;
        const double thr=PI*170/180;
        P.copy_size(x);
        P(0)=x(0);
        for(unsigned int r=1; r<x.n_rows; r++)
        {
            pdiff = x(r)-x(r-1);
            if( pdiff >= thr ) pacc += -PI_2;
            if( pdiff <= -thr) pacc +=  PI_2;
            P(r)=pacc+x(r);
        }
        return P;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief 1D FFT shift.
    /// @returns Circular shifted FFT
    /// @param Pxx Complex FFT
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    arma::Col<T> fftshift(const arma::Col<T>& Pxx)
    {
        arma::Col<T> x(Pxx.n_elem);
        x = shift(Pxx, floor(Pxx.n_elem / 2));
        return x;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief 1D FFT inverse/reverse shift.
    /// @returns Circular shifted FFT
    /// @param Pxx Complex FFT
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    arma::Col<T> ifftshift(const arma::Col<T>& Pxx)
    {
        arma::Col<T> x(Pxx.n_elem);
        x = shift(Pxx, -ceil(Pxx.n_elem / 2));
        return x;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief 2D FFT shift.
    /// @returns Circular shifted FFT
    /// @param Pxx FFT
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    arma::Mat<T> fftshift(const arma::Mat<T>& Pxx)
    {
        arma::uword R = Pxx.n_rows;
        arma::uword C = Pxx.n_cols;
        arma::Mat<T> x(R, C);
        x = shift(Pxx, floor(R / 2), 0);
        x = shift(x, floor(C / 2), 1);
        return x;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief 2D FFT inverse/reverse shift.
    /// @returns Circular shifted FFT
    /// @param Pxx FFT
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    arma::Mat<T> ifftshift(const arma::Mat<T>& Pxx)
    {
        arma::uword R = Pxx.n_rows;
        arma::uword C = Pxx.n_cols;
        arma::Mat<T> x(R, C);
        x = shift(Pxx, -ceil(R / 2), 0);
        x = shift(x, -ceil(C / 2), 1);
        return x;
    }
    /// @}


    ///
    /// @defgroup misc Misc
    /// \brief Misc functions, such as error handling etc.
    /// @{

    ///////////////////////////////////
    // err_handler("Error string")
    //      Prints an error message, waits for input and
    //      then exits with error
#define err_handler(msg) \
    { \
        std::cout << "SigPack Error [" << __FILE__  << "@" << __LINE__ << "]: " << msg << std::endl; \
        std::cin.get(); \
        exit(EXIT_FAILURE);\
    }

    ///////////////////////////////////
    // wrn_handler("Warning string")
    //      Prints an warning message
#define wrn_handler(msg)  \
    { \
        std::cout << "SigPack warning [" << __FILE__ << "@" << __LINE__ << "]: " << msg << std::endl;\
    }
    /// @}

} // end namespace
#endif

