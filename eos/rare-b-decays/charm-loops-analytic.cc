/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Viktor Kuschke
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59.0 Temple
 * Place, Suite 330, Boston, MA  02111-1307.0  USA
 */

#include <eos/maths/power-of.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/multiplepolylog-li22.hh>

#include <eos/rare-b-decays/charm-loops.hh>

#include <eos/utils/exception.hh>
#include <eos/utils/log.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <complex>

namespace eos
{
    using std::atan;
    using std::complex;
    using std::log;
    using std::pow;
    using std::sqrt;
    using std::abs;
    using std::conj;
    using std::min;

    namespace agv_2019a
    {
        // Often used constants

        complex<double> myi(0.0, 1.0);
        const double ln2        = std::log(2.0);
        const double ln2squ     = power_of<2>(std::log(2.0));
        const double ln2cube    = power_of<3>(std::log(2.0));
        const double pisqu      = power_of<2>(M_PI);
        const double zeta3      = 1.2020569031595943;

        // Loop variables from [AGV:2019A]

        inline complex<double> s_eps(const complex<double> & s, const double & feynepsilonhat)
        {
            return s / (1.0 - feynepsilonhat * myi);
        }

        inline complex<double> z_eps(const double & z, const double & feynepsilonhat)
        {
            return (z - feynepsilonhat * myi) / (1.0 - feynepsilonhat * myi);
        }

        inline complex<double> x_a(const double & z, const double & feynepsilonhat)
        {
            if (z == 0.25)
            {
                throw InternalError("x_a becomes singular for m_c = m_b / 2");
            }
            else
            {
                return 1.0 / std::sqrt(1.0 - 4.0 * z_eps(z, feynepsilonhat));
            }
        }

        inline complex<double> x_b(const double & z, const double & feynepsilonhat)
        {
            return std::sqrt(4.0 * z_eps(z, feynepsilonhat)) - std::sqrt(4.0 * z_eps(z, feynepsilonhat) - 1.0);
        }

        /* Not needed
        complex<double> x_c(const double & z, const double & feynepsilonhat)
        {
            return x_a(z, feynepsilonhat);
        }

        complex<double> x_d(const double & z, const double & feynepsilonhat)
        {
            return x_b(z, feynepsilonhat);
        }

        complex<double> x_e(const double & z, const double & feynepsilonhat)
        {
            return x_a(z, feynepsilonhat);
        }
        */

        inline complex<double> y_a(const complex<double> & s, const double & z, const double & feynepsilonhat)
        {
            if (1.0 - s_eps(s, feynepsilonhat) == 0.0)
            {
                return 0.0;
            }
            else
            {
                return 1.0 / std::sqrt(1.0 - 4.0 * z_eps(z, feynepsilonhat) / (1.0 - s_eps(s, feynepsilonhat)));
            }
        }

        complex<double> y_b(const complex<double> & s, const double & feynepsilonhat)
        {
            if (s_eps(s, feynepsilonhat) == 0.0)
            {
                return 0.0;
            }
            else
            {
                return 1.0 / std::sqrt(1.0 - 4.0 / s_eps(s, feynepsilonhat));
            }
        }

        complex<double> y_c(const complex<double> & s, const double & z, const double & feynepsilonhat)
        {
            if (s_eps(s, feynepsilonhat) == 0.0)
            {
                return 0.0;
            }
            else
            {
                return 1.0 / std::sqrt(1.0 - 4.0 * z_eps(z, feynepsilonhat) / s_eps(s, feynepsilonhat));
            }
        }

        /* Not needed
        inline complex<double> y_d(const complex<double> & s, const double & z, const double & feynepsilonhat)
        {
            return y_c(s, z, feynepsilonhat);
        }

        inline complex<double> y_e(const complex<double> & s, const double & z, const double & feynepsilonhat)
        {
            return y_c(s, z, feynepsilonhat);
        }
        */

        // Constructor of the loop parameters
        CharmLoopsParameters::CharmLoopsParameters(const double & muhat, const complex<double> & s, const double & z, const double & feynepsilonhat) :
            muhat(muhat)
        {
            s_eps = agv_2019a::s_eps(s, feynepsilonhat);
            z_eps = agv_2019a::z_eps(z, feynepsilonhat);

            x_a = agv_2019a::x_a(z, feynepsilonhat);
            x_b = agv_2019a::x_b(z, feynepsilonhat);
            x_c = x_a;
            x_d = x_b;
            x_e = x_a;

            y_a = agv_2019a::y_a(s, z, feynepsilonhat);
            y_b = agv_2019a::y_b(s, feynepsilonhat);
            y_c = agv_2019a::y_c(s, z, feynepsilonhat);
            y_d = y_c;
            y_e = y_c;
        }

        // Signum
        inline double my_sign(const double x)
        {
            if (x < 0.0)
            {
                return - 1.0;
            }
            if (x > 0.0)
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        }

        // Functions depending on Heaviside step function

        // Triangle function T(a,b;x) from [FTW:2016A]
        inline double T(const complex<double> a, const complex<double> b, const complex<double> x)
        {
            const complex<double> amb = a - b;
            const complex<double> xconj = std::conj(x);
            const complex<double> bconj = std::conj(b);
            const double denom = (xconj * amb).imag();

            if (denom == 0.0)
            {
                return 0.0;
            }
            else
            {
                const double arg1 = (xconj * a).imag() / denom;
                const double arg2 = 1.0 - (xconj * a).imag() / denom;
                const double arg3 = - 1.0 - (bconj * a).imag() / denom;

                if ( arg1 < 0.0 || arg2 < 0.0 || arg1 > 0.0)
                {
                    return 0.0;
                }
                if (arg1 == 0.0 || arg2 == 0.0 || arg3 == 0.0)
                {
                    throw InternalError("Ill-defined Theta(0.0)");
                }
                else
                {
                    return 1.0;
                }
            }
        }

        // point projected from x in side the triangle 0ab onto the edge ab
        inline complex<double> p(const complex<double> a, const complex<double> b, const complex<double> x)
        {
            const complex<double> amb = a - b;
            const complex<double> aconj = std::conj(a);
            const complex<double> xconj = std::conj(x);
            const complex<double> aconjb = aconj * b;
            const double denom = (xconj * amb).imag();

            if (denom == 0.0) // is this the best criterion or should T(a,b,x)!=1.0 be taken?
            {
                return 0.0; // should there be thrown an error, if the point x does not lie inside the triangel 0ab?
            }
            else
            {
                return x * aconjb.imag() / denom;
            }
        }

        // r function from [FTW:2016A]
        inline double r(const complex<double> & a, const complex<double> & b)
        {
            const complex<double> aconjb = std::conj(a) * b;
            const double denom = aconjb.imag();
            const double num = power_of<2>(std::abs(a)) * b.imag() - power_of<2>(std::abs(b)) * a.imag();

            if (denom == 0.0)
            {
                throw InternalError("0 in denominator");
            }
            else
            {
                return num / denom;
            }
        }

        // H1 and H2 function from [FTW:2016A]
        inline double H1(const complex<double> & a, const complex<double> & b)
        {
            const complex<double> aconjb = std::conj(a) * b;

            if (aconjb.imag() == 0.0)
            {
                return 0.0;
            }
            else
            {
                double min;

                if (1.0 <= power_of<2>(std::abs(a)) * b.imag() / aconjb.imag())
                {
                    min = 1.0;
                }
                else
                {
                   min = power_of<2>(std::abs(a)) * b.imag() / aconjb.imag();
                }

                const double r_val = agv_2019a::r(a,b);
                const double arg = min - r_val;

                if (r_val < 0.0 || arg < 0.0)
                {
                    return 0.0;
                }
                if (r_val == 0.0 || arg == 0.0)
                {
                    throw InternalError("Ill-defined Theta(0.0)");
                }
                else
                {
                    return 1.0;
                }
            }
        }

        inline double H2(const complex<double> & a, const complex<double> & b)
        {
            const complex<double> aconjb = std::conj(a) * b;

            if (aconjb.imag() == 0.0)
            {
                return 0.0;
            }
            else
            {
                const double r_val = agv_2019a::r(a,b);
                const double imaimb = a.imag() * b.imag();

                if (r_val < 0.0 || r_val > 1.0 || imaimb > 0.0)
                {
                    return 0.0;
                }
                if (r_val == 0.0 || r_val == 1.0 || imaimb == 0.0)
                {
                    throw InternalError("Ill-defined Theta(0.0)");
                }
                else
                {
                    return 1.0;
                }
            }
        }

        // Leading order one-loop functions

        complex<double> f290(const CharmLoopsParameters & clp)
        {
            const double q_c = 2.0 / 3.0; // charm quark charge

            const complex<double> result = 2.0 / 3.0 + 4.0 * clp.s_eps / clp.z_eps + std::log(power_of<2>(clp.muhat) / clp.z_eps)
                + (1.0 - 3.0 * power_of<2>(clp.y_c)) / (2.0 * power_of<3>(clp.y_c)) * (std::log(1.0 + clp.y_c) - std::log(1.0 - clp.y_c));

            return 2.0 / 3.0 * q_c * result;
        }

        complex<double> f190(const CharmLoopsParameters & clp)
        {
            const double c_F = 4.0 / 3.0; // SU(3.0) color factor

            return c_F * f290(clp);
        }

        // Counterterms to the two-loop functions

        complex<double> f17ctQs(const CharmLoopsParameters & clp)
        {
            return 0.0;
        }

        complex<double> f17ctQc(const CharmLoopsParameters & clp)
        {
            return 0.0;
        }

        complex<double> f17ctQb(const CharmLoopsParameters & clp)
        {
            const double muhat = clp.muhat;
            const complex<double> yb = clp.y_b;

            const complex<double> result = - 8.0 + 4.0 * (std::log(1.0 + yb) - std::log(1.0 - yb)) / yb - 8.0 * std::log(muhat);

            return result / 81.0;
        }

        complex<double> f27ctQs(const CharmLoopsParameters & clp)
        {
            return 0.0;
        }

        complex<double> f27ctQc(const CharmLoopsParameters & clp)
        {
            return 0.0;
        }

        complex<double> f27ctQb(const CharmLoopsParameters & clp)
        {
            return - 6.0 * f17ctQb(clp);
        }

        complex<double> f19ctQs(const CharmLoopsParameters & clp)
        {
            const complex<double> s = clp.s_eps;
            const complex<double> lnms = std::log(-s);
            const complex<double> lnmuhat = std::log(clp.muhat);

            const complex<double> result = - 104.0 / 9.0 + 2.0 * pisqu / 3.0 - 32.0 * lnmuhat / 3.0 - 16.0 * power_of<2>(lnmuhat) + 16.0 * lnms / 3.0
                + 16.0 * lnmuhat * lnms - 4.0 * power_of<2>(lnms);

            return result / 243.0;
        }

        complex<double> f19ctQc(const CharmLoopsParameters & clp)
        {
            const complex<double> muhat = clp.muhat;
            const complex<double> xa    = clp.x_a;
            const complex<double> yc    = clp.y_c;

            const complex<double> yc2 = power_of<2>(yc);
            const complex<double> yc3 = power_of<3>(yc);

            const complex<double> lnmuhat   = std::log(muhat);
            const complex<double> lnz       = std::log(clp.z_eps);
            const complex<double> lnxa      = std::log(xa);
            const complex<double> ln1pxa    = std::log(1.0 + xa);
            const complex<double> ln1mxa    = std::log(1.0 - xa);
            const complex<double> ln1pyc    = std::log(1.0 + yc);
            const complex<double> ln1myc    = std::log(1.0 - yc);

            const complex<double> result = 1792.0 / 3.0 + (256.0 * myi * M_PI) - (20.0 * power_of<2>(M_PI)) / 3.0 - (512.0 + 208.0 * myi * M_PI) / yc2 + (512.0 * ln2) + (32.0 * myi * M_PI * ln2)
                - (416.0 * ln2) / yc2 + (32.0 * power_of<2>(ln2)) + (896.0 * lnmuhat) + (32.0 * myi * M_PI * lnmuhat ) - (800.0 * lnmuhat) / yc2
                + (64.0 * ln2 * lnmuhat) + (32.0 * power_of<2>(lnmuhat)) - (192.0 * lnz) + (192.0 * lnz) / yc2 + (32.0 * dilog(1.0 / 2.0)) - (16.0 * dilog((1.0 - xa) / 2.0))
                - (16.0 * dilog((1.0 + xa) / 2.0)) + (-104.0 / yc3 + 216.0 / yc - (96.0 * yc)) * dilog((1.0 - yc) / 2.0) + (104.0 / yc3 - 216.0 / yc + (96.0 * yc)) * dilog((1.0 + yc) / 2.0)
                + (8.0 * power_of<2>(ln1mxa)) + ln1mxa * (-256.0 - 16.0 * myi * M_PI + 208.0 / yc2 - (16.0 * ln2) - (32.0 * lnmuhat) - (32.0 * lnxa))
                + (32.0 * power_of<2>(lnxa)) + lnxa * (512.0 + (32.0 * myi * M_PI) - 416.0 / yc2 + (64.0 * ln2) + (64.0 * lnmuhat) - (32.0 * ln1pxa))
                + (-256.0 - (16.0 * myi * M_PI) + 208.0 / yc2 - (16.0 * ln2) - (32.0 * lnmuhat)) * ln1pxa + (8.0 * power_of<2>(ln1pxa))
                + (-256.0 / yc3 - (104.0 * myi * M_PI) / yc3 + 336.0 / yc + (216.0 * myi * M_PI) / yc - (32.0 * yc) - (96.0 * myi * M_PI * yc) / 9.0 - (208.0 * ln2) / yc3 + (432.0 * ln2) / yc
                    - (192.0 * yc * ln2) - (400.0 * lnmuhat) / yc3 + (816.0 * lnmuhat) / yc - (384.0 * yc * lnmuhat) + (96.0 * lnz) / yc3 - (192.0 * lnz) / yc + (96.0 * yc * lnz)
                    + (104.0 * ln2) / yc3 - (216.0 * ln2) / yc + (96.0 * yc * ln2) + (104.0 / yc3 - 216.0 / yc + 32.0 * yc) * ln1mxa
                    + (-208.0 / yc3 + 16.0 / yc - 192.0 * yc) * lnxa + (104.0 / yc3 - 216.0 / yc + 96.0 * yc) * ln1pxa) * ln1myc
                + (-52.0 / yc3 + 108.0 / yc - (48.0 * yc)) * power_of<2>(ln1myc)
                + (256.0 / yc3 + (104.0 * myi * M_PI) / yc3 - 336.0 / yc - (216.0 * myi * M_PI) / yc + (32.0 + 96.0 * myi * M_PI) * yc
                    + (208.0 * ln2) / yc3 - (432.0 * ln2) / yc + (192.0 * yc * ln2) + (400.0 * lnmuhat) / yc3 - (816.0 * lnmuhat) / yc + (384.0 * yc * lnmuhat)
                    - (96.0 * lnz) / yc3 + (192.0 * lnz) / yc - (96.0 * yc * lnz) - (104.0 * ln2) / yc3 + (216.0 * ln2) / yc - (96.0 * yc * ln2)
                    + (-104.0 / yc3 + 8.0 / yc - 96.0 * yc) * ln1mxa + (208.0 / yc3 - 432.0 / yc + 192.0 * yc) * lnxa + (-104.0 / yc3 + 216.0 / yc - 96.0 * yc) * ln1pxa) * ln1pyc
                + (52.0 / yc3 - 108.0 / yc + (48.0 * yc)) * power_of<2>(ln1pyc);

            return result / 27.0;
        }

        complex<double> f19ctQb(const CharmLoopsParameters & clp)
        {
            const complex<double> muhat = clp.muhat;
            const complex<double> yb    = clp.y_b;

            const complex<double> yb3 = power_of<3>(yb);

            const complex<double> lnmub     = std::log(muhat);
            const complex<double> ln1pyb    = std::log(1.0 + yb);
            const complex<double> ln1myb    = std::log(1.0 - yb);
            const complex<double> dilog1pybhalf = dilog((1.0 + yb) / 2.0);
            const complex<double> dilog1mybhalf = dilog((1.0 - yb) / 2.0);

            const complex<double> result = -224.0 / 27.0 - (2.0 * power_of<2>(M_PI)) / 9.0 + 40.0 / (9.0 * power_of<2>(yb))
                - (80.0 * lnmub) / 9.0 + (16.0 * lnmub) / (3.0 * power_of<2>(yb)) - (16.0 * power_of<2>(lnmub)) / 3.0
                + (4.0 / (3.0 * yb3) - 4.0 / yb) * dilog1mybhalf + (-4.0 / (3.0 * yb3) + 4.0 / yb) * dilog1pybhalf
                + (20.0 / (9.0 * yb3) + (8.0 * lnmub) / (3.0 * yb3) - (4.0 + 8.0 * lnmub) / yb - (4.0 * ln2) / (3.0 * yb3) + (4.0 * ln2) /  yb) * ln1myb
                + (2.0 / (3.0 * yb3) - 2.0 / yb) * power_of<2>(ln1myb)
                + (-20.0 / (9.0 * yb3) + 4.0 / yb - (8.0 * lnmub) / (3.0 * yb3) + (8.0 * lnmub) / yb + (4.0 * ln2) / (3.0 * yb3) - (4.0 * ln2) / yb) * ln1pyb
                + (-2.0 / (3.0 * yb3) + 2.0 / yb) * power_of<2>(ln1pyb);

            return result / 81.0;
        }

        complex<double> f29ctQs(const CharmLoopsParameters & clp)
        {
            return - 6.0 * f19ctQs(clp);
        }

        complex<double> f29ctQc(const CharmLoopsParameters & clp)
        {
            const complex<double> muhat = clp.muhat;
            const complex<double> xa    = clp.x_a;
            const complex<double> yc    = clp.y_c;

            const complex<double> yc2 = power_of<2>(yc);
            const complex<double> yc3 = power_of<3>(yc);

            const complex<double> lnmuhat   = std::log(muhat);
            const complex<double> lnz       = std::log(clp.z_eps);
            const complex<double> lnxa      = std::log(xa);
            const complex<double> ln1pxa    = std::log(1.0 + xa);
            const complex<double> ln1mxa    = std::log(1.0 - xa);
            const complex<double> ln1pyc    = std::log(1.0 + yc);
            const complex<double> ln1myc    = std::log(1.0 - yc);
            const complex<double> ln1pyc2   = power_of<2>(ln1pyc);
            const complex<double> ln1myc2   = power_of<2>(ln1myc);

            const complex<double> result = 208.0 + (48.0 * myi *  M_PI) + 40.0 * power_of<2>(M_PI) - 240.0 / yc2 - (48.0 * myi * M_PI) / yc2
                - 192.0 * dilog(1.0 / 2.0) + 96.0 * dilog((1.0 - xa) / 2.0) + 96.0 * dilog((1.0 + xa) / 2.0)
                + (-24.0 / yc3 - 72.0 * yc) * dilog((1.0 - yc) / 2.0) + (24.0 / yc3 + 72.0 * yc) * dilog((1.0 + yc) / 2.0) + 96.0 * ln2 - (192.0 * myi * M_PI *  ln2)
                - (96.0 * ln2) / yc2 - 192.0 * power_of<2>(ln2) + 384.0 * lnmuhat - (192.0 * myi * M_PI *  lnmuhat) - (384.0 * lnmuhat) / yc2 - (384.0 * ln2 *  lnmuhat)
                - 192.0 * power_of<2>(lnmuhat) - 48.0 * power_of<2>(ln1mxa) - 192.0 * power_of<2>(lnxa) + ln1mxa * (-48.0 + (96.0 * myi *  M_PI) + 48.0 / yc2 + 96.0 * ln2 + 192.0 * lnmuhat
                + 192.0 * lnxa) + (-48.0 + (96.0 * myi *  M_PI) + 48.0 / yc2 + 96.0 * ln2 + 192.0 * lnmuhat) * ln1pxa - 48.0 * power_of<2>(ln1pxa)
                + lnxa * (96.0 - (192.0 * myi *  M_PI) - 96.0 / yc2 - 384.0 * ln2 - 384.0 * lnmuhat + 192.0 * ln1pxa)
                + (-12.0 / yc3 - 36.0 * yc) * ln1myc2 + (12.0 / yc3 + 36.0 * yc) * ln1pyc2 - 144.0 * lnz + (144.0 * lnz) / yc2
                + ln1pyc * (120.0 / yc3 + (24.0 * myi * M_PI) / yc3 - 144.0 / yc + 24.0 * yc + (72.0 * myi * M_PI *  yc) + (24.0 * ln2) / yc3 + (72.0 * yc *  ln2)
                    + (192.0 * lnmuhat) / yc3 - (288.0 * lnmuhat) / yc + (288.0 * yc *  lnmuhat) + (-24.0 / yc3 - 72.0 * yc) * ln1mxa
                    + (48.0 / yc3 + 144.0 * yc) * lnxa + (-24.0 / yc3 - 72.0 * yc) * ln1pxa - (72.0 * lnz) / yc3 + (144.0 * lnz) / yc - (72.0 * yc *  lnz))
                + ln1myc * (-120.0 / yc3 - (24.0 * myi * M_PI) / yc3 + 144.0 / yc - 24.0 * yc - (72.0 * myi * M_PI *  yc) - (24.0 * ln2) / yc3 - (72.0 * yc *  ln2)
                    - (192.0 * lnmuhat) / yc3 + (288.0 * lnmuhat) / yc - (288.0 * yc *  lnmuhat) + (24.0 / yc3 + 72.0 * yc) * ln1mxa
                    + (-48.0 / yc3 - 144.0 * yc) * lnxa + (24.0 / yc3 + 72.0 * yc) * ln1pxa + (72.0 * lnz) / yc3 - (144.0 * lnz) / yc + (72.0 * yc *  lnz));

            return result / 27.0;
        }

        complex<double> f29ctQb(const CharmLoopsParameters & clp)
        {
            return -6.0 * f19ctQb(clp);
        }

        // NLO two-loop function

        complex<double> f17e(const CharmLoopsParameters & clp)
        {
            return 0.0;
        }

        complex<double> f27e(const CharmLoopsParameters & clp)
        {
            return 0.0;
        }

        complex<double> f29e(const CharmLoopsParameters & clp)
        {
            const complex<double> muhat = clp.muhat;
            const complex<double> xe    = clp.x_e;
            const complex<double> ye    = clp.y_e;

            const complex<double> ye2 = power_of<2>(ye);
            const complex<double> ye3 = power_of<3>(ye);
            const complex<double> ye4 = power_of<4>(ye);

            const complex<double> lnmuhat   = std::log(muhat);
            const complex<double> lnxe      = std::log(xe);
            const complex<double> lnye      = std::log(ye);
            const complex<double> lnmye     = std::log(-ye);
            const complex<double> ln1pxe    = std::log(1.0 + xe);
            const complex<double> ln1mxe    = std::log(1.0 - xe);
            const complex<double> ln1pye    = std::log(1.0 + ye);
            const complex<double> ln1mye    = std::log(1.0 - ye);

            const complex<double> result = -636.0 - (336.0 * myi * M_PI) + 440.0 / ye2 + (288.0 * myi * M_PI) / ye2 + (-264.0 + 56.0 / ye4 - 176.0 / ye2) * dilog(1.0 / 2.0) + (576.0 - 192.0 / ye4 + 384.0 / ye2) * trilog(1.0 / 2.0)
                + (696.0 - 232.0 / ye4 + 464.0 / ye2) * zeta3 + (-576.0 + 192.0 / ye4 - 384.0 / ye2) * trilog((1.0 - ye) / 2.0) + (192.0 - 64.0 / ye4 + 128.0 / ye2) * trilog(1.0 - ye)
                + (192.0 - 64.0 / ye4 + 128.0 / ye2) * trilog(-ye) + (192.0 - 64.0 / ye4 + 128.0 / ye2) * trilog(ye) + (192.0 - 64.0 / ye4 + 128.0 / ye2) * trilog(ye / (-1.0 + ye))
                + (-288.0 + 96.0 / ye4 - 192.0 / ye2) * trilog((2.0 * ye) / (-1.0 + ye)) + (192.0 - 64.0 / ye4 + 128.0 / ye2) * trilog(ye / (1.0 + ye)) + (-288.0 + 96.0 / ye4 - 192.0 / ye2) * trilog((2.0 * ye) / (1.0 + ye))
                + (-576.0 + 192.0 / ye4 - 384.0 / ye2) * trilog((1.0 + ye) / 2.0) + (192.0 - 64.0 / ye4 + 128.0 / ye2) * trilog(1.0 + ye) - 672.0 * ln2 - 80.0 * pisqu * ln2 + (80.0 * pisqu * ln2) / (3.0 * ye4)
                + (576.0 * ln2) / ye2 - (160.0 * pisqu * ln2) / (3.0 * ye2) + 96.0 * ln2cube - (32.0 * ln2cube) / ye4 + (64.0 * ln2cube) / ye2 + 32.0 * pisqu * ln2
                - (32.0 * pisqu * ln2) / (3.0 * ye4) + (64.0 * pisqu * ln2) / (3.0 * ye2) - 672.0 * lnmuhat + (576.0 * lnmuhat) / ye2 + (336.0 - 288.0 / ye2) * ln1mxe
                + (-672.0 + 576.0 / ye2) * lnxe + (336.0 - 288.0 / ye2) * ln1pxe + (-288.0 + 96.0 / ye4 - 192.0 / ye2) * dilog(1.0 - ye) * ln1mye
                + power_of<2>(ln1mye) * (66.0 - 14.0 / ye4 + 84.0 / ye3 + 44.0 / ye2 - 216.0 / ye + 36.0 * ye - 48.0 * ln2 + (16.0 * ln2) / ye4 - (32.0 * ln2) / ye2 - 48.0 * ln2
                    + (48.0 * ln2) / (3.0 * ye4) - (96.0 * ln2) / (3.0 * ye2) + (16.0 - 16.0 / (3.0 * ye4) + 32.0 / (3.0 * ye2)) * std::log((1.0 - ye) / 8.0) + (-192.0 + 64.0 / ye4 - 128.0 / ye2) * lnye
                    + (432.0 - 144.0 / ye4 + 288.0 / ye2) * std::log((1.0 + ye) / 2.0))
                + (-288.0 + 96.0 / ye4 - 192.0 / ye2) * dilog(1.0 + ye) * ln1pye + (-32.0 * pisqu + (32.0 * pisqu) / (3.0 * ye4) - 240.0 / ye3
                    - (144.0 * myi * M_PI) / ye3 - (64.0 * pisqu) / (3.0 * ye2) + 312.0 / ye + (288.0 * myi * M_PI) / ye - 24.0 * ye - (144.0 * myi * M_PI * ye) - 132.0 * ln2 + (28.0 * ln2) / ye4 + (168.0 * ln2) / ye3
                    - (88.0 * ln2) / ye2 + (144.0 * ln2) / ye + 72.0 * ye * ln2 + 192.0 * ln2squ - (64.0 * ln2squ) / ye4 + (128.0 * ln2squ) / ye2 - 96.0 * (pisqu / 12.0 - ln2squ / 2.0)
                    + (32.0 * (pisqu / 12.0 - ln2squ / 2.0)) / ye4 - (64.0 * (pisqu / 12.0 - ln2squ / 2.0)) / ye2 + 8.0 * (pisqu + 6.0 * ln2squ) - (8.0 * (pisqu + 6.0 * ln2squ)) / (3.0 * ye4)
                    + (16.0 * (pisqu + 6.0 * ln2squ)) / (3.0 * ye2) - (288.0 * ln2) / ye3 - 288.0 * ye * ln2 - (288.0 * lnmuhat) / ye3 + (576.0 * lnmuhat) / ye - 288.0 * ye * lnmuhat
                    + (144.0 / ye3 - 288.0 / ye + 144.0 * ye) * ln1mxe + (-288.0 / ye3 + 576.0 / ye - 288.0 * ye) * lnxe + (144.0 / ye3 - 288.0 / ye + 144.0 * ye) * ln1pxe) * ln1pye
                + (66.0 - 14.0 / ye4 - 84.0 / ye3 + 44.0 / ye2 + 216.0 / ye - 36.0 * ye - 48.0 * ln2 + (16.0 * ln2) / ye4 - (32.0 * ln2) / ye2 - 48.0 * ln2 + (48.0 * ln2) / (3.0 * ye4)
                    - (96.0 * ln2) / (3.0 * ye2) + (432.0 - 144.0 / ye4 + 288.0 / ye2) * std::log((1.0 - ye) / 2.0) + (-192.0 + 64.0 / ye4 - 128.0 / ye2) * lnmye
                    + (16.0 - 16.0 / (3.0 * ye4) + 32.0 / (3.0 * ye2)) * std::log((1.0 + ye) / 8.0)) * power_of<2>(ln1pye) + dilog(ye) * (64.0 / ye3 - 192.0 / ye + (-96.0 + 32.0 / ye4 - 64.0 / ye2) * ln1mye
                    + (-192.0 + 64.0 / ye4 - 128.0 / ye2) * ln1pye) + dilog(-ye) * (-64.0 / ye3 + 192.0 / ye + (-192.0 + 64.0 / ye4 - 128.0 / ye2) * ln1mye + (-96.0 + 32.0 / ye4 - 64.0 / ye2) * ln1pye)
                + dilog((1.0 - ye) / 2.0) * (132.0 - 28.0 / ye4 + 168.0 / ye3 + 88.0 / ye2 - 432.0 / ye + 72.0 * ye + (864.0 - 288.0 / ye4 + 576.0 / ye2) * ln1mye + (288.0 - 96.0 / ye4 + 192.0 / ye2) * ln1pye)
                + dilog((1.0 + ye) / 2.0) * (132.0 - 28.0 / ye4 - 168.0 / ye3 + 88.0 / ye2 + 432.0 / ye - 72.0 * ye + (288.0 - 96.0 / ye4 + 192.0 / ye2) * ln1mye + (864.0 - 288.0 / ye4 + 576.0 / ye2) * ln1pye)
                + ln1mye * (-32.0 * pisqu + (32.0 * pisqu) / (3.0 * ye4) + 240.0 / ye3 + (144.0 * myi * M_PI) / ye3 - (64.0 * pisqu) / (3.0 * ye2) - 312.0 / ye - (288.0 * myi * M_PI) / ye + 24.0 * ye + (144.0 * myi * M_PI * ye) - 132.0 * ln2
                    + (28.0 * ln2) / ye4 - (168.0 * ln2) / ye3 - (88.0 * ln2) / ye2 - (144.0 * ln2) / ye - 72.0 * ye * ln2 + 192.0 * ln2squ - (64.0 * ln2squ) / ye4 + (128.0 * ln2squ) / ye2
                    - 96.0 * (pisqu / 12.0 - ln2squ / 2.0) + (32.0 * (pisqu / 12.0 - ln2squ / 2.0)) / ye4 - (64.0 * (pisqu / 12.0 - ln2squ / 2.0)) / ye2 + 8.0 * (pisqu + 6.0 * ln2squ) - (8.0 * (pisqu + 6.0 * ln2squ)) / (3.0 * ye4)
                    + (16.0 * (pisqu + 6.0 * ln2squ)) / (3.0 * ye2) + (288.0 * ln2) / ye3 + 288.0 * ye * ln2 + (288.0 * lnmuhat) / ye3 - (576.0 * lnmuhat) / ye + 288.0 * ye * lnmuhat
                    + (-144.0 / ye3 + 288.0 / ye - 144.0 * ye) * ln1mxe + (288.0 / ye3 - 576.0 / ye + 288.0 * ye) * lnxe + (-144.0 / ye3 + 288.0 / ye - 144.0 * ye) * ln1pxe
                    + (-576.0 * ln2 + (192.0 * ln2) / ye4 - (384.0 * ln2) / ye2) * ln1pye)
                - 576.0 * zeta3 + (192.0 * zeta3) / ye4 - (384.0 * zeta3) / ye2;

            return result / 27.0;
        }

        complex<double> f19e(const CharmLoopsParameters & clp)
        {
            const double c_F = 4.0 / 3.0; // SU(3.0) color factor

            return c_F * f29e(clp);
        }

    }
}
