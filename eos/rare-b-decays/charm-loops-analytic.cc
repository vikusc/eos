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
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/maths/power-of.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/multiplepolylog_li22.hh>
#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/rare-b-decays/long-distance.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/log.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <complex>
#include <cstring>
#include <vector>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_dilog.h>

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

    namespace charmLoopsAnalytic
    {
        // Loop variables from [AGV:2019A]
        complex<double> s(const complex<double> q2, const complex<double> m_b)
        {
            return q2 / power_of<2>(m_b);
        }

        complex<double> z(const complex<double> m_c, const complex<double> m_b)
        {
            return power_of<2>(m_c) / power_of<2>(m_b);
        }

        complex<double> xa(const complex<double> m_c, const complex<double> m_b)
        {
            if (z(m_c, m_b) == 0.25)
            {
                throw InternalError("m_c is not half of m_b");
            }
            else{
                return 1.0 / std::sqrt(1 - 4.0 * z(m_c, m_b));#
            }
        }

        complex<double> xc(const complex<double> m_c, const complex<double> m_b)
        {
            return xa(m_c, m_b);
        }

        complex<double> xe(const complex<double> m_c, const complex<double> m_b)
        {
            return  xa(m_c, m_b);
        }

        complex<double> xb(const complex<double> m_c, const complex<double> m_b)
        {
            return std::sqrt(4.0 * z(m_c, m_b)) - std::sqrt(4.0 * z(m_c, m_b) - 1.0);
        }

        complex<double> xd(const complex<double> m_c, const complex<double> m_b)
        {
            return xb(m_c, m_b);
        }

        complex<double> ya(const complex<double> q2, const complex<double> m_c, const complex<double> m_b)
        {
            if (1.0 - s(q2, m_b) == 0.0)
            {
                return 0.0;
            }
            else
            {
                return 1.0 / std::sqrt(1.0 - 4.0 * z(m_c, m_b) / (1.0 - s(q2, m_b)));
            }
        }

        complex<double> yb(const complex<double> q2, const complex<double> m_b)
        {
            if (s(q2, m_b) == 0.0)
            {
                return 0.0;
            }
            else
            {
                return 1.0 / std::sqrt(1.0 - 4.0 / s(q2, m_b));
            }
        }

        complex<double> yc(const complex<double> q2, const complex<double> m_c)
        {
            if (q2 == 0.0)
            {
                return 0.0;
            }
            else
            {
                return 1.0 / std::sqrt(1.0 - 4.0 * power_of<2>(mc_2) / q2);
            }
        }

        complex<double> yd(const complex<double> q2, const complex<double> m_c)
        {
            return yc(q2, m_c);
        }

        complex<double> ye(const complex<double> q2, const complex<double> m_c)
        {
            return yc(q2, m_c);
        }

        // Functions depending on Heaviside step function

        // Triangle function T(a,b;x) from [FTW:2016A]
        double T(const complex<double> a, const complex<double> b, const complex<double> x)
        {
            const complex<double> amb = a-b;
            const complex<double> xconj = std::conj(x);
            const complex<double> bconj = std::conj(b);
            const double denom = (xconj * amb).imag();

            if (denom == 0.0)
            {
                return 0.0;
            }
            else
            {
                const complex<double> arg1 = (xconj * a).imag() / denom;
                const complex<double> arg2 = 1.0 - (xconj * a).imag() / denom;
                const complex<double> arg3 = - 1.0 - (bconj * a).imag() / denom;

                if ( arg1 < 0.0 || arg2 < 0.0 || arg1 > 0.0)
                {
                    return 0.0;
                }
                if (arg1 == 0.0 || arg2 == 0.0 || arg3 == 0.0)
                {
                    throw InternalError("Ill-defined Theta(0)");
                }
                else
                {
                    return 1.0;
                }
            }
        }

        // r function from [FTW:2016A]
        double r(const complex<double> a, const complex<double> b)
        {
            const complex<double> aconjb = std::conj(a) * b;
            const double denom = aconjb.imag();
            const double num = power_of<2>(std::abs(a)) * b.imag() - power_of<2>(std::abs(b)) * a.imag()

            if (denom == 0)
            {
                throw InternalError("0 in denominator");
            }
            else
            {
                return num / denom;
            }
        }

        double H1(const complex<double> a, const complex<double> b)
        {
            const complex<double> aconjb = std::conj(a) * b;

            if (aconjb.imag() == 0)
            {
                return 0.0;
            }
            else
            {
                double arg = std::min(1, power_of<2>(std::abs(a)) * b.imag() / aconjb.imag()) - r(a,b);

                if (r < 0.0 || arg < 0.0)
                {
                    return 0.0;
                }
                if (r == 0.0 || arg == 0.0)
                {
                    throw InternalError("Ill-defined Theta(0)");
                }
                else
                {
                    return 1.0;
                }
            }
        }

        // Counterterms to the two-loop functions
    }
}
