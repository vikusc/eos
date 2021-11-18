/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
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

#include <eos/form-factors/k-lcdas.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <>
    struct Implementation<KaonLCDAs>
    {
        std::shared_ptr<Model> model;

        // twist 2 Gegenbauer coefficients at mu = 1 GeV
        UsedParameter a1K_0;
        UsedParameter a2K_0;

        // twist 3 parameters
        UsedParameter f3K_0;
        UsedParameter lambda3K_0;
        UsedParameter omega3K_0;

        // mass and decay constant of the pion
        UsedParameter m_K;
        UsedParameter f_K;

        // matching scales for the individual n-flavor effective QCDs
        UsedParameter _mu_c;
        UsedParameter _mu_b;
        UsedParameter _mu_t;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            a1K_0(p["K::a1@1GeV"], u),
            a2K_0(p["K::a2@1GeV"], u),
            f3K_0(p["K::f3@1GeV"], u),
            lambda3K_0(p["K::lambda3@1GeV"], u),
            omega3K_0(p["K::omega3@1GeV"], u),
            m_K(p["mass::K_u"], u),
            f_K(p["decay-constant::K_u"], u),
            _mu_c(p["QCD::mu_c"], u),
            _mu_b(p["QCD::mu_b"], u),
            _mu_t(p["QCD::mu_t"], u)
        {
        }

        inline double c_rge(const double & _mu) const
        {
            /*
             * RGE coefficient, basically
             *
             *     (alpha_s/alpha_s_0)^(1/beta_0),
             *
             * with matching between the individual n-flavor QCDs.
             */

            double mu = _mu, alpha_s_mu = model->alpha_s(mu);
            double mu_0 = 1.0, alpha_s_0 = model->alpha_s(mu_0);

            if (mu < _mu_c)
                return std::pow(alpha_s_mu / alpha_s_0, 1.0 / QCD::beta_function_nf_3[0]);

            double alpha_s_c = model->alpha_s(_mu_c);
            double result = std::pow(alpha_s_c / alpha_s_0, 1.0 / QCD::beta_function_nf_3[0]);

            if (mu < _mu_b)
                return result * std::pow(alpha_s_mu / alpha_s_c, 1.0 / QCD::beta_function_nf_4[0]);

            double alpha_s_b = model->alpha_s(_mu_b);
            result *= std::pow(alpha_s_b / alpha_s_c, 1.0 / QCD::beta_function_nf_4[0]);

            if (mu < _mu_t)
                return result * std::pow(alpha_s_mu / alpha_s_b, 1.0 / QCD::beta_function_nf_5[0]);

            throw InternalError("Implementation<KaonLCDAs>: RGE coefficient must not be evolved above mu_t = " + stringify(_mu_t()));
        }

        inline double a1K(const double & mu) const
        {
            return a1K_0 * std::pow(c_rge(mu), 50.0 / 9.0);
        }

        inline double a2K(const double & mu) const
        {
            return a2K_0 * std::pow(c_rge(mu), 364.0 / 45.0);
        }

        inline double muK(const double & mu) const
        {
            return m_K * m_K / model->m_s_msbar(mu);
        }

        double f3K(const double & mu) const
        {
            const double c_rge = this->c_rge(mu);
            const double mu_0  = 1.0; // initial state is fixed at 1 GeV
            const double m_s_0 = model->m_s_msbar(mu_0);

            return f3K_0 * std::pow(c_rge, 55.0 / 9.0)
                + 2.0 / 19.0 * (std::pow(c_rge, 4.0)        - std::pow(c_rge, 55.0 / 9.0)) * f_K * m_s_0
                + 6.0 / 65.0 * (std::pow(c_rge, 55.0 / 9.0) - std::pow(c_rge, 68.0 / 9.0)) * f_K * m_s_0 * a1K_0;
        }

        double omega3K(const double & mu) const
        {
            const double c_rge = this->c_rge(mu);
            const double mu_0  = 1.0; // initial state is fixed at 1 GeV
            const double m_s_0 = model->m_s_msbar(mu_0);

            return f3K_0 * omega3K_0 * std::pow(c_rge, 105.0 / 9.0)
                + 1.0 / 170.0 * (std::pow(c_rge, 4.0)        - std::pow(c_rge, 104.0 / 9.0)) * f_K * m_s_0
                + 1.0 /  10.0 * (std::pow(c_rge, 68.0 / 9.0) - std::pow(c_rge, 104.0 / 9.0)) * f_K * m_s_0 * a1K_0
                + 2.0 /  15.0 * (std::pow(c_rge, 86.0 / 9.0) - std::pow(c_rge, 104.0 / 9.0)) * f_K * m_s_0 * a2K_0;
        }

        double lambda3K(const double & mu) const
        {
            const double c_rge = this->c_rge(mu);
            const double mu_0  = 1.0; // initial state is fixed at 1 GeV
            const double m_s_0 = model->m_s_msbar(mu_0);

            return f3K_0 * lambda3K_0 * std::pow(c_rge, 139.0 / 18.0)
                + 14.0 / 67.0 * (std::pow(c_rge, 4.0)        - std::pow(c_rge, 139.0 / 18.0)) * f_K * m_s_0
                + 14.0 /  5.0 * (std::pow(c_rge, 68.0 / 9.0) - std::pow(c_rge, 139.0 / 18.0)) * f_K * m_s_0 * a1K_0
                - 4.0  / 11.0 * (std::pow(c_rge, 86.0 / 9.0) - std::pow(c_rge, 139.0 / 18.0)) * f_K * m_s_0 * a2K_0;
        }

        inline double eta3K(const double & mu) const
        {
            return f3K(mu) / (f_K() * muK(mu));
        }

    };

    KaonLCDAs::KaonLCDAs(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<KaonLCDAs>(new Implementation<KaonLCDAs>(p, o, *this))
    {
    }

    KaonLCDAs::~KaonLCDAs()
    {
    }

    double
    KaonLCDAs::a1(const double & mu) const
    {
        return _imp->a1K(mu);
    }

    double
    KaonLCDAs::a2(const double & mu) const
    {
        return _imp->a2K(mu);
    }

    double
    KaonLCDAs::mu(const double & mu) const
    {
        return _imp->muK(mu);
    }

    double
    KaonLCDAs::f3(const double & mu) const
    {
        return _imp->f3K(mu);
    }

    double
    KaonLCDAs::eta3(const double & mu) const
    {
        return _imp->eta3K(mu);
    }

    double
    KaonLCDAs::lambda3(const double & mu) const
    {
        return _imp->lambda3K(mu);
    }

    double
    KaonLCDAs::omega3(const double & mu) const
    {
        return _imp->omega3K(mu);
    }

    double
    KaonLCDAs::phi(const double & u, const double & mu) const
    {
        // Gegenbauer polynomials C_n^(3/2)
        const double x = 2.0 * u - 1.0, x2 = x * x;
        const double c1 = 3.0 * x;
        const double c2 = (15.0 * x2 - 3.0) / 2.0;

        return 6.0 * u * (1.0 - u) * (1.0 + _imp->a1K(mu) * c1 + _imp->a2K(mu) * c2);
    }

    double
    KaonLCDAs::phi3p(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s = _imp->model->m_s_msbar(mu);

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>(m_s / _imp->m_K); // EOM constaints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = rhopK;                        // identical in the limit m_q -> 0
        const double eta3K    = _imp->eta3K(mu);
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);

        // Gegenbauer polynomials C_n^(1/2)
        const double x = 2.0 * u - 1.0, x2 = x * x, x3 = x2 * x, x4 = x2 * x2;
        const double c1 = x;
        const double c2 = (3.0 * x2 - 1.0) / 2.0;
        const double c3 = (5.0 * x3 - 3.0 * x) / 2.0;
        const double c4 = (35.0 * x4 - 30.0 * x2 + 3.0) / 8.0;

        return 1.0 + 3.0 * rhopK * (1.0 + 6.0 * a2K) - 9.0 * rhomK * a1K
            + c1 * (27.0 / 2.0 * rhopK * a1K - rhomK * (3.0 / 2.0 + 27.0 * a2K))
            + c2 * (30.0 * eta3K + 15.0 * rhopK * a2K - 3.0 * rhomK * a1K)
            + c3 * (10.0 * eta3K * lambda3K - 9.0 / 2.0 * rhomK * a2K)
            + c4 * (-3.0 * eta3K * omega3K)
            + 3.0 / 2.0 * (rhopK + rhomK) * (1.0 - 3.0 * a1K + 6.0 * a2K) * std::log(u)
            + 3.0 / 2.0 * (rhopK - rhomK) * (1.0 + 3.0 * a1K + 6.0 * a2K) * std::log(1.0 - u);
    }

    double
    KaonLCDAs::phi3s(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s = _imp->model->m_s_msbar(mu);

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>(m_s / _imp->m_K); // EOM constaints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = rhopK;                        // identical in the limit m_q -> 0
        const double eta3K    = _imp->eta3K(mu);
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);

        // Gegenbauer polynomials C_n^(1/2)
        const double x = 2.0 * u - 1.0, x2 = x * x, x3 = x2 * x, x4 = x2 * x2;
        const double c1 = x;
        const double c2 = (3.0 * x2 - 1.0) / 2.0;
        const double c3 = (5.0 * x3 - 3.0 * x) / 2.0;
        const double c4 = (35.0 * x4 - 30.0 * x2 + 3.0) / 8.0;

        const double ubar = 1.0 - u;

        return 6.0 * u * ubar * (
            1.0 + 3.0 * rhopK * (1.0 + 6.0 * a2K) - 9.0 * rhomK * a1K
            + c1 * (27.0 / 2.0 * rhopK * a1K - rhomK * (3.0 / 2.0 + 27.0 * a2K))
            + c2 * (30.0 * eta3K + 15.0 * rhopK * a2K - 3.0 * rhomK * a1K)
            + c3 * (10.0 * eta3K * lambda3K - 9.0 / 2.0 * rhomK * a2K)
            + c4 * (-3.0 * eta3K * omega3K)
            + 3.0 / 2.0 * (rhopK + rhomK) * (1.0 - 3.0 * a1K + 6.0 * a2K) * std::log(u)
            + 3.0 / 2.0 * (rhopK - rhomK) * (1.0 + 3.0 * a1K + 6.0 * a2K) * std::log(1.0 - u)
        );
    }

    double
    KaonLCDAs::phi3s_d1(const double & u, const double & mu) const
    {
        // strange quark mass
        const double m_s = _imp->model->m_s_msbar(mu);

        // Twist 2 Gegenbauer coefficients
        const double a1K = _imp->a1K(mu);
        const double a2K = _imp->a2K(mu);

        // Twist 3 coefficients
        const double rhopK    = power_of<2>(m_s / _imp->m_K); // EOM constaints, cf. [BBL:2006A], cf. eq. (3.12)
        const double rhomK    = rhopK;                        // identical in the limit m_q -> 0
        const double eta3K    = _imp->eta3K(mu);
        const double omega3K  = _imp->omega3K(mu);
        const double lambda3K = _imp->lambda3K(mu);

        // Gegenbauer polynomials C_n^(1/2)
        const double x = 2.0 * u - 1.0, x2 = x * x, x3 = x2 * x, x4 = x2 * x2;
        const double c1  = x;
        const double c2  = (3.0 * x2 - 1.0) / 2.0;
        const double c3  = (5.0 * x3 - 3.0 * x) / 2.0;
        const double c4  = (35.0 * x4 - 30.0 * x2 + 3.0) / 8.0;
        const double c1p = 4.0 / 3.0;
        const double c2p = 32.0 / 9.0 * x;
        const double c3p = (224.0 * x2 - 48.0) / 27.0;
        const double c4p = (4480.0 * x3 - 2016.0 * x) / 243.0;

        const double ubar = 1.0 - u;

        return 6.0 * (1.0 - 2.0 * u) * ( // term 1: d/du of the prefactor
            1.0 + 3.0 * rhopK * (1.0 + 6.0 * a2K) - 9.0 * rhomK * a1K
            + c1 * (27.0 / 2.0 * rhopK * a1K - rhomK * (3.0 / 2.0 + 27.0 * a2K))
            + c2 * (30.0 * eta3K + 15.0 * rhopK * a2K - 3.0 * rhomK * a1K)
            + c3 * (10.0 * eta3K * lambda3K - 9.0 / 2.0 * rhomK * a2K)
            + c4 * (-3.0 * eta3K * omega3K)
            + 3.0 / 2.0 * (rhopK + rhomK) * (1.0 - 3.0 * a1K + 6.0 * a2K) * std::log(u)
            + 3.0 / 2.0 * (rhopK - rhomK) * (1.0 + 3.0 * a1K + 6.0 * a2K) * std::log(1.0 - u)
        ) + 6.0 * u * (1.0 - u) * ( // term 2: d/du of the Gegenbauer polynomials
            0.0
            + c1p * (27.0 / 2.0 * rhopK * a1K - rhomK * (3.0 / 2.0 + 27.0 * a2K))
            + c2p * (30.0 * eta3K + 15.0 * rhopK * a2K - 3.0 * rhomK * a1K)
            + c3p * (10.0 * eta3K * lambda3K - 9.0 / 2.0 * rhomK * a2K)
            + c4p * (-3.0 * eta3K * omega3K)
        );
    }

    Diagnostics
    KaonLCDAs::diagnostics() const
    {
        Diagnostics results;

        results.add(Diagnostics::Entry{ _imp->c_rge(1.0), "RGE coefficient C(mu = 1.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(2.0), "RGE coefficient C(mu = 2.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(3.0), "RGE coefficient C(mu = 3.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(4.0), "RGE coefficient C(mu = 4.0 GeV)" });
        results.add(Diagnostics::Entry{ _imp->c_rge(5.0), "RGE coefficient C(mu = 5.0 GeV)" });

        return results;
    }
}
