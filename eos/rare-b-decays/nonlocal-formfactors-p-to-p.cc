/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2020 Danny van Dyk
 * Copyright (c) 2020 Nico Gubernari
 * Copyright (c) 2021 Méril Reboud
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

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>

#include <map>
#include <numeric>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

namespace eos
{
    using std::abs;

    namespace nff
    {
        struct BToK
        {
            constexpr static const char * label = "B->K";
        };
        constexpr const char * BToK::label;
    }

    namespace nff_p_to_p
    {
        class Naive :
            public NonlocalFormFactor<nff::PToP>
        {
            public:
                Naive(const Parameters &, const Options &)
                {
                }

                ~Naive() = default;

                virtual complex<double> H_plus(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_plus(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> Hhat_plus(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_plus(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_plus(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> F_ratio_plus(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> P_ratio_plus(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> get_orthonormal_coefficients(const unsigned &) const
                {
                    return 0.0;
                }

                virtual double weak_bound() const
                {
                    return 0.0;
                }

                virtual double strong_bound() const
                {
                    return 0.0;
                }

                virtual double weak_bound_log_likelihood() const
                {
                    return 0.0;
                }

                virtual double strong_bound_log_likelihood() const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_A(const double &) const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<nff::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToP>(new Naive(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    return {};
                }
        };

        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GvDV:2020].
         */
        template <typename Process_>
        class GvDV2020 :
            public NonlocalFormFactor<nff::PToP>
        {
            public:
                std::shared_ptr<FormFactors<PToP>> form_factors;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_plus;
                UsedParameter im_alpha_0_plus;
                UsedParameter re_alpha_1_plus;
                UsedParameter im_alpha_1_plus;
                UsedParameter re_alpha_2_plus;
                UsedParameter im_alpha_2_plus;
                UsedParameter re_alpha_3_plus;
                UsedParameter im_alpha_3_plus;
                UsedParameter re_alpha_4_plus;
                UsedParameter im_alpha_4_plus;
                UsedParameter re_alpha_5_plus;
                UsedParameter im_alpha_5_plus;
                UsedParameter re_alpha_6_plus;
                UsedParameter im_alpha_6_plus;

                // Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_P;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                // Orthogonal polynomials on an arc of the unit circle
                const SzegoPolynomial<6u> polynomials;

                GvDV2020(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToP>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),

                    re_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_0^plus}@GvDV2020"], *this),
                    im_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_0^plus}@GvDV2020"], *this),
                    re_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_1^plus}@GvDV2020"], *this),
                    im_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_1^plus}@GvDV2020"], *this),
                    re_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_2^plus}@GvDV2020"], *this),
                    im_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_2^plus}@GvDV2020"], *this),
                    re_alpha_3_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_3^plus}@GvDV2020"], *this),
                    im_alpha_3_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_3^plus}@GvDV2020"], *this),
                    re_alpha_4_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_4^plus}@GvDV2020"], *this),
                    im_alpha_4_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_4^plus}@GvDV2020"], *this),
                    re_alpha_5_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_5^plus}@GvDV2020"], *this),
                    im_alpha_5_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_5^plus}@GvDV2020"], *this),
                    re_alpha_6_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_6^plus}@GvDV2020"], *this),
                    im_alpha_6_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_6^plus}@GvDV2020"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_d"], *this),

                    m_P(p["mass::K_d"], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),
                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this),

                    // The parameters of the polynomial expension are computed using t0 = 4.0 and
                    // the masses are set to mB = 5.279 and mK = 0.492 (same values as for local form-factors)
                    polynomials(2.48247,
                                {0.762292, -0.798241, 0.807153, -0.810097, 0.811376, -0.812046})
                {
                    this->uses(*form_factors);
                }

                ~GvDV2020() = default;

                inline complex<double> phi(const complex<double> & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d
                    // 0(P->P) aka plus          3    3    2    2
                    // perp(P->V) = par(P->V)    3    1    3    0
                    // 0(P->V) aka long          3    1    2    2

                    const double m_P2  = power_of<2>(m_P);
                    const double m_B2  = power_of<2>(m_B),  m_B4 =  power_of<4>(m_B);
                    const double m_D02 = power_of<2>(m_D0), m_D04 = power_of<4>(m_D0);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * power_of<2>(m_D0), s_0);
                    const double Q2    = this->t_s();
                    const double chi   = this->chiOPE();

                    const double a = phi_parameters[0], b = phi_parameters[1], c = phi_parameters[2], d = phi_parameters[3];

                    const double Nlambda = 4 * M_PI * pow(m_B2, 0.5 * (a - b + c + d) - 1.) * pow(2 * (4 * m_D02 - s_0) / 3 / chi, 0.5); //(C6)
                    const complex<double> phi1 = -pow(2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 - Q2 - s_0, 0.5) /
                                                (2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 + Q2 * (z - 1.) - s_0 * (z + 1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4 * power_of<4>(z - 1.) - 2 * m_B2 * power_of<2>(z - 1.) * (-16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) + s_0 * pow(z + 1., 2)) +
                                                power_of<2>(16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) - s_0 * power_of<2>(z + 1.)), 0.5);//(C8)
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0*(z+1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * power_of<2>(z + 1.) - 16. * z * m_D02, -0.5); //(C10)

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z, a-b+c+d-1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                inline complex<double> phi(const double & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    return phi(complex<double>(q2, 0.0), phi_parameters);
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, 7> & alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const auto & polynomials_at_z = polynomials(z_Jpsi);
                    const complex<double> p_at_z = std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);


                    return p_at_z / phi(m_Jpsi2, phi_parameters) * (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, 7> & alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const auto & polynomials_at_z = polynomials(z_psi2S);
                    const complex<double> p_at_z = std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_psi2S2, phi_parameters) * (1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                }


                virtual complex<double> H_plus(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 7> alpha{
                        complex<double>(re_alpha_0_plus, im_alpha_0_plus),
                        complex<double>(re_alpha_1_plus, im_alpha_1_plus),
                        complex<double>(re_alpha_2_plus, im_alpha_2_plus),
                        complex<double>(re_alpha_3_plus, im_alpha_3_plus),
                        complex<double>(re_alpha_4_plus, im_alpha_4_plus),
                        complex<double>(re_alpha_5_plus, im_alpha_5_plus),
                        complex<double>(re_alpha_6_plus, im_alpha_6_plus),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                   s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    const auto & polynomials_at_z = polynomials(z);
                    const complex<double> p_at_z = std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_plus(const double & q2) const
                {
                    return H_plus(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    const std::array<complex<double>, 7> alpha{
                        complex<double>(re_alpha_0_plus, im_alpha_0_plus),
                        complex<double>(re_alpha_1_plus, im_alpha_1_plus),
                        complex<double>(re_alpha_2_plus, im_alpha_2_plus),
                        complex<double>(re_alpha_3_plus, im_alpha_3_plus),
                        complex<double>(re_alpha_4_plus, im_alpha_4_plus),
                        complex<double>(re_alpha_5_plus, im_alpha_5_plus),
                        complex<double>(re_alpha_6_plus, im_alpha_6_plus),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    const auto & polynomials_at_z = polynomials(z);

                    return std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));
                }


                virtual complex<double> H_plus_residue_jpsi() const
                {
                    const std::array<complex<double>, 7> alpha{
                        complex<double>(re_alpha_0_plus, im_alpha_0_plus),
                        complex<double>(re_alpha_1_plus, im_alpha_1_plus),
                        complex<double>(re_alpha_2_plus, im_alpha_2_plus),
                        complex<double>(re_alpha_3_plus, im_alpha_3_plus),
                        complex<double>(re_alpha_4_plus, im_alpha_4_plus),
                        complex<double>(re_alpha_5_plus, im_alpha_5_plus),
                        complex<double>(re_alpha_6_plus, im_alpha_6_plus),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return H_residue_jpsi(phi_parameters, alpha);
                }

                virtual complex<double> H_plus_residue_psi2s() const
                {
                    const std::array<complex<double>, 7> alpha{
                        complex<double>(re_alpha_0_plus, im_alpha_0_plus),
                        complex<double>(re_alpha_1_plus, im_alpha_1_plus),
                        complex<double>(re_alpha_2_plus, im_alpha_2_plus),
                        complex<double>(re_alpha_3_plus, im_alpha_3_plus),
                        complex<double>(re_alpha_4_plus, im_alpha_4_plus),
                        complex<double>(re_alpha_5_plus, im_alpha_5_plus),
                        complex<double>(re_alpha_6_plus, im_alpha_6_plus),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return H_residue_psi2s(phi_parameters, alpha);
                }

                virtual complex<double> normalized_moment_A(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_plus(const complex<double> & q2) const
                {
                    const complex<double> F_plus = form_factors->f_p(q2);

                    return H_plus(q2) / F_plus;
                }

                virtual complex<double> ratio_plus(const double & q2) const
                {
                    return ratio_plus(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_plus(const complex<double> & q2) const
                {
                    return form_factors->f_t(q2) * q2 / m_B() / (m_B + m_P) / form_factors->f_p(q2);
                }

                virtual complex<double> P_ratio_plus(const double & q2) const
                {
                    const std::array<complex<double>, 7> alpha{
                        complex<double>(re_alpha_0_plus, im_alpha_0_plus),
                        complex<double>(re_alpha_1_plus, im_alpha_1_plus),
                        complex<double>(re_alpha_2_plus, im_alpha_2_plus),
                        complex<double>(re_alpha_3_plus, im_alpha_3_plus),
                        complex<double>(re_alpha_4_plus, im_alpha_4_plus),
                        complex<double>(re_alpha_5_plus, im_alpha_5_plus),
                        complex<double>(re_alpha_6_plus, im_alpha_6_plus),
                    };

                    const complex<double> F_plus = form_factors->f_p(q2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    const auto & polynomials_at_z = polynomials(z);
                    const complex<double> p_at_z = std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return p_at_z / phi(q2, phi_parameters) / F_plus;
                }

                virtual complex<double> get_orthonormal_coefficients(const unsigned & i) const
                {
                    const std::array<complex<double>, 7> alpha{
                        complex<double>(re_alpha_0_plus, im_alpha_0_plus),
                        complex<double>(re_alpha_1_plus, im_alpha_1_plus),
                        complex<double>(re_alpha_2_plus, im_alpha_2_plus),
                        complex<double>(re_alpha_3_plus, im_alpha_3_plus),
                        complex<double>(re_alpha_4_plus, im_alpha_4_plus),
                        complex<double>(re_alpha_5_plus, im_alpha_5_plus),
                        complex<double>(re_alpha_6_plus, im_alpha_6_plus),
                    };

                    return alpha[i];
                }

                virtual double weak_bound() const
                {
                    return 0.;
                }

                virtual double strong_bound() const
                {
                    return 0.;
                }

                virtual double weak_bound_log_likelihood() const
                {
                    return 0.0;
                }

                virtual double strong_bound_log_likelihood() const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<nff::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToP>(new GvDV2020<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2}; //plus polarization

                    results.add({ real(1./this->phi(0.0, phi_parameters)), "Re{1/phi_+(q2 = 0.0)}" });
                    results.add({ imag(1./this->phi(0.0, phi_parameters)), "Im{1/phi_+(q2 = 0.0)}" });
                    results.add({ real(this->phi(16.0, phi_parameters)), "Re{phi_+(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phi_parameters)), "Im{phi_+(q2 = 16.0)}" });

                    const double s_0   = this->t_0();
                    const auto z1 = eos::nff_utils::z(1.0, 4.0 * power_of<2>(m_D0), s_0);
                    const std::array<complex<double>, 6> alpha = {2.0, 3.0, 4.0, 5.0, 0.0};

                    SzegoPolynomial<5u> p{
                        1.854590436,
                        {0.8627241729, -0.8863992518, 0.8911183885, -0.8926229103, 0.8932870967}
                    };

                    const auto & polynomials_at_z = p(z1);

                    const complex<double> p_at_z = std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    results.add({ std::real(p_at_z), "Re{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, {2.0, 3.0, 4.0, 5.0})}" });
                    results.add({ std::imag(p_at_z), "Im{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, {2.0, 3.0, 4.0, 5.0})}" });

                    return results;
                }
        };


        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GvDV:2020],
         * but using an ad hoc Lagrange polynomial.
         */
        template <typename Process_>
        class GRvDV2022order5 :
            public NonlocalFormFactor<nff::PToP>
        {
            public:
                std::shared_ptr<FormFactors<PToP>> form_factors;

                // Polynomial expansion parameters
                UsedParameter re_at_m7_plus;
                UsedParameter im_at_m7_plus;
                UsedParameter re_at_m5_plus;
                UsedParameter im_at_m5_plus;
                UsedParameter re_at_m3_plus;
                UsedParameter im_at_m3_plus;
                UsedParameter re_at_m1_plus;
                UsedParameter im_at_m1_plus;
                UsedParameter abs_at_Jpsi_plus;
                UsedParameter arg_at_Jpsi_plus;
                UsedParameter abs_at_psi2S_plus;
                UsedParameter arg_at_psi2S_plus;

                // Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_P;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;
                UsedParameter bound;
                UsedParameter bound_uncertainty;

                // Lagrange interpolating polynomial
                const static unsigned interpolation_order = 5;
                const LagrangePolynomial<interpolation_order> lagrange;

                // Orthogonal polynomials on an arc of the unit circle used for the computation of dispersive bounds
                const SzegoPolynomial<interpolation_order> orthonormal_polynomials;

                GRvDV2022order5(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToP>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),

                    re_at_m7_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m7_plus@GRvDV2022"], *this),
                    im_at_m7_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m7_plus@GRvDV2022"], *this),
                    re_at_m5_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m5_plus@GRvDV2022"], *this),
                    im_at_m5_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m5_plus@GRvDV2022"], *this),
                    re_at_m3_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m3_plus@GRvDV2022"], *this),
                    im_at_m3_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m3_plus@GRvDV2022"], *this),
                    re_at_m1_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m1_plus@GRvDV2022"], *this),
                    im_at_m1_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m1_plus@GRvDV2022"], *this),
                    abs_at_Jpsi_plus(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_Jpsi_plus@GRvDV2022"], *this),
                    arg_at_Jpsi_plus(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_Jpsi_plus@GRvDV2022"], *this),
                    abs_at_psi2S_plus(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_psi2S_plus@GRvDV2022"], *this),
                    arg_at_psi2S_plus(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_psi2S_plus@GRvDV2022"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_d"], *this),

                    m_P(p["mass::K_d"], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),
                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this),
                    bound(p["b->sccbar::bound@GvDV2020"], *this),
                    bound_uncertainty(p["b->sccbar::bound_uncertainty@GvDV2020"], *this),

                    lagrange({eos::nff_utils::z(-7.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(-5.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(-3.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(-1.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(power_of<2>(m_Jpsi),  4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(power_of<2>(m_psi2S), 4.0 * power_of<2>(m_D0), t_0())}),

                    // The parameters of the polynomial expension are computed using t0 = 4.0 and
                    // the masses are set to mB = 5.279 and mK = 0.492 (same values as for local form-factors)
                    orthonormal_polynomials(2.48247,
                                            {0.762292, -0.798241, 0.807153, -0.810097, 0.811376})
                {
                    this->uses(*form_factors);
                }

                ~GRvDV2022order5() = default;

                inline complex<double> phi(const complex<double> & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d
                    // 0(P->P) aka plus          3    3    2    2
                    // perp(P->V) = par(P->V)    3    1    3    0
                    // 0(P->V) aka long          3    1    2    2

                    const double m_P2  = power_of<2>(m_P);
                    const double m_B2  = power_of<2>(m_B),  m_B4 =  power_of<4>(m_B);
                    const double m_D02 = power_of<2>(m_D0), m_D04 = power_of<4>(m_D0);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * m_D02, s_0);
                    const double Q2    = this->t_s();
                    const double chi   = this->chiOPE();

                    const double a = phi_parameters[0], b = phi_parameters[1], c = phi_parameters[2], d = phi_parameters[3];

                    const double Nlambda = 4 * M_PI * pow(m_B2, 0.5 * (a - b + c + d) - 1.) * pow(2 * (4 * m_D02 - s_0) / 3 / chi, 0.5); //(C6)
                    const complex<double> phi1 = -pow(2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 - Q2 - s_0, 0.5) /
                                                (2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 + Q2 * (z - 1.) - s_0 * (z + 1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4 * power_of<4>(z - 1.) - 2 * m_B2 * power_of<2>(z - 1.) * (-16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) + s_0 * pow(z + 1., 2)) +
                                                power_of<2>(16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) - s_0 * power_of<2>(z + 1.)), 0.5);//(C8)
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0*(z+1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * power_of<2>(z + 1.) - 16. * z * m_D02, -0.5); //(C10)

                    return Nlambda * pow(1. + z, 0.5) * pow(1. - z, a - b + c + d - 1.5) * pow(phi1, a) * pow(phi2, 0.5 * b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                inline complex<double> phi(const double & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    return phi(complex<double>(q2, 0.0), phi_parameters);
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, interpolation_order + 1> & interpolation_values) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> p_at_z = lagrange(interpolation_values, z_Jpsi);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);


                    return p_at_z / phi(m_Jpsi2, phi_parameters) * (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, interpolation_order + 1> & interpolation_values) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> p_at_z = lagrange(interpolation_values, z_psi2S);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_psi2S2, phi_parameters) * (1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                }


                virtual complex<double> H_plus(const complex<double> & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    const complex<double> p_at_z = lagrange(interpolation_values, z);

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_plus(const double & q2) const
                {
                    return H_plus(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return lagrange(interpolation_values, z);
                }

                virtual complex<double> H_plus_residue_jpsi() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return H_residue_jpsi(phi_parameters, interpolation_values);
                }

                virtual complex<double> H_plus_residue_psi2s() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return H_residue_psi2s(phi_parameters, interpolation_values);
                }

                virtual complex<double> normalized_moment_A(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_plus(const complex<double> & q2) const
                {
                    const complex<double> F_plus = form_factors->f_p(q2);

                    return H_plus(q2) / F_plus;
                }

                virtual complex<double> ratio_plus(const double & q2) const
                {
                    return ratio_plus(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_plus(const complex<double> & q2) const
                {
                    return form_factors->f_t(q2) * q2 / m_B() / (m_B + m_P) / form_factors->f_p(q2);
                }

                virtual complex<double> P_ratio_plus(const double & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);
                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};
                    const complex<double> F_plus = form_factors->f_p(q2);

                    const complex<double> p_at_z = lagrange(interpolation_values, z);

                    return p_at_z / phi(q2, phi_parameters) / F_plus;
                }

                inline std::pair<gsl_vector *, gsl_vector *> orthonormal_coefficients() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    std::array<complex<double>, interpolation_order + 1> dL = lagrange.get_coefficients(interpolation_values);

                    // Split array of derivatives to real and imaginary parts
                    gsl_vector * dL_real_part = gsl_vector_calloc(interpolation_order + 1);
                    gsl_vector * dL_imag_part = gsl_vector_calloc(interpolation_order + 1);

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        gsl_vector_set(dL_real_part, i, real(dL[i]));
                        gsl_vector_set(dL_imag_part, i, imag(dL[i]));
                    }

                    const gsl_matrix * coefficient_matrix = orthonormal_polynomials.coefficient_matrix();

                    // Solve the system by computing (coefficient_matrix)^(-1) . dL_real_part and idem for imag
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, dL_real_part);
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, dL_imag_part);

                    return std::make_pair(dL_real_part, dL_imag_part);
                }

                virtual complex<double> get_orthonormal_coefficients(const unsigned & i) const
                {
                    auto coefficients = orthonormal_coefficients();

                    return complex<double>(gsl_vector_get(coefficients.first,  i),
                                           gsl_vector_get(coefficients.second, i));
                }

                virtual double weak_bound() const
                {
                    auto coefficients = orthonormal_coefficients();

                    double largest_absolute_coeff = 0.0, coeff;

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        coeff =   power_of<2>(gsl_vector_get(coefficients.first,  i))
                                + power_of<2>(gsl_vector_get(coefficients.second, i));
                        if (coeff > largest_absolute_coeff)
                        {
                            largest_absolute_coeff = coeff;
                        }
                    }

                    return largest_absolute_coeff;
                }

                virtual double strong_bound() const
                {
                    auto coefficients = orthonormal_coefficients();

                    double coefficient_sum = 0.0;

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        coefficient_sum +=   power_of<2>(gsl_vector_get(coefficients.first,  i))
                                           + power_of<2>(gsl_vector_get(coefficients.second, i));
                    }

                    return coefficient_sum;
                }

                virtual double weak_bound_log_likelihood() const
                {
                    const double saturation = weak_bound();
                    if (saturation < this->bound)
                    {
                        return 0.;
                    }
                    else
                    {
                        // Halfnormal constraint
                        return -0.5 * power_of<2>( (saturation - this->bound) / this->bound_uncertainty );
                    }
                }

                virtual double strong_bound_log_likelihood() const
                {
                    const double saturation = strong_bound();
                    if (saturation < this->bound)
                    {
                        return 0.;
                    }
                    else
                    {
                        // Halfnormal constraint
                        return -0.5 * power_of<2>( (saturation - this->bound) / this->bound_uncertainty );
                    }
                }

                static NonlocalFormFactorPtr<nff::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToP>(new GRvDV2022order5<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    return results;
                }
        };


        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GvDV:2020],
         * but using an ad hoc Lagrange polynomial.
         */
        template <typename Process_, std::array s_constraints_ >
        class GRvDV2022order3_c :
            public NonlocalFormFactor<nff::PToP>
        {
            public:
                std::shared_ptr<FormFactors<PToP>> form_factors;

                // Polynomial expansion parameters
                UsedParameter re_at_s0_plus;
                UsedParameter im_at_s0_plus;
                UsedParameter re_at_s1_plus;
                UsedParameter im_at_s1_plus;
                UsedParameter re_at_s2_plus;
                UsedParameter im_at_s2_plus;
                UsedParameter re_at_s3_plus;
                UsedParameter im_at_s3_plus;

                // Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_P;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;
                UsedParameter bound;
                UsedParameter bound_uncertainty;

                // Lagrange interpolating polynomial
                const static unsigned interpolation_order = 3;
                const LagrangePolynomial<interpolation_order> lagrange;

                // Orthogonal polynomials on an arc of the unit circle used for the computation of dispersive bounds
                const SzegoPolynomial<interpolation_order> orthonormal_polynomials;

                GRvDV2022order3_c(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToP>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),

                    re_at_s0_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_s0_plus@GRvDV2022"], *this),
                    im_at_s0_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_s0_plus@GRvDV2022"], *this),
                    re_at_s1_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_s1_plus@GRvDV2022"], *this),
                    im_at_s1_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_s1_plus@GRvDV2022"], *this),
                    re_at_s2_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_s2_plus@GRvDV2022"], *this),
                    im_at_s2_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_s2_plus@GRvDV2022"], *this),
                    re_at_s3_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_s3_plus@GRvDV2022"], *this),
                    im_at_s3_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_s3_plus@GRvDV2022"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_d"], *this),

                    m_P(p["mass::K_d"], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),
                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this),
                    bound(p["b->sccbar::bound@GvDV2020"], *this),
                    bound_uncertainty(p["b->sccbar::bound_uncertainty@GvDV2020"], *this),

                    lagrange({eos::nff_utils::z(complex<double>(s_constraints_[0], s_constraints_[1]), 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(complex<double>(s_constraints_[2], s_constraints_[3]), 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(complex<double>(s_constraints_[4], s_constraints_[5]), 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(complex<double>(s_constraints_[6], s_constraints_[7]), 4.0 * power_of<2>(m_D0), t_0())}),

                    // The parameters of the polynomial expension are computed using t0 = 4.0 and
                    // the masses are set to mB = 5.279 and mK = 0.492 (same values as for local form-factors)
                    orthonormal_polynomials(2.48247,
                                            {0.762292, -0.798241, 0.807153})
                {
                    this->uses(*form_factors);
                }

                ~GRvDV2022order3_c() = default;

                inline complex<double> phi(const complex<double> & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d
                    // 0(P->P) aka plus          3    3    2    2
                    // perp(P->V) = par(P->V)    3    1    3    0
                    // 0(P->V) aka long          3    1    2    2

                    const double m_P2  = power_of<2>(m_P);
                    const double m_B2  = power_of<2>(m_B),  m_B4 =  power_of<4>(m_B);
                    const double m_D02 = power_of<2>(m_D0), m_D04 = power_of<4>(m_D0);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * m_D02, s_0);
                    const double Q2    = this->t_s();
                    const double chi   = this->chiOPE();

                    const double a = phi_parameters[0], b = phi_parameters[1], c = phi_parameters[2], d = phi_parameters[3];

                    const double Nlambda = 4 * M_PI * pow(m_B2, 0.5 * (a - b + c + d) - 1.) * pow(2 * (4 * m_D02 - s_0) / 3 / chi, 0.5); //(C6)
                    const complex<double> phi1 = -pow(2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 - Q2 - s_0, 0.5) /
                                                (2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 + Q2 * (z - 1.) - s_0 * (z + 1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4 * power_of<4>(z - 1.) - 2 * m_B2 * power_of<2>(z - 1.) * (-16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) + s_0 * pow(z + 1., 2)) +
                                                power_of<2>(16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) - s_0 * power_of<2>(z + 1.)), 0.5);//(C8)
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0*(z+1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * power_of<2>(z + 1.) - 16. * z * m_D02, -0.5); //(C10)

                    return Nlambda * pow(1. + z, 0.5) * pow(1. - z, a - b + c + d - 1.5) * pow(phi1, a) * pow(phi2, 0.5 * b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                inline complex<double> phi(const double & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    return phi(complex<double>(q2, 0.0), phi_parameters);
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, interpolation_order + 1> & interpolation_values) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> p_at_z = lagrange(interpolation_values, z_Jpsi);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);


                    return p_at_z / phi(m_Jpsi2, phi_parameters) * (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, interpolation_order + 1> & interpolation_values) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> p_at_z = lagrange(interpolation_values, z_psi2S);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_psi2S2, phi_parameters) * (1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                }


                virtual complex<double> H_plus(const complex<double> & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_s0_plus, im_at_s0_plus),
                        complex<double>(re_at_s1_plus, im_at_s1_plus),
                        complex<double>(re_at_s2_plus, im_at_s2_plus),
                        complex<double>(re_at_s3_plus, im_at_s3_plus)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    const complex<double> p_at_z = lagrange(interpolation_values, z);

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_plus(const double & q2) const
                {
                    return H_plus(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_s0_plus, im_at_s0_plus),
                        complex<double>(re_at_s1_plus, im_at_s1_plus),
                        complex<double>(re_at_s2_plus, im_at_s2_plus),
                        complex<double>(re_at_s3_plus, im_at_s3_plus)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return lagrange(interpolation_values, z);
                }

                virtual complex<double> H_plus_residue_jpsi() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_s0_plus, im_at_s0_plus),
                        complex<double>(re_at_s1_plus, im_at_s1_plus),
                        complex<double>(re_at_s2_plus, im_at_s2_plus),
                        complex<double>(re_at_s3_plus, im_at_s3_plus)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return H_residue_jpsi(phi_parameters, interpolation_values);
                }

                virtual complex<double> H_plus_residue_psi2s() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_s0_plus, im_at_s0_plus),
                        complex<double>(re_at_s1_plus, im_at_s1_plus),
                        complex<double>(re_at_s2_plus, im_at_s2_plus),
                        complex<double>(re_at_s3_plus, im_at_s3_plus)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return H_residue_psi2s(phi_parameters, interpolation_values);
                }

                virtual complex<double> normalized_moment_A(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_plus(const complex<double> & q2) const
                {
                    const complex<double> F_plus = form_factors->f_p(q2);

                    return H_plus(q2) / F_plus;
                }

                virtual complex<double> ratio_plus(const double & q2) const
                {
                    return ratio_plus(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_plus(const complex<double> & q2) const
                {
                    return form_factors->f_t(q2) * q2 / m_B() / (m_B + m_P) / form_factors->f_p(q2);
                }

                virtual complex<double> P_ratio_plus(const double & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_s0_plus, im_at_s0_plus),
                        complex<double>(re_at_s1_plus, im_at_s1_plus),
                        complex<double>(re_at_s2_plus, im_at_s2_plus),
                        complex<double>(re_at_s3_plus, im_at_s3_plus)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);
                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};
                    const complex<double> F_plus = form_factors->f_p(q2);

                    const complex<double> p_at_z = lagrange(interpolation_values, z);

                    return p_at_z / phi(q2, phi_parameters) / F_plus;
                }

                inline std::pair<gsl_vector *, gsl_vector *> orthonormal_coefficients() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_s0_plus, im_at_s0_plus),
                        complex<double>(re_at_s1_plus, im_at_s1_plus),
                        complex<double>(re_at_s2_plus, im_at_s2_plus),
                        complex<double>(re_at_s3_plus, im_at_s3_plus)
                    };

                    std::array<complex<double>, interpolation_order + 1> dL = lagrange.get_coefficients(interpolation_values);

                    // Split array of derivatives to real and imaginary parts
                    gsl_vector * dL_real_part = gsl_vector_calloc(interpolation_order + 1);
                    gsl_vector * dL_imag_part = gsl_vector_calloc(interpolation_order + 1);

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        gsl_vector_set(dL_real_part, i, real(dL[i]));
                        gsl_vector_set(dL_imag_part, i, imag(dL[i]));
                    }

                    const gsl_matrix * coefficient_matrix = orthonormal_polynomials.coefficient_matrix();

                    // Solve the system by computing (coefficient_matrix)^(-1) . dL_real_part and idem for imag
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, dL_real_part);
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, dL_imag_part);

                    return std::make_pair(dL_real_part, dL_imag_part);
                }

                virtual complex<double> get_orthonormal_coefficients(const unsigned & i) const
                {
                    auto coefficients = orthonormal_coefficients();

                    return complex<double>(gsl_vector_get(coefficients.first,  i),
                                           gsl_vector_get(coefficients.second, i));
                }

                virtual double weak_bound() const
                {
                    auto coefficients = orthonormal_coefficients();

                    double largest_absolute_coeff = 0.0, coeff;

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        coeff =   power_of<2>(gsl_vector_get(coefficients.first,  i))
                                + power_of<2>(gsl_vector_get(coefficients.second, i));
                        if (coeff > largest_absolute_coeff)
                        {
                            largest_absolute_coeff = coeff;
                        }
                    }

                    return largest_absolute_coeff;
                }

                virtual double strong_bound() const
                {
                    auto coefficients = orthonormal_coefficients();

                    double coefficient_sum = 0.0;

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        coefficient_sum +=   power_of<2>(gsl_vector_get(coefficients.first,  i))
                                           + power_of<2>(gsl_vector_get(coefficients.second, i));
                    }

                    return coefficient_sum;
                }

                virtual double weak_bound_log_likelihood() const
                {
                    const double saturation = weak_bound();
                    if (saturation < this->bound)
                    {
                        return 0.;
                    }
                    else
                    {
                        // Halfnormal constraint
                        return -0.5 * power_of<2>( (saturation - this->bound) / this->bound_uncertainty );
                    }
                }

                virtual double strong_bound_log_likelihood() const
                {
                    const double saturation = strong_bound();
                    if (saturation < this->bound)
                    {
                        return 0.;
                    }
                    else
                    {
                        // Halfnormal constraint
                        return -0.5 * power_of<2>( (saturation - this->bound) / this->bound_uncertainty );
                    }
                }

                static NonlocalFormFactorPtr<nff::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToP>(new GRvDV2022order3_c<Process_, s_constraints_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    return results;
                }
        };


        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GvDV:2020],
         * but using an ad hoc Lagrange polynomial.
         */
        template <typename Process_>
        class GRvDV2022order6 :
            public NonlocalFormFactor<nff::PToP>
        {
            public:
                std::shared_ptr<FormFactors<PToP>> form_factors;

                // Polynomial expansion parameters
                UsedParameter re_at_m7_plus;
                UsedParameter im_at_m7_plus;
                UsedParameter re_at_m5_plus;
                UsedParameter im_at_m5_plus;
                UsedParameter re_at_m3_plus;
                UsedParameter im_at_m3_plus;
                UsedParameter re_at_m1_plus;
                UsedParameter im_at_m1_plus;
                UsedParameter re_at_t0_plus;
                UsedParameter im_at_t0_plus;
                UsedParameter abs_at_Jpsi_plus;
                UsedParameter arg_at_Jpsi_plus;
                UsedParameter abs_at_psi2S_plus;
                UsedParameter arg_at_psi2S_plus;

                // Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_P;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;
                UsedParameter bound;
                UsedParameter bound_uncertainty;

                // Lagrange interpolating polynomial
                const static unsigned interpolation_order = 6;
                const LagrangePolynomial<interpolation_order> lagrange;

                // Orthogonal polynomials on an arc of the unit circle used for the computation of dispersive bounds
                const SzegoPolynomial<interpolation_order> orthonormal_polynomials;

                GRvDV2022order6(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToP>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),

                    re_at_m7_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m7_plus@GRvDV2022"], *this),
                    im_at_m7_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m7_plus@GRvDV2022"], *this),
                    re_at_m5_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m5_plus@GRvDV2022"], *this),
                    im_at_m5_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m5_plus@GRvDV2022"], *this),
                    re_at_m3_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m3_plus@GRvDV2022"], *this),
                    im_at_m3_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m3_plus@GRvDV2022"], *this),
                    re_at_m1_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m1_plus@GRvDV2022"], *this),
                    im_at_m1_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m1_plus@GRvDV2022"], *this),
                    re_at_t0_plus(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_t0_plus@GRvDV2022"], *this),
                    im_at_t0_plus(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_t0_plus@GRvDV2022"], *this),
                    abs_at_Jpsi_plus(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_Jpsi_plus@GRvDV2022"], *this),
                    arg_at_Jpsi_plus(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_Jpsi_plus@GRvDV2022"], *this),
                    abs_at_psi2S_plus(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_psi2S_plus@GRvDV2022"], *this),
                    arg_at_psi2S_plus(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_psi2S_plus@GRvDV2022"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_d"], *this),

                    m_P(p["mass::K_d"], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),
                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this),
                    bound(p["b->sccbar::bound@GvDV2020"], *this),
                    bound_uncertainty(p["b->sccbar::bound_uncertainty@GvDV2020"], *this),

                    lagrange({eos::nff_utils::z(-7.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(-5.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(-3.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(-1.0, 4.0 * power_of<2>(m_D0), t_0()),
                              0., // z(t_0) = 0. by construction
                              eos::nff_utils::z(power_of<2>(m_Jpsi),  4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(power_of<2>(m_psi2S), 4.0 * power_of<2>(m_D0), t_0())}),

                    // The parameters of the polynomial expension are computed using t0 = 4.0 and
                    // the masses are set to mB = 5.279 and mK = 0.492 (same values as for local form-factors)
                    orthonormal_polynomials(2.48247,
                                {0.762292, -0.798241, 0.807153, -0.810097, 0.811376, -0.812046})
                {
                    this->uses(*form_factors);
                }

                ~GRvDV2022order6() = default;

                inline complex<double> phi(const complex<double> & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d
                    // 0(P->P) aka plus          3    3    2    2
                    // perp(P->V) = par(P->V)    3    1    3    0
                    // 0(P->V) aka long          3    1    2    2

                    const double m_P2  = power_of<2>(m_P);
                    const double m_B2  = power_of<2>(m_B),  m_B4 =  power_of<4>(m_B);
                    const double m_D02 = power_of<2>(m_D0), m_D04 = power_of<4>(m_D0);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * m_D02, s_0);
                    const double Q2    = this->t_s();
                    const double chi   = this->chiOPE();

                    const double a = phi_parameters[0], b = phi_parameters[1], c = phi_parameters[2], d = phi_parameters[3];

                    const double Nlambda = 4 * M_PI * pow(m_B2, 0.5 * (a - b + c + d) - 1.) * pow(2 * (4 * m_D02 - s_0) / 3 / chi, 0.5); //(C6)
                    const complex<double> phi1 = -pow(2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 - Q2 - s_0, 0.5) /
                                                (2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 + Q2 * (z - 1.) - s_0 * (z + 1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4 * power_of<4>(z - 1.) - 2 * m_B2 * power_of<2>(z - 1.) * (-16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) + s_0 * pow(z + 1., 2)) +
                                                power_of<2>(16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) - s_0 * power_of<2>(z + 1.)), 0.5);//(C8)
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0*(z+1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * power_of<2>(z + 1.) - 16. * z * m_D02, -0.5); //(C10)

                    return Nlambda * pow(1. + z, 0.5) * pow(1. - z, a - b + c + d - 1.5) * pow(phi1, a) * pow(phi2, 0.5 * b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                inline complex<double> phi(const double & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    return phi(complex<double>(q2, 0.0), phi_parameters);
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, interpolation_order + 1> & interpolation_values) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> p_at_z = lagrange(interpolation_values, z_Jpsi);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);


                    return p_at_z / phi(m_Jpsi2, phi_parameters) * (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, interpolation_order + 1> & interpolation_values) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> p_at_z = lagrange(interpolation_values, z_psi2S);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_psi2S2, phi_parameters) * (1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                }


                virtual complex<double> H_plus(const complex<double> & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        complex<double>(re_at_t0_plus, im_at_t0_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    const complex<double> p_at_z = lagrange(interpolation_values, z);

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_plus(const double & q2) const
                {
                    return H_plus(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        complex<double>(re_at_t0_plus, im_at_t0_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return lagrange(interpolation_values, z);
                }

                virtual complex<double> H_plus_residue_jpsi() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        complex<double>(re_at_t0_plus, im_at_t0_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return H_residue_jpsi(phi_parameters, interpolation_values);
                }

                virtual complex<double> H_plus_residue_psi2s() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        complex<double>(re_at_t0_plus, im_at_t0_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return H_residue_psi2s(phi_parameters, interpolation_values);
                }

                virtual complex<double> normalized_moment_A(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_plus(const complex<double> & q2) const
                {
                    const complex<double> F_plus = form_factors->f_p(q2);

                    return H_plus(q2) / F_plus;
                }

                virtual complex<double> ratio_plus(const double & q2) const
                {
                    return ratio_plus(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_plus(const complex<double> & q2) const
                {
                    return form_factors->f_t(q2) * q2 / m_B() / (m_B + m_P) / form_factors->f_p(q2);
                }

                virtual complex<double> P_ratio_plus(const double & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        complex<double>(re_at_t0_plus, im_at_t0_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);
                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};
                    const complex<double> F_plus = form_factors->f_p(q2);

                    const complex<double> p_at_z = lagrange(interpolation_values, z);

                    return p_at_z / phi(q2, phi_parameters) / F_plus;
                }

                inline std::pair<gsl_vector *, gsl_vector *> orthonormal_coefficients() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_plus, im_at_m7_plus),
                        complex<double>(re_at_m5_plus, im_at_m5_plus),
                        complex<double>(re_at_m3_plus, im_at_m3_plus),
                        complex<double>(re_at_m1_plus, im_at_m1_plus),
                        complex<double>(re_at_t0_plus, im_at_t0_plus),
                        polar<double>(abs_at_Jpsi_plus, arg_at_Jpsi_plus),
                        polar<double>(abs_at_psi2S_plus, arg_at_psi2S_plus)
                    };

                    std::array<complex<double>, interpolation_order + 1> dL = lagrange.get_coefficients(interpolation_values);

                    // Split array of derivatives to real and imaginary parts
                    gsl_vector * dL_real_part = gsl_vector_calloc(interpolation_order + 1);
                    gsl_vector * dL_imag_part = gsl_vector_calloc(interpolation_order + 1);

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        gsl_vector_set(dL_real_part, i, real(dL[i]));
                        gsl_vector_set(dL_imag_part, i, imag(dL[i]));
                    }

                    const gsl_matrix * coefficient_matrix = orthonormal_polynomials.coefficient_matrix();

                    // Solve the system by computing (coefficient_matrix)^(-1) . dL_real_part and idem for imag
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, dL_real_part);
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, dL_imag_part);

                    return std::make_pair(dL_real_part, dL_imag_part);
                }

                virtual complex<double> get_orthonormal_coefficients(const unsigned & i) const
                {
                    auto coefficients = orthonormal_coefficients();

                    return complex<double>(gsl_vector_get(coefficients.first,  i),
                                           gsl_vector_get(coefficients.second, i));
                }

                virtual double weak_bound() const
                {
                    auto coefficients = orthonormal_coefficients();

                    double largest_absolute_coeff = 0.0, coeff;

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        coeff =  power_of<2>(gsl_vector_get(coefficients.first,  i))
                               + power_of<2>(gsl_vector_get(coefficients.second, i));
                        if (coeff > largest_absolute_coeff)
                        {
                            largest_absolute_coeff = coeff;
                        }
                    }

                    return largest_absolute_coeff;
                }

                virtual double strong_bound() const
                {
                    auto coefficients = orthonormal_coefficients();

                    double coefficient_sum = 0.0;

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        coefficient_sum +=  power_of<2>(gsl_vector_get(coefficients.first,  i))
                                          + power_of<2>(gsl_vector_get(coefficients.second, i));
                    }

                    return coefficient_sum;
                }

                virtual double weak_bound_log_likelihood() const
                {
                    const double saturation = weak_bound();
                    if (saturation < this->bound)
                    {
                        return 0.;
                    }
                    else
                    {
                        // Halfnormal constraint
                        return -0.5 * power_of<2>( (saturation - this->bound) / this->bound_uncertainty );
                    }
                }

                virtual double strong_bound_log_likelihood() const
                {
                    const double saturation = strong_bound();
                    if (saturation < this->bound)
                    {
                        return 0.;
                    }
                    else
                    {
                        // Halfnormal constraint
                        return -0.5 * power_of<2>( (saturation - this->bound) / this->bound_uncertainty );
                    }
                }

                static NonlocalFormFactorPtr<nff::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToP>(new GRvDV2022order6<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    return results;
                }
        };

    }

    NonlocalFormFactorPtr<nff::PToP>
    NonlocalFormFactor<nff::PToP>::make(const QualifiedName & name, const Parameters & p, const Options & o)
    {
        typedef QualifiedName KeyType;
        typedef std::function<NonlocalFormFactorPtr<nff::PToP> (const Parameters &, const Options &)> ValueType;
        std::map<KeyType, ValueType> entries
        {
            // trivial
            std::make_pair("B->K::naive",                &nff_p_to_p::Naive::make),
            // parametrizations
            std::make_pair("B->K::GvDV2020",             &nff_p_to_p::GvDV2020<nff::BToK>::make),
            std::make_pair("B->K::GRvDV2022order5",      &nff_p_to_p::GRvDV2022order5<nff::BToK>::make),
            std::make_pair("B->K::GRvDV2022order6",      &nff_p_to_p::GRvDV2022order6<nff::BToK>::make),
            std::make_pair("B->K::GRvDV2022order3_c",    &nff_p_to_p::GRvDV2022order3_c<nff::BToK, {-7.0, 0.0, -5.0, 0.5, -3.0, -0.5, -1.0, 0.0}>::make),
        };

        auto i = entries.find(name);
        if (entries.end() == i)
        {
            return NonlocalFormFactorPtr<nff::PToP>(nullptr);
        }

        return i->second(p, o);
    }

    template <typename Process_>
    struct Implementation<NonlocalFormFactorObservable<Process_, nff::PToP>>
    {
        NameOption opt_formfactor;
        NonlocalFormFactorPtr<nff::PToP> nff;
        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_formfactor(o, "nonlocal-formfactor", qnp::Name("GvDV2020")),
            nff(NonlocalFormFactor<nff::PToP>::make(QualifiedName(qnp::Prefix(Process_::label), opt_formfactor.value()), p, o))
        {
            u.uses(*nff);
        }

        ~Implementation() = default;
    };

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nff::PToP>::NonlocalFormFactorObservable(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<NonlocalFormFactorObservable<Process_, nff::PToP>>(new Implementation<NonlocalFormFactorObservable<Process_, nff::PToP>>(p, o, *this))
    {
    }

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nff::PToP>::~NonlocalFormFactorObservable() = default;

    template <typename Process_>
    const std::vector<OptionSpecification>
    Implementation<NonlocalFormFactorObservable<Process_, nff::PToP>>::options
    {
    };

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_H_plus(const double & q2) const
    {
        return real(this->_imp->nff->H_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::im_H_plus(const double & q2) const
    {
        return imag(this->_imp->nff->H_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::abs_H_plus(const double & q2) const
    {
        return abs(this->_imp->nff->H_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_Hhat_plus(const double & q2) const
    {
        return real(this->_imp->nff->Hhat_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::im_Hhat_plus(const double & q2) const
    {
        return imag(this->_imp->nff->Hhat_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::abs_Hhat_plus(const double & q2) const
    {
        return abs(this->_imp->nff->Hhat_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_ratio_plus(const double & q2) const
    {
        return real(this->_imp->nff->ratio_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::im_ratio_plus(const double & q2) const
    {
        return imag(this->_imp->nff->ratio_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::abs_ratio_plus(const double & q2) const
    {
        return abs(this->_imp->nff->ratio_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::abs_P_ratio_plus(const double & q2) const
    {
        return abs(this->_imp->nff->P_ratio_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_ratio_plus_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->ratio_plus(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::im_ratio_plus_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->ratio_plus(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_F_ratio_plus_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->F_ratio_plus(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::im_F_ratio_plus_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->F_ratio_plus(complex<double>(re_q2, im_q2)));
    }


    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::get_real_alpha(const unsigned & i) const
    {
        return real(this->_imp->nff->get_orthonormal_coefficients(i));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::get_imag_alpha(const unsigned & i) const
    {
        return imag(this->_imp->nff->get_orthonormal_coefficients(i));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::weak_bound() const
    {
        return this->_imp->nff->weak_bound();
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::weak_bound_log_likelihood() const
    {
        return this->_imp->nff->weak_bound_log_likelihood();
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::strong_bound() const
    {
        return this->_imp->nff->strong_bound();
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::strong_bound_log_likelihood() const
    {
        return this->_imp->nff->strong_bound_log_likelihood();
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_normalized_moment_A(const double & q2) const
    {
        return real(this->_imp->nff->normalized_moment_A(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::im_normalized_moment_A(const double & q2) const
    {
        return imag(this->_imp->nff->normalized_moment_A(q2));
    }

    template <typename Process_>
    const std::set<ReferenceName>
    NonlocalFormFactorObservable<Process_, nff::PToP>::references
    {
        "GvDV:2020A"_rn,
        "GRvDV:2022A"_rn,
    };

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    NonlocalFormFactorObservable<Process_, nff::PToP>::begin_options()
    {
        return Implementation<NonlocalFormFactorObservable<Process_, nff::PToP>>::options.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    NonlocalFormFactorObservable<Process_, nff::PToP>::end_options()
    {
        return Implementation<NonlocalFormFactorObservable<Process_, nff::PToP>>::options.cend();
    }

    template class NonlocalFormFactorObservable<nff::BToK, nff::PToP>;
}
