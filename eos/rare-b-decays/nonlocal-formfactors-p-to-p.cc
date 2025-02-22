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
                    polynomials(2.487638017,
                                {0.7613788603, -0.7974181049, 0.8063703241, -0.8093292634, 0.8106139436, -0.8112876795})
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
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GvDV:2020].
         */
        template <typename Process_>
        class GRvDV2021 :
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

                GRvDV2021(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToP>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),

                    re_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_0^plus}@GRvDV2021"], *this),
                    im_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_0^plus}@GRvDV2021"], *this),
                    re_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_1^plus}@GRvDV2021"], *this),
                    im_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_1^plus}@GRvDV2021"], *this),
                    re_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_2^plus}@GRvDV2021"], *this),
                    im_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_2^plus}@GRvDV2021"], *this),
                    re_alpha_3_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_3^plus}@GRvDV2021"], *this),
                    im_alpha_3_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_3^plus}@GRvDV2021"], *this),
                    re_alpha_4_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_4^plus}@GRvDV2021"], *this),
                    im_alpha_4_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_4^plus}@GRvDV2021"], *this),
                    re_alpha_5_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_5^plus}@GRvDV2021"], *this),
                    im_alpha_5_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_5^plus}@GRvDV2021"], *this),
                    re_alpha_6_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_6^plus}@GRvDV2021"], *this),
                    im_alpha_6_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_6^plus}@GRvDV2021"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_d"], *this),

                    m_P(p["mass::K_d"], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),
                    chiOPE(p["b->sccbar::chiOPE@GRvDV2021"], *this)
                {
                    this->uses(*form_factors);
                }

                ~GRvDV2021() = default;

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
                    const complex<double> phi2 = pow(m_B4 * power_of<4>(z - 1.) - 2 * m_B2 * power_of<2>(z - 1.) * (-16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) + s_0 * power_of<2>(z + 1.)) +
                                                power_of<2>(16 * m_D02 * z + m_P2 * power_of<2>(z - 1.) - s_0 * power_of<2>(z + 1.)), 0.5);//(C8)
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0*(z + 1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * power_of<2>(z + 1.) - 16. * z * m_D02, -0.5); //(C10)

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z, a-b+c+d-1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                inline complex<double> phi(const double & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    return phi(complex<double>(q2, 0.0), phi_parameters);
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, 7> alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return eos::nff_utils::P<6u>(z_Jpsi, alpha) / phi(m_Jpsi2, phi_parameters) *
                            (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, 7> alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return eos::nff_utils::P<6u>(z_psi2S, alpha) / phi(m_psi2S2, phi_parameters) *
                            (1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
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
                    const auto z       = eos::nff_utils::z(q2,              s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return eos::nff_utils::P<6u>(z, alpha) / phi(q2, phi_parameters) / blaschke_factor;
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

                    return eos::nff_utils::P<6u>(z, alpha);
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
                };

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

                    const std::array<unsigned, 4> phi_parameters = {3, 3, 2, 2};

                    return eos::nff_utils::P<6u>(z, alpha) / phi(q2, phi_parameters) / F_plus;
                }

                static NonlocalFormFactorPtr<nff::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToP>(new GRvDV2021<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const double s_0   = this->t_0();
                    const auto z1 = eos::nff_utils::z(1.0, 4.0 * power_of<2>(m_D0), s_0);
                    const std::array<complex<double>, 7> alpha{2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0};

                    results.add({ real(eos::nff_utils::P<6u>(z1, alpha)), "Re{P(q2 = 1.0, {2.0, 3.0, 4.0})}" });
                    results.add({ imag(eos::nff_utils::P<6u>(z1, alpha)), "Im{P(q2 = 1.0, {2.0, 3.0, 4.0})}" });

                    const std::array<complex<double>, 7> alpha2{
                        complex<double>(2.0,5.0),
                        complex<double>(3.0,6.0),
                        complex<double>(4.0,7.0),
                        complex<double>(0.0,8.0),
                        0., 0., 0.
                        };

                    results.add({ real(eos::nff_utils::P<6u>(z1, alpha2)),
                                "Re{P(q2 = 1.0, {(2.0,5.0), (3.0,6.0), (4.0,7.0), (0.0,8.0))}}" });
                    results.add({ imag(eos::nff_utils::P<6u>(z1, alpha2)),
                                "Im{P(q2 = 1.0, {(2.0,5.0), (3.0,6.0), (4.0,7.0), (0.0,8.0))}}" });

                    return results;

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
            std::make_pair("B->K::naive",         &nff_p_to_p::Naive::make),
            // parametrizations
            std::make_pair("B->K::GvDV2020",      &nff_p_to_p::GvDV2020<nff::BToK>::make),
            std::make_pair("B->K::GRvDV2021",     &nff_p_to_p::GRvDV2021<nff::BToK>::make),
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
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_normalized_moment_A(const double & q2) const
    {
        return real(this->_imp->nff->normalized_moment_A(q2));
    }

    template <typename Process_>
    const std::set<ReferenceName>
    NonlocalFormFactorObservable<Process_, nff::PToP>::references
    {
        "GvDV:2020A"_rn
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
