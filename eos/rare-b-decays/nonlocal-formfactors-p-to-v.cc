/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2019 Danny van Dyk
 * Copyright (c) 2019 Nico Gubernari
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
#include <eos/maths/complex.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>

#include <map>
#include <numeric>

namespace eos
{
    using std::abs;

    namespace nff
    {
        struct BToKstar
        {
            static constexpr const char * label = "B->K^*";
        };
        constexpr const char * BToKstar::label;

        struct BsToPhi
        {
            static constexpr const char * label = "B_s->phi";
        };
        constexpr const char * BsToPhi::label;
    }

    namespace nff_p_to_v
    {
        class Naive :
            public NonlocalFormFactor<nff::PToV>
        {
            public:
                Naive(const Parameters &, const Options &)
                {
                }

                ~Naive() = default;

                virtual complex<double> H_perp(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_para(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_long(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_perp(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_para(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_long(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> Hhat_perp(const double &) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_para(const double &) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_long(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_perp(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_para(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_long(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_perp(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_para(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_long(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> F_ratio_perp(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> F_ratio_para(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> F_ratio_long(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V1(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V2(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V23(const double &) const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<nff::PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToV>(new Naive(p, o));
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
            public NonlocalFormFactor<nff::PToV>
        {
            private:
                std::shared_ptr<FormFactors<PToV>> form_factors;

                // spectator quark option
                SwitchOption opt_q;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_perp;
                UsedParameter im_alpha_0_perp;
                UsedParameter re_alpha_1_perp;
                UsedParameter im_alpha_1_perp;
                UsedParameter re_alpha_2_perp;
                UsedParameter im_alpha_2_perp;
                UsedParameter re_alpha_3_perp;
                UsedParameter im_alpha_3_perp;
                UsedParameter re_alpha_4_perp;
                UsedParameter im_alpha_4_perp;
                UsedParameter re_alpha_5_perp;
                UsedParameter im_alpha_5_perp;

                UsedParameter re_alpha_0_para;
                UsedParameter im_alpha_0_para;
                UsedParameter re_alpha_1_para;
                UsedParameter im_alpha_1_para;
                UsedParameter re_alpha_2_para;
                UsedParameter im_alpha_2_para;
                UsedParameter re_alpha_3_para;
                UsedParameter im_alpha_3_para;
                UsedParameter re_alpha_4_para;
                UsedParameter im_alpha_4_para;
                UsedParameter re_alpha_5_para;
                UsedParameter im_alpha_5_para;

                UsedParameter re_alpha_0_long;
                UsedParameter im_alpha_0_long;
                UsedParameter re_alpha_1_long;
                UsedParameter im_alpha_1_long;
                UsedParameter re_alpha_2_long;
                UsedParameter im_alpha_2_long;
                UsedParameter re_alpha_3_long;
                UsedParameter im_alpha_3_long;
                UsedParameter re_alpha_4_long;
                UsedParameter im_alpha_4_long;
                UsedParameter re_alpha_5_long;
                UsedParameter im_alpha_5_long;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_V;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                // Orthogonal polynomials on an arc of the unit circle
                std::shared_ptr<SzegoPolynomial<5u>> polynomials;

                std::string _final_state() const
                {
                    switch (opt_q.value()[0])
                    {
                        case 's':
                            return "phi";
                            break;

                        default:
                            return "K_d^*";
                    }
                }

            public:
                GvDV2020(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToV>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),
                    opt_q(o, "q", { "u", "d", "s" }),

                    re_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_0^perp}@GvDV2020"], *this),
                    im_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_0^perp}@GvDV2020"], *this),
                    re_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_1^perp}@GvDV2020"], *this),
                    im_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_1^perp}@GvDV2020"], *this),
                    re_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_2^perp}@GvDV2020"], *this),
                    im_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_2^perp}@GvDV2020"], *this),
                    re_alpha_3_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_3^perp}@GvDV2020"], *this),
                    im_alpha_3_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_3^perp}@GvDV2020"], *this),
                    re_alpha_4_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_4^perp}@GvDV2020"], *this),
                    im_alpha_4_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_4^perp}@GvDV2020"], *this),
                    re_alpha_5_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_5^perp}@GvDV2020"], *this),
                    im_alpha_5_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_5^perp}@GvDV2020"], *this),

                    re_alpha_0_para(p[stringify(Process_::label) + "ccbar::Re{alpha_0^para}@GvDV2020"], *this),
                    im_alpha_0_para(p[stringify(Process_::label) + "ccbar::Im{alpha_0^para}@GvDV2020"], *this),
                    re_alpha_1_para(p[stringify(Process_::label) + "ccbar::Re{alpha_1^para}@GvDV2020"], *this),
                    im_alpha_1_para(p[stringify(Process_::label) + "ccbar::Im{alpha_1^para}@GvDV2020"], *this),
                    re_alpha_2_para(p[stringify(Process_::label) + "ccbar::Re{alpha_2^para}@GvDV2020"], *this),
                    im_alpha_2_para(p[stringify(Process_::label) + "ccbar::Im{alpha_2^para}@GvDV2020"], *this),
                    re_alpha_3_para(p[stringify(Process_::label) + "ccbar::Re{alpha_3^para}@GvDV2020"], *this),
                    im_alpha_3_para(p[stringify(Process_::label) + "ccbar::Im{alpha_3^para}@GvDV2020"], *this),
                    re_alpha_4_para(p[stringify(Process_::label) + "ccbar::Re{alpha_4^para}@GvDV2020"], *this),
                    im_alpha_4_para(p[stringify(Process_::label) + "ccbar::Im{alpha_4^para}@GvDV2020"], *this),
                    re_alpha_5_para(p[stringify(Process_::label) + "ccbar::Re{alpha_5^para}@GvDV2020"], *this),
                    im_alpha_5_para(p[stringify(Process_::label) + "ccbar::Im{alpha_5^para}@GvDV2020"], *this),

                    re_alpha_0_long(p[stringify(Process_::label) + "ccbar::Re{alpha_0^long}@GvDV2020"], *this),
                    im_alpha_0_long(p[stringify(Process_::label) + "ccbar::Im{alpha_0^long}@GvDV2020"], *this),
                    re_alpha_1_long(p[stringify(Process_::label) + "ccbar::Re{alpha_1^long}@GvDV2020"], *this),
                    im_alpha_1_long(p[stringify(Process_::label) + "ccbar::Im{alpha_1^long}@GvDV2020"], *this),
                    re_alpha_2_long(p[stringify(Process_::label) + "ccbar::Re{alpha_2^long}@GvDV2020"], *this),
                    im_alpha_2_long(p[stringify(Process_::label) + "ccbar::Im{alpha_2^long}@GvDV2020"], *this),
                    re_alpha_3_long(p[stringify(Process_::label) + "ccbar::Re{alpha_3^long}@GvDV2020"], *this),
                    im_alpha_3_long(p[stringify(Process_::label) + "ccbar::Im{alpha_3^long}@GvDV2020"], *this),
                    re_alpha_4_long(p[stringify(Process_::label) + "ccbar::Re{alpha_4^long}@GvDV2020"], *this),
                    im_alpha_4_long(p[stringify(Process_::label) + "ccbar::Im{alpha_4^long}@GvDV2020"], *this),
                    re_alpha_5_long(p[stringify(Process_::label) + "ccbar::Re{alpha_5^long}@GvDV2020"], *this),
                    im_alpha_5_long(p[stringify(Process_::label) + "ccbar::Im{alpha_5^long}@GvDV2020"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_" + opt_q.value()], *this),

                    m_V(p["mass::" + _final_state()], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),

                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this),

                    // The parameters of the polynomial expension are computed using t0 = 4.0 and
                    // the masses are set to the same values as for local form-factors
                    polynomials(PolynomialsFactory::create(opt_q.value()))
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

                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B),  m_B4 =  power_of<4>(m_B);
                    const double m_D02 = power_of<2>(m_D0), m_D04 = power_of<4>(m_D0);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * power_of<2>(m_D0), s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phi_parameters[0], b = phi_parameters[1], c = phi_parameters[2], d = phi_parameters[3];

                    const double Nlambda = 4 * M_PI * pow(m_B2, 0.5 * (a - b + c + d) - 1.) * pow(2 * (4 * m_D02 - s_0) / 3 / chi, 0.5); //(C6)
                    const complex<double> phi1 = -pow(2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 - Q2 - s_0, 0.5) /
                                                (2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 + Q2 * (z - 1.) - s_0 * (z + 1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4 * power_of<4>(z - 1.) - 2 * m_B2 * power_of<2>(z - 1.) * (-16 * m_D02 * z + m_V2 * power_of<2>(z - 1.) + s_0 * power_of<2>(z + 1.)) +
						power_of<2>(16 * m_D02 * z + m_V2 * power_of<2>(z - 1.) - s_0 * power_of<2>(z + 1.)), 0.5);//(C8)
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
                inline complex<double> H_residue_jpsi(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, 6> & alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z_Jpsi);
                    const complex<double> p_at_z = std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_Jpsi2, phi_parameters) * (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, 6> & alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z_psi2S);
                    const complex<double> p_at_z = std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_psi2S2, phi_parameters) *(1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                }

                virtual complex<double> H_perp(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                     s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),    s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S),   s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z);
                    const complex<double> p_at_z = std::inner_product(alpha_perp.begin(), alpha_perp.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_perp(const double & q2) const
                {
                    return H_perp(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_perp(const double & q2) const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z);

                    return std::inner_product(alpha_perp.begin(), alpha_perp.end(), polynomials_at_z.begin(), complex<double>(0, 0));
                }

                virtual complex<double> H_para(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                     s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),    s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S),   s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    const auto & polynomials_at_z = (*polynomials)(z);
                    const complex<double> p_at_z = std::inner_product(alpha_para.begin(), alpha_para.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_para(const double & q2) const
                {
                    return H_para(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_para(const double & q2) const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z);

                    return std::inner_product(alpha_para.begin(), alpha_para.end(), polynomials_at_z.begin(), complex<double>(0, 0));
                }

                virtual complex<double> H_long(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                     s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),    s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S),   s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 2, 2};

                    const auto & polynomials_at_z = (*polynomials)(z);
                    const complex<double> p_at_z = std::inner_product(alpha_long.begin(), alpha_long.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_long(const double & q2) const
                {
                    return H_long(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_long(const double & q2) const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z);

                    return std::inner_product(alpha_long.begin(), alpha_long.end(), polynomials_at_z.begin(), complex<double>(0, 0));
                }


                virtual complex<double> H_perp_residue_jpsi() const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_jpsi(phi_parameters, alpha_perp);
                }

                virtual complex<double> H_perp_residue_psi2s() const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_psi2s(phi_parameters, alpha_perp);
                }

                virtual complex<double> H_para_residue_jpsi() const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_jpsi(phi_parameters, alpha_para);
                }

                virtual complex<double> H_para_residue_psi2s() const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_psi2s(phi_parameters, alpha_para);
                }

                virtual complex<double> H_long_residue_jpsi() const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 2, 2};

                    return H_residue_jpsi(phi_parameters, alpha_long);
                }

                virtual complex<double> H_long_residue_psi2s() const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 2, 2};

                    return H_residue_psi2s(phi_parameters, alpha_long);
                }

                virtual complex<double> ratio_perp(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_perp = pow(2.0 * lambda, 0.5) / (m_B + m_V) / m_B() * form_factors->v(q2);

                    return H_perp(q2) / F_perp;
                }

                virtual complex<double> ratio_perp(const double & q2) const
                {
                    return ratio_perp(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_perp(const complex<double> & q2) const
                {
                    return (m_B + m_V) / m_B * form_factors->t_1(q2) / form_factors->v(q2);
                }

                virtual complex<double> ratio_para(const complex<double> & q2) const
                {
                    const complex<double> F_para = sqrt(2) * (m_B + m_V) / m_B * form_factors->a_1(q2);

                    return H_para(q2) / F_para;
                }

                virtual complex<double> ratio_para(const double & q2) const
                {
                    return ratio_para(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_para(const complex<double> & q2) const
                {
                    return (m_B - m_V) / m_B * form_factors->t_2(q2) / form_factors->a_1(q2);
                }

                virtual complex<double> ratio_long(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_long = ((m_B2 - m_V2 - q2) * power_of<2>(m_B + m_V) * form_factors->a_1(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * (m_B + m_V));

                    return H_long(q2) / F_long;
                }

                virtual complex<double> ratio_long(const double & q2) const
                {
                    return ratio_long(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_long(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_long = ((m_B2 - m_V2 - q2) * power_of<2>(m_B + m_V) * form_factors->a_1(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * (m_B + m_V));
                    const complex<double> F_T_long = q2 * ((m_B2 + 3 * m_V2 - q2) * (m_B2 - m_V2) * form_factors->t_2(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * m_B * (m_B2 - m_V2));

                    return F_T_long / F_long;
                }

                virtual complex<double> normalized_moment_V1(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V2(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V23(const double &) const
                {
                    return 0.0;
                }


                static NonlocalFormFactorPtr<nff::PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToV>(new GvDV2020<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const std::array<unsigned, 4> phi_parameters_long = {3, 1, 2, 2}; //long polarization
                    results.add({ real(1./this->phi(0.0, phi_parameters_long)), "Re{1/phi_long(q2 = 0.0)}" });
                    results.add({ imag(1./this->phi(0.0, phi_parameters_long)), "Im{1/phi_long(q2 = 0.0)}" });
                    results.add({ real(this->phi(16.0, phi_parameters_long)), "Re{phi_long(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phi_parameters_long)), "Im{phi_long(q2 = 16.0)}" });

                    const std::array<unsigned, 4> phi_parameters_perp = {3, 1, 3, 0}; //perp or para polarization
                    results.add({ real(this->phi(16.0, phi_parameters_perp)), "Re{phi_perp(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phi_parameters_perp)), "Im{phi_perp(q2 = 16.0)}" });

                    return results;
                }
        };



        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GRvDV:2021].
         */
        template <typename Process_>
        class GRvDV2021 :
            public NonlocalFormFactor<nff::PToV>
        {
            private:
                std::shared_ptr<FormFactors<PToV>> form_factors;

                // spectator quark option
                SwitchOption opt_q;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_perp;
                UsedParameter im_alpha_0_perp;
                UsedParameter re_alpha_1_perp;
                UsedParameter im_alpha_1_perp;
                UsedParameter re_alpha_2_perp;
                UsedParameter im_alpha_2_perp;
                UsedParameter re_alpha_3_perp;
                UsedParameter im_alpha_3_perp;
                UsedParameter re_alpha_4_perp;
                UsedParameter im_alpha_4_perp;
                UsedParameter re_alpha_5_perp;
                UsedParameter im_alpha_5_perp;

                UsedParameter re_alpha_0_para;
                UsedParameter im_alpha_0_para;
                UsedParameter re_alpha_1_para;
                UsedParameter im_alpha_1_para;
                UsedParameter re_alpha_2_para;
                UsedParameter im_alpha_2_para;
                UsedParameter re_alpha_3_para;
                UsedParameter im_alpha_3_para;
                UsedParameter re_alpha_4_para;
                UsedParameter im_alpha_4_para;
                UsedParameter re_alpha_5_para;
                UsedParameter im_alpha_5_para;

                UsedParameter re_alpha_0_long;
                UsedParameter im_alpha_0_long;
                UsedParameter re_alpha_1_long;
                UsedParameter im_alpha_1_long;
                UsedParameter re_alpha_2_long;
                UsedParameter im_alpha_2_long;
                UsedParameter re_alpha_3_long;
                UsedParameter im_alpha_3_long;
                UsedParameter re_alpha_4_long;
                UsedParameter im_alpha_4_long;
                UsedParameter re_alpha_5_long;
                UsedParameter im_alpha_5_long;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;
                UsedParameter m_Bsst;

                // final state meson parameters
                UsedParameter m_V;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                std::string _final_state() const
                {
                    switch (opt_q.value()[0])
                    {
                        case 's':
                            return "phi";
                            break;

                        default:
                            return "K_d^*";
                    }
                }

            public:
                GRvDV2021(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToV>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),
                    opt_q(o, "q", { "u", "d", "s" }),

                    re_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_0^perp}@GRvDV2021"], *this),
                    im_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_0^perp}@GRvDV2021"], *this),
                    re_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_1^perp}@GRvDV2021"], *this),
                    im_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_1^perp}@GRvDV2021"], *this),
                    re_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_2^perp}@GRvDV2021"], *this),
                    im_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_2^perp}@GRvDV2021"], *this),
                    re_alpha_3_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_3^perp}@GRvDV2021"], *this),
                    im_alpha_3_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_3^perp}@GRvDV2021"], *this),
                    re_alpha_4_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_4^perp}@GRvDV2021"], *this),
                    im_alpha_4_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_4^perp}@GRvDV2021"], *this),
                    re_alpha_5_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_5^perp}@GRvDV2021"], *this),
                    im_alpha_5_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_5^perp}@GRvDV2021"], *this),

                    re_alpha_0_para(p[stringify(Process_::label) + "ccbar::Re{alpha_0^para}@GRvDV2021"], *this),
                    im_alpha_0_para(p[stringify(Process_::label) + "ccbar::Im{alpha_0^para}@GRvDV2021"], *this),
                    re_alpha_1_para(p[stringify(Process_::label) + "ccbar::Re{alpha_1^para}@GRvDV2021"], *this),
                    im_alpha_1_para(p[stringify(Process_::label) + "ccbar::Im{alpha_1^para}@GRvDV2021"], *this),
                    re_alpha_2_para(p[stringify(Process_::label) + "ccbar::Re{alpha_2^para}@GRvDV2021"], *this),
                    im_alpha_2_para(p[stringify(Process_::label) + "ccbar::Im{alpha_2^para}@GRvDV2021"], *this),
                    re_alpha_3_para(p[stringify(Process_::label) + "ccbar::Re{alpha_3^para}@GRvDV2021"], *this),
                    im_alpha_3_para(p[stringify(Process_::label) + "ccbar::Im{alpha_3^para}@GRvDV2021"], *this),
                    re_alpha_4_para(p[stringify(Process_::label) + "ccbar::Re{alpha_4^para}@GRvDV2021"], *this),
                    im_alpha_4_para(p[stringify(Process_::label) + "ccbar::Im{alpha_4^para}@GRvDV2021"], *this),
                    re_alpha_5_para(p[stringify(Process_::label) + "ccbar::Re{alpha_5^para}@GRvDV2021"], *this),
                    im_alpha_5_para(p[stringify(Process_::label) + "ccbar::Im{alpha_5^para}@GRvDV2021"], *this),

                    re_alpha_0_long(p[stringify(Process_::label) + "ccbar::Re{alpha_0^long}@GRvDV2021"], *this),
                    im_alpha_0_long(p[stringify(Process_::label) + "ccbar::Im{alpha_0^long}@GRvDV2021"], *this),
                    re_alpha_1_long(p[stringify(Process_::label) + "ccbar::Re{alpha_1^long}@GRvDV2021"], *this),
                    im_alpha_1_long(p[stringify(Process_::label) + "ccbar::Im{alpha_1^long}@GRvDV2021"], *this),
                    re_alpha_2_long(p[stringify(Process_::label) + "ccbar::Re{alpha_2^long}@GRvDV2021"], *this),
                    im_alpha_2_long(p[stringify(Process_::label) + "ccbar::Im{alpha_2^long}@GRvDV2021"], *this),
                    re_alpha_3_long(p[stringify(Process_::label) + "ccbar::Re{alpha_3^long}@GRvDV2021"], *this),
                    im_alpha_3_long(p[stringify(Process_::label) + "ccbar::Im{alpha_3^long}@GRvDV2021"], *this),
                    re_alpha_4_long(p[stringify(Process_::label) + "ccbar::Re{alpha_4^long}@GRvDV2021"], *this),
                    im_alpha_4_long(p[stringify(Process_::label) + "ccbar::Im{alpha_4^long}@GRvDV2021"], *this),
                    re_alpha_5_long(p[stringify(Process_::label) + "ccbar::Re{alpha_5^long}@GRvDV2021"], *this),
                    im_alpha_5_long(p[stringify(Process_::label) + "ccbar::Im{alpha_5^long}@GRvDV2021"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_" + opt_q.value()], *this),
                    m_Bsst(p["mass::B_s^*"], *this),

                    m_V(p["mass::" + _final_state()], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),

                    chiOPE(p["b->sccbar::chiOPE@GRvDV2021"], *this)
                {
                    this->uses(*form_factors);
                }

                ~GRvDV2021() = default;

                inline complex<double> phi(const complex<double> & q2, const std::array<unsigned, 5> & phi_parameters) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d    e
                    // 0(P->P) aka plus          5    3    2    2    2
                    // perp(P->V) = par(P->V)    5    1    3    0    2
                    // 0(P->V) aka long          5    1    2    2    2

                    const double m_V2  = power_of<2>(m_V);
                    const double m_Bsst2 = power_of<2>(m_Bsst);
                    const double m_B2  = power_of<2>(m_B),  m_B4 =  power_of<4>(m_B);
                    const double m_D02 = power_of<2>(m_D0), m_D04 = power_of<4>(m_D0);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * power_of<2>(m_D0), s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phi_parameters[0], b = phi_parameters[1], c = phi_parameters[2], d = phi_parameters[3], e = phi_parameters[4];

                    const complex<double> Nlambda = 4. * M_PI * pow(m_B2, 0.5 * (a - b + c + d - e) - 1.) * pow(2. * (4. * m_D02 - s_0) / 3. / chi, 0.5);
                    const complex<double> phi1 = -pow(2. * pow((4. * m_D02 - Q2) * (4. * m_D02 - s_0), 0.5) + 8. * m_D02 - Q2 - s_0, 0.5) /
                                                (2. * pow((4. * m_D02 - Q2) * (4. * m_D02 - s_0), 0.5) + 8. * m_D02 + Q2 * (z - 1.) - s_0*(z + 1.));
                    const complex<double> phi2 = pow(m_B4 * power_of<4>(z - 1.) - 2. * m_B2 * power_of<2>(z - 1.) * (-16 * m_D02 * z + m_V2 * power_of<2>(z - 1.) +
                                                s_0 * power_of<2>(z + 1.)) + power_of<2>(16 * m_D02 * z + m_V2 * power_of<2>(z - 1.) - s_0 * power_of<2>(z + 1.)), 0.5);
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0 * (z + 1.));
                    const complex<double> phi4 = pow(s_0 * power_of<2>(z + 1.) - 16. * z * m_D02, -0.5);
                    const complex<double> phi5 = pow(s_0 * power_of<2>(z + 1.) - 16. * z * m_D02 - m_Bsst2 * power_of<2>(-z + 1.), 0.5);

                    return Nlambda * pow(1. + z, 0.5) * pow(1. - z, a - b + c + d - e - 1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d) * pow(phi5, e);
                }

                inline complex<double> phi(const double & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    return phi(complex<double>(q2, 0.0), phi_parameters);
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const std::array<unsigned, 5> & phi_parameters, const std::array<complex<double>, 6> alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return eos::nff_utils::P<5u>(z_Jpsi, alpha) / phi(m_Jpsi2, phi_parameters) *
                            (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const std::array<unsigned, 5> & phi_parameters, const std::array<complex<double>, 6> alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return eos::nff_utils::P<5u>(z_psi2S, alpha) / phi(m_psi2S2, phi_parameters) *
                            (1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                };

                virtual complex<double> H_perp(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,              s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 5> phi_parameters = {5, 1, 3, 0, 2};

                    return eos::nff_utils::P<5u>(z, alpha_perp) / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_perp(const double & q2) const
                {
                    return H_perp(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_perp(const double & q2) const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return eos::nff_utils::P<5u>(z, alpha_perp);
                }

                virtual complex<double> H_para(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,              s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 5> phi_parameters = {5, 1, 3, 0, 2};

                    return eos::nff_utils::P<5u>(z, alpha_para) / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_para(const double & q2) const
                {
                    return H_para(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_para(const double & q2) const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return eos::nff_utils::P<5u>(z, alpha_para);
                }

                virtual complex<double> H_long(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,              s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 5> phi_parameters = {5, 1, 2, 2, 2};

                    return eos::nff_utils::P<5u>(z, alpha_long) / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_long(const double & q2) const
                {
                    return H_long(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_long(const double & q2) const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return eos::nff_utils::P<5u>(z, alpha_long);
                }

                virtual complex<double> H_perp_residue_jpsi() const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const std::array<unsigned, 5> phi_parameters = {5, 1, 3, 0, 2};

                    return H_residue_jpsi(phi_parameters, alpha_perp);
                }

                virtual complex<double> H_perp_residue_psi2s() const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const std::array<unsigned, 5> phi_parameters = {5, 1, 3, 0, 2};

                    return H_residue_psi2s(phi_parameters, alpha_perp);
                }

                virtual complex<double> H_para_residue_jpsi() const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const std::array<unsigned, 5> phi_parameters = {5, 1, 3, 0, 2};

                    return H_residue_jpsi(phi_parameters, alpha_para);
                }

                virtual complex<double> H_para_residue_psi2s() const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const std::array<unsigned, 5> phi_parameters = {5, 1, 3, 0, 2};

                    return H_residue_psi2s(phi_parameters, alpha_para);
                }

                   virtual complex<double> H_long_residue_jpsi() const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const std::array<unsigned, 5> phi_parameters = {5, 1, 2, 2, 2};

                    return H_residue_jpsi(phi_parameters, alpha_long);
                }

                virtual complex<double> H_long_residue_psi2s() const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const std::array<unsigned, 5> phi_parameters = {5, 1, 2, 2, 2};

                    return H_residue_psi2s(phi_parameters, alpha_long);
                }

                virtual complex<double> ratio_perp(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_perp = pow(2.0 * lambda, 0.5) / (m_B + m_V) / m_B() * form_factors->v(q2);

                    return H_perp(q2) / F_perp;
                }

                virtual complex<double> ratio_perp(const double & q2) const
                {
                    return ratio_perp(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_perp(const complex<double> & q2) const
                {
                    return (m_B + m_V) / m_B * form_factors->t_1(q2) / form_factors->v(q2);
                }

                virtual complex<double> ratio_para(const complex<double> & q2) const
                {
                    const complex<double> F_para = sqrt(2) * (m_B + m_V) / m_B * form_factors->a_1(q2);

                    return H_para(q2) / F_para;
                }

                virtual complex<double> ratio_para(const double & q2) const
                {
                    return ratio_para(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_para(const complex<double> & q2) const
                {
                    return (m_B - m_V) / m_B * form_factors->t_2(q2) / form_factors->a_1(q2);
                }

                virtual complex<double> ratio_long(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_long = ((m_B2 - m_V2 - q2) * power_of<2>(m_B + m_V) * form_factors->a_1(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * (m_B + m_V));

                    return H_long(q2) / F_long;
                }

                virtual complex<double> ratio_long(const double & q2) const
                {
                    return ratio_long(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_long(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_long = ((m_B2 - m_V2 - q2) * power_of<2>(m_B + m_V) * form_factors->a_1(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * (m_B + m_V));
                    const complex<double> F_T_long = q2 * ((m_B2 + 3 * m_V2 - q2) * (m_B2 - m_V2) * form_factors->t_2(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * m_B * (m_B2 - m_V2));

                    return F_T_long / F_long;
                }

                virtual complex<double> normalized_moment_V1(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V2(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V23(const double &) const
                {
                    return 0.0;
                }


                static NonlocalFormFactorPtr<nff::PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToV>(new GRvDV2021<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const std::array<unsigned, 5> & phi_parameters_long = {5, 1, 2, 2, 2}; //long polarization
                    results.add({ real(this->phi(16.0, phi_parameters_long)), "Re{phi_long(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phi_parameters_long)), "Im{phi_long(q2 = 16.0)}" });

                    const std::array<unsigned, 5> & phi_parameters_perp = {5, 1, 3, 0, 2}; //perp or para polarization
                    results.add({ real(this->phi(16.0, phi_parameters_perp)), "Re{phi_perp(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phi_parameters_perp)), "Im{phi_perp(q2 = 16.0)}" });

                    return results;
                }
        };
    }

    NonlocalFormFactorPtr<nff::PToV>
    NonlocalFormFactor<nff::PToV>::make(const QualifiedName & name, const Parameters & p, const Options & o)
    {
        typedef QualifiedName KeyType;
        typedef std::function<NonlocalFormFactorPtr<nff::PToV> (const Parameters &, const Options &)> ValueType;
        std::map<KeyType, ValueType> entries
        {
            // trivial
            std::make_pair("B->K^*::naive",         &nff_p_to_v::Naive::make),
            // parametrizations
            std::make_pair("B->K^*::GvDV2020",      &nff_p_to_v::GvDV2020<nff::BToKstar>::make),
            std::make_pair("B->K^*::GRvDV2021",     &nff_p_to_v::GRvDV2021<nff::BToKstar>::make),
            std::make_pair("B_s->phi::GvDV2020",    &nff_p_to_v::GvDV2020<nff::BsToPhi>::make),
            std::make_pair("B_s->phi::GRvDV2021",   &nff_p_to_v::GRvDV2021<nff::BsToPhi>::make),
        };

        auto i = entries.find(name);
        if (entries.end() == i)
        {
            return NonlocalFormFactorPtr<nff::PToV>(nullptr);
        }

        return i->second(p, o);
    }

    template <typename Process_>
    struct Implementation<NonlocalFormFactorObservable<Process_, nff::PToV>>
    {
        NameOption opt_formfactor;
        NonlocalFormFactorPtr<nff::PToV> nff;
        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_formfactor(o, "nonlocal-formfactor", qnp::Name("GvDV2020")),
            nff(NonlocalFormFactor<nff::PToV>::make(QualifiedName(qnp::Prefix(Process_::label), opt_formfactor.value()), p, o))
        {
            u.uses(*nff);
        }

        ~Implementation() = default;
    };

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nff::PToV>::NonlocalFormFactorObservable(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<NonlocalFormFactorObservable<Process_, nff::PToV>>(new Implementation<NonlocalFormFactorObservable<Process_, nff::PToV>>(p, o, *this))
    {
    }

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nff::PToV>::~NonlocalFormFactorObservable() = default;

    template <typename Process_>
    const std::vector<OptionSpecification>
    Implementation<NonlocalFormFactorObservable<Process_, nff::PToV>>::options
    {
    };

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_H_perp(const double & q2) const
    {
        return real(this->_imp->nff->H_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_H_perp(const double & q2) const
    {
        return imag(this->_imp->nff->H_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::abs_H_perp(const double & q2) const
    {
        return abs(this->_imp->nff->H_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_Hhat_perp(const double & q2) const
    {
        return real(this->_imp->nff->Hhat_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_Hhat_perp(const double & q2) const
    {
        return imag(this->_imp->nff->Hhat_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::abs_Hhat_perp(const double & q2) const
    {
        return abs(this->_imp->nff->Hhat_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_H_para(const double & q2) const
    {
        return real(this->_imp->nff->H_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_H_para(const double & q2) const
    {
        return imag(this->_imp->nff->H_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::abs_H_para(const double & q2) const
    {
        return abs(this->_imp->nff->H_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_Hhat_para(const double & q2) const
    {
        return real(this->_imp->nff->Hhat_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_Hhat_para(const double & q2) const
    {
        return imag(this->_imp->nff->Hhat_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::abs_Hhat_para(const double & q2) const
    {
        return abs(this->_imp->nff->Hhat_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_H_long(const double & q2) const
    {
        return real(this->_imp->nff->H_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_H_long(const double & q2) const
    {
        return imag(this->_imp->nff->H_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::abs_H_long(const double & q2) const
    {
        return abs(this->_imp->nff->H_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_Hhat_long(const double & q2) const
    {
        return real(this->_imp->nff->Hhat_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_Hhat_long(const double & q2) const
    {
        return imag(this->_imp->nff->Hhat_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::abs_Hhat_long(const double & q2) const
    {
        return abs(this->_imp->nff->Hhat_long(q2));
    }


    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_ratio_perp(const double & q2) const
    {
        return real(this->_imp->nff->ratio_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_ratio_perp(const double & q2) const
    {
        return imag(this->_imp->nff->ratio_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::abs_ratio_perp(const double & q2) const
    {
        return abs(this->_imp->nff->ratio_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_ratio_para(const double & q2) const
    {
        return real(this->_imp->nff->ratio_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_ratio_para(const double & q2) const
    {
        return imag(this->_imp->nff->ratio_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::abs_ratio_para(const double & q2) const
    {
        return abs(this->_imp->nff->ratio_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_ratio_long(const double & q2) const
    {
        return real(this->_imp->nff->ratio_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_ratio_long(const double & q2) const
    {
        return imag(this->_imp->nff->ratio_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::abs_ratio_long(const double & q2) const
    {
        return abs(this->_imp->nff->ratio_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_ratio_perp_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->ratio_perp(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_ratio_perp_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->ratio_perp(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_ratio_para_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->ratio_para(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_ratio_para_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->ratio_para(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_ratio_long_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->ratio_long(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_ratio_long_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->ratio_long(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_F_ratio_perp_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->F_ratio_perp(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_F_ratio_perp_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->F_ratio_perp(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_F_ratio_para_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->F_ratio_para(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_F_ratio_para_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->F_ratio_para(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_F_ratio_long_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->F_ratio_long(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::im_F_ratio_long_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->F_ratio_long(complex<double>(re_q2, im_q2)));
    }


    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_normalized_moment_V1(const double & q2) const
    {
        return real(this->_imp->nff->normalized_moment_V1(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_normalized_moment_V2(const double & q2) const
    {
        return real(this->_imp->nff->normalized_moment_V2(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToV>::re_normalized_moment_V23(const double & q2) const
    {
        return real(this->_imp->nff->normalized_moment_V23(q2));
    }

    template <typename Process_>
    const std::set<ReferenceName>
    NonlocalFormFactorObservable<Process_, nff::PToV>::references
    {
        "GvDV:2020A"_rn
    };

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    NonlocalFormFactorObservable<Process_, nff::PToV>::begin_options()
    {
        return Implementation<NonlocalFormFactorObservable<Process_, nff::PToV>>::options.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    NonlocalFormFactorObservable<Process_, nff::PToV>::end_options()
    {
        return Implementation<NonlocalFormFactorObservable<Process_, nff::PToV>>::options.cend();
    }

    template class NonlocalFormFactorObservable<nff::BToKstar, nff::PToV>;

    template class NonlocalFormFactorObservable<nff::BsToPhi, nff::PToV>;
}
