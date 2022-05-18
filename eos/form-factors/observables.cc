/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2019, 2020 Danny van Dyk
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

#include <eos/observable-impl.hh>
#include <eos/form-factors/form-factor-adapter.hh>
#include <eos/form-factors/analytic-b-to-pi-pi.hh>
#include <eos/form-factors/analytic-b-to-v-lcsr.hh>
#include <eos/form-factors/analytic-b-to-p-lcsr.hh>
#include <eos/form-factors/b-lcdas.hh>
#include <eos/form-factors/baryonic-impl.hh>
#include <eos/form-factors/mesonic-impl.hh>
#include <eos/form-factors/mesonic-hqet.hh>
#include <eos/form-factors/observables.hh>
#include <eos/form-factors/parametric-bmrvd2022.hh>
#include <eos/form-factors/unitarity-bounds.hh>
#include <eos/form-factors/zero-recoil-sum-rule.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    /* form factors as observables */
    template <typename Transition_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_form_factor_adapter(const char * name,
            const char * latex,
            double (FormFactors<Transition_>::* _function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        QualifiedName qn(name);
        qnp::Prefix pp = qn.prefix_part();
        std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> function(_function);

        auto result = std::make_pair(qn, std::make_shared<FormFactorAdapterEntry<Transition_, Args_ ...>>(qn, latex, Unit::None(), pp, function, kinematics_names));

        impl::observable_entries.insert(result);

        return result;
    }

    template <typename Transition_, typename Tuple_, typename ... Args_>
    std::pair<QualifiedName, ObservableEntryPtr> make_form_factor_adapter(const char * name,
            double (FormFactors<Transition_>::* _function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        QualifiedName qn(name);
        qnp::Prefix pp = qn.prefix_part();
        std::function<double (const FormFactors<Transition_> *, const Args_ & ...)> function(_function);

        auto result = std::make_pair(qn, std::make_shared<FormFactorAdapterEntry<Transition_, Args_ ...>>(qn, "", Unit::None(), pp, function, kinematics_names));

        impl::observable_entries.insert(result);

        return result;
    }

    // B -> P(seudoscalar)
    // {{{

    // B -> pi
    // {{{
    ObservableGroup
    make_b_to_pi_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \pi$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to \pi$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->pi::f_+(q2)", R"(f_+^{B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_+'(q2)", R"(f_+^{',B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_p_d1, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_+''(q2)", R"(f_+^{'',B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_p_d2, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_T(q2)", R"(f_T^{B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_0(q2)", R"(f_0^{B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->pi::f_-(q2)", R"(f_-^{B\to\pi}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_expression_observable("B->pi::f_0(q2)/f_+(q2)", R"(f_0(q^2)/f_+(q^2))",
                        Unit::None(),
                        R"( <<B->pi::f_0(q2)>> / <<B->pi::f_+(q2)>> )"),

                // auxiliary variables, e.g. for determining the B-LCSR threshold parameters
                make_observable("B->pi::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BToPi>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B->pi::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BToPi>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B->pi::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BToPi>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                // auxiliary variables, e.g. for determining the pi-LCSR/SVZ threshold parameters
                make_observable("B->pi::M_B(f_+,LCSR)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPiDKMMO2008::MBp_lcsr,
                        std::make_tuple("q2")),

                make_observable("B->pi::M_B(f_0,LCSR)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPiDKMMO2008::MB0_lcsr,
                        std::make_tuple("q2")),

                make_observable("B->pi::M_B(f_T,LCSR)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPiDKMMO2008::MBT_lcsr,
                        std::make_tuple("q2")),

                make_observable("B->pi::M_B(SVZ)@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPiDKMMO2008::MB_svz),

                make_observable("B->pi::f_B@DKMMO2008",
                        Unit::GeV(),
                        &AnalyticFormFactorBToPiDKMMO2008::decay_constant),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> K
    // {{{
    ObservableGroup
    make_b_to_k_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to K$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to K$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->K::f_+(q2)", R"(f_+^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->K::f_0(q2)", R"(f_0^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->K::f_T(q2)", R"(f_T^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B->K::f_-(q2)", R"(f_-^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_observable("B->K::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BToK>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B->K::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BToK>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B->K::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BToK>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                make_form_factor_adapter("B->K::F_plus(q2)", R"(F_+^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->K::F_plus_T(q2)", R"(F_T^{B\to K}(q^2))",
                        &FormFactors<PToP>::f_plus_T, std::make_tuple("q2")),

                make_expression_observable("B->K::F_T(q2)/F_plus(q2)", R"(F_T(q^2)/F_+(q^2))",
                        Unit::None(),
                        R"( <<B->K::f_T(q2)>> / <<B->K::F_plus(q2)>> )"),

                make_expression_observable("B->K::F_plus_T(q2)/F_plus(q2)", R"(F_{+,T}(q^2)/F_+(q^2))",
                        Unit::None(),
                        R"( <<B->K::F_plus_T(q2)>> / <<B->K::F_plus(q2)>> )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> D
    // {{{
    ObservableGroup
    make_b_to_d_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \bar{D}$ transitions and related pseudo-observables)",
            R"(Pseudo observables representing the full basis of $B\to\bar{D}$ form factors. )"
            R"(For most pseudo-observables, the specific parametrization can be chosen via the "form-factors" option.)",
            {
                // B -> D Form Factors
                make_form_factor_adapter("B->D::f_+(q2)", R"(f_+^{B\to \bar{D}}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B->D::f_0(q2)", R"(f_0^{B\to \bar{D}}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->D::f_T(q2)", R"(f_T^{B\to \bar{D}}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B->D::f_-(q2)", R"(f_-^{B\to \bar{D}}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_observable("B->D::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BToD>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B->D::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BToD>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B->D::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BToD>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                make_observable("B->D::a_0[S_1]@HQE", R"(a_0^{S_1})",
                        Unit::None(),
                        &BGLCoefficients::S1_a0),

                make_observable("B->D::a_1[S_1]@HQE", R"(a_1^{S_1})",
                        Unit::None(),
                        &BGLCoefficients::S1_a1),

                make_observable("B->D::a_2[S_1]@HQE", R"(a_2^{S_1})",
                        Unit::None(),
                        &BGLCoefficients::S1_a2),

                make_expression_observable("B->D^*::a_1/a_0[S_1]@HQE", R"(a_1^{S_1}/a_0^{S_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_1[S_1]@HQE>>
                        /
                        <<B->D::a_0[S_1]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[S_1]@HQE", R"(a_2^{S_1}/a_0^{S_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_2[S_1]@HQE>>
                        /
                        <<B->D::a_0[S_1]@HQE>>
                        )"),

                make_observable("B->D::a_0[V_1]@HQE", R"(a_0^{V_1})",
                        Unit::None(),
                        &BGLCoefficients::V1_a0),

                make_observable("B->D::a_1[V_1]@HQE", R"(a_1^{V_1})",
                        Unit::None(),
                        &BGLCoefficients::V1_a1),

                make_observable("B->D::a_2[V_1]@HQE", R"(a_2^{V_1})",
                        Unit::None(),
                        &BGLCoefficients::V1_a2),

                make_expression_observable("B->D^*::a_1/a_0[V_1]@HQE", R"(a_1^{V_1}/a_0^{V_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_1[V_1]@HQE>>
                        /
                        <<B->D::a_0[V_1]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[V_1]@HQE", R"(a_2^{V_1}/a_0^{V_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_2[V_1]@HQE>>
                        /
                        <<B->D::a_0[V_1]@HQE>>
                        )"),

                make_expression_observable("B->D::f_T(q2)/f_+(q2)", R"(f_T(q^2)/f_+(q^2))",
                        Unit::None(),
                        R"( <<B->D::f_T(q2)>> / <<B->D::f_+(q2)>> )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B_s -> P(seudoscalar)
    // {{{

    // B_s -> K
    // {{{
    ObservableGroup
    make_bs_to_k_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \bar{K}$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \bar{K}$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B_s->K::f_+(q2)", R"(f_+^{B_s\to \bar{K}}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K::f_0(q2)", R"(f_0^{B_s\to \bar{K}}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K::f_T(q2)", R"(f_T^{B_s\to \bar{K}}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K::f_-(q2)", R"(f_-^{B_s\to \bar{K}}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_observable("B_s->K::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToK>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B_s->K::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToK>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B_s->K::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToK>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> D_s
    // {{{
    ObservableGroup
    make_bs_to_ds_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \bar{D_s}$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to\bar{D}_s$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                // B -> D Form Factors
                make_form_factor_adapter("B_s->D_s::f_+(q2)", R"(f_+^{B_s\to \bar{D}_s}(q^2))",
                        &FormFactors<PToP>::f_p, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s::f_0(q2)", R"(f_0^{B_s\to \bar{D}_s}(q^2))",
                        &FormFactors<PToP>::f_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s::f_T(q2)", R"(f_T^{B_s\to \bar{D}_s}(q^2))",
                        &FormFactors<PToP>::f_t, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s::f_-(q2)", R"(f_-^{B_s\to \bar{D}_s}(q^2))",
                        &FormFactors<PToP>::f_m, std::make_tuple("q2")),

                make_observable("B_s->D_s::f_+[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToDs>::normalized_moment_1_f_p,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s::f_0[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToDs>::normalized_moment_1_f_pm,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s::f_T[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToPLCSR<lcsr::BsToDs>::normalized_moment_1_f_t,
                        std::make_tuple("q2")),

                make_expression_observable("B(_s)->D(_s)::f_0(q2_num)/f_0(q2_denom)", R"(f_0(q^2_\mathrm{num})/f_+(q^2_\mathrm{denom}))",
                        Unit::None(),
                        R"(
                        <<B_s->D_s::f_0(q2)>>[q2=>q2_num]
                        /
                        <<B->D::f_0(q2)>>[q2=>q2_denom]
                        )"),

                make_expression_observable("B_s->D_s::f_T(q2)/f_+(q2)", R"(f_T(q^2)/f_+(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s::f_T(q2)>> / <<B_s->D_s::f_+(q2)>> )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B -> V(ector)
    // {{{

    // B -> omega
    // {{{
    ObservableGroup
    make_b_to_omega_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \omega$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to \omega$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->omega::V(q2)", R"(V^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::A_0(q2)", R"(A_0^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::A_1(q2)", R"(A_1^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::A_2(q2)", R"(A_2^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::A_12(q2)", R"(A_{12}^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::T_1(q2)", R"(T_1^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::T_2(q2)", R"(T_2^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::T_3(q2)", R"(T_3^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->omega::T_23(q2)", R"(T_{23}^{B\to \omega}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> rho
    // {{{
    ObservableGroup
    make_b_to_rho_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \rho$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to \rho$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->rho::V(q2)", R"(V^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_0(q2)", R"(A_0^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_1(q2)", R"(A_1^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_2(q2)", R"(A_2^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::A_12(q2)", R"(A_{12}^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::T_1(q2)", R"(T_1^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::T_2(q2)", R"(T_2^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::T_3(q2)", R"(T_3^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->rho::T_23(q2)", R"(T_{23}^{B\to \rho}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_expression_observable("B->rho::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->rho::V(q2)>> / <<B->rho::A_1(q2)>> )"),

                make_expression_observable("B->rho::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->rho::A_2(q2)>> / <<B->rho::A_1(q2)>> )"),

                make_expression_observable("B->rho::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->rho::A_12(q2)>> / <<B->rho::A_1(q2)>> )"),

                make_expression_observable("B->rho::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B->rho::T_23(q2)>> / <<B->rho::T_2(q2)>> )"),

                make_observable("B->rho::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B->rho::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B->rho::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B->rho::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B->rho::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B->rho::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B->rho::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> K^*
    // {{{
    ObservableGroup
    make_b_to_kstar_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to K^*$ transitions)",
            R"(Pseudo observables representing the full basis of $B\to K^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->K^*::V(q2)", R"(V^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_0(q2)", R"(A_0^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_1(q2)", R"(A_1^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_2(q2)", R"(A_2^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::A_12(q2)", R"(A_{12}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_1(q2)", R"(T_1^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_2(q2)", R"(T_2^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_3(q2)", R"(T_3^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::T_23(q2)", R"(T_{23}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_perp(q2)", R"(\mathcal{F}_\perp^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_perp, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_para(q2)",  R"(\mathcal{F}_\parallel^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_para, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_long(q2)",  R"(\mathcal{F}_0^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_long, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_perp_T(q2)", R"(\mathcal{F}_{\perp,T}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_perp_T, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_para_T(q2)", R"(\mathcal{F}_{\parallel,T}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_para_T, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_long_T(q2)", R"(\mathcal{F}_{0,T}^{B\to K^*}(q^2))",
                        &FormFactors<PToV>::f_long_T, std::make_tuple("q2")),

                make_form_factor_adapter("B->K^*::F_long_T(q2)_Normalized",
                        &FormFactors<PToV>::f_long_T_Normalized, std::make_tuple("q2")),

                make_expression_observable("B->K^*::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::V(q2)>> / <<B->K^*::A_1(q2)>> )"),

                make_expression_observable("B->K^*::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::A_2(q2)>> / <<B->K^*::A_1(q2)>> )"),

                make_expression_observable("B->K^*::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::A_12(q2)>> / <<B->K^*::A_1(q2)>> )"),

                make_expression_observable("B->K^*::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::T_23(q2)>> / <<B->K^*::T_2(q2)>> )"),

                make_observable("B->K^*::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B->K^*::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B->K^*::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B->K^*::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B->K^*::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B->K^*::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B->K^*::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),

                make_expression_observable("B->K^*::F_perp_T(q2)/F_perp(q2)", R"(\mathcal{F}_{\perp,T}(q^2)/\mathcal{F}_\perp(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::F_perp_T(q2)>> / <<B->K^*::F_perp(q2)>> )"),

                make_expression_observable("B->K^*::F_para_T(q2)/F_para(q2)", R"(\mathcal{F}_{\parallel,T}(q^2)/\mathcal{F}_\parallel(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::F_para_T(q2)>> / <<B->K^*::F_para(q2)>> )"),

                make_expression_observable("B->K^*::F_long_T(q2)/F_long(q2)", R"(\mathcal{F}_{0,T}(q^2)/\mathcal{F}_0(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::F_long_T(q2)>> / <<B->K^*::F_long(q2)>> )"),

                make_expression_observable("B->K^*::F_long_T(q2)/F_long(q2)_Normalized", R"(\mathcal{N} \mathcal{F}_{0,T}(q^2)/\mathcal{F}_0(q^2))",
                        Unit::None(),
                        R"( <<B->K^*::F_long_T(q2)_Normalized>> / <<B->K^*::F_long(q2)>> )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> D^*
    // {{{
    ObservableGroup
    make_b_to_dstar_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B\to \bar{D}^*$ transitions and related pseudo-observables)",
            R"(Pseudo observables representing the full basis of $B\to \bar{D}^*$ form factors. )"
            R"(For most pseudo-observables, the specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->D^*::V(q2)", R"(V^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_0(q2)", R"(V^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_1(q2)", R"(A_1^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_2(q2)", R"(A_2^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::A_12(q2)", R"(A_{12}^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_1(q2)", R"(T_1^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_2(q2)", R"(T_2^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_3(q2)", R"(T_3^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B->D^*::T_23(q2)", R"(T_{23}^{B\to D^*}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_expression_observable("B->D^*::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->D^*::V(q2)>> / <<B->D^*::A_1(q2)>> )"),

                make_expression_observable("B->D^*::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->D^*::A_2(q2)>> / <<B->D^*::A_1(q2)>> )"),

                make_expression_observable("B->D^*::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B->D^*::A_12(q2)>> / <<B->D^*::A_1(q2)>> )"),

                make_expression_observable("B->D^*::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B->D^*::T_23(q2)>> / <<B->D^*::T_2(q2)>> )"),

                make_observable("B->D^*::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B->D^*::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B->D^*::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                         &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B->D^*::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B->D^*::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B->D^*::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B->D^*::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),

                make_observable("B->D^*::a_0[A_1]@HQE", R"(a_0^{A_1})",
                        Unit::None(),
                        &BGLCoefficients::A1_a0),

                make_observable("B->D^*::a_1[A_1]@HQE", R"(a_1^{A_1})",
                        Unit::None(),
                        &BGLCoefficients::A1_a1),

                make_observable("B->D^*::a_2[A_1]@HQE", R"(a_2^{A_1})",
                        Unit::None(),
                        &BGLCoefficients::A1_a2),

                make_expression_observable("B->D^*::a_1/a_0[A_1]@HQE", R"(a_1^{A_1}/a_0^{A_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_1[A_1]@HQE>>
                        /
                        <<B->D::a_0[A_1]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[A_1]@HQE", R"(a_2^{A_1}/a_0^{A_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_2[A_1]@HQE>>
                        /
                        <<B->D::a_0[A_1]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[A_5]@HQE", R"(a_0^{A_5})",
                        Unit::None(),
                        &BGLCoefficients::A5_a0),

                make_observable("B->D^*::a_1[A_5]@HQE", R"(a_1^{A_5})",
                        Unit::None(),
                        &BGLCoefficients::A5_a1),

                make_observable("B->D^*::a_2[A_5]@HQE", R"(a_2^{A_5})",
                        Unit::None(),
                        &BGLCoefficients::A5_a2),

                make_expression_observable("B->D^*::a_1/a_0[A_5]@HQE", R"(a_1^{A_5}/a_0^{A_5})",
                        Unit::None(),
                        R"(
                        <<B->D::a_1[A_5]@HQE>>
                        /
                        <<B->D::a_0[A_5]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[A_5]@HQE", R"(a_2^{A_5}/a_0^{A_5})",
                        Unit::None(),
                        R"(
                        <<B->D::a_2[A_5]@HQE>>
                        /
                        <<B->D::a_0[A_5]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[P_1]@HQE", R"(a_0^{P_1})",
                        Unit::None(),
                        &BGLCoefficients::P1_a0),

                make_observable("B->D^*::a_1[P_1]@HQE", R"(a_1^{P_1})",
                        Unit::None(),
                        &BGLCoefficients::P1_a1),

                make_observable("B->D^*::a_2[P_1]@HQE", R"(a_2^{P_1})",
                        Unit::None(),
                        &BGLCoefficients::P1_a2),

                make_expression_observable("B->D^*::a_1/a_0[P_1]@HQE", R"(a_1^{P_1}/a_0^{P_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_1[P_1]@HQE>>
                        /
                        <<B->D::a_0[P_1]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[P_1]@HQE", R"(a_2^{P_1}/a_0^{P_1})",
                        Unit::None(),
                        R"(
                        <<B->D::a_2[P_1]@HQE>>
                        /
                        <<B->D::a_0[P_1]@HQE>>
                        )"),

                make_observable("B->D^*::a_0[V_4]@HQE", R"(a_0^{V_4})",
                        Unit::None(),
                        &BGLCoefficients::V4_a0),

                make_observable("B->D^*::a_1[V_4]@HQE", R"(a_1^{V_4})",
                        Unit::None(),
                        &BGLCoefficients::V4_a1),

                make_observable("B->D^*::a_2[V_4]@HQE", R"(a_2^{V_4})",
                        Unit::None(),
                        &BGLCoefficients::V4_a2),

                make_expression_observable("B->D^*::a_1/a_0[V_4]@HQE", R"(a_1^{V_4}/a_0^{V_4})",
                        Unit::None(),
                        R"(
                        <<B->D::a_1[V_4]@HQE>>
                        /
                        <<B->D::a_0[V_4]@HQE>>
                        )"),

                make_expression_observable("B->D^*::a_2/a_0[V_4]@HQE", R"(a_2^{V_4}/a_0^{V_4})",
                        Unit::None(),
                        R"(
                        <<B->D::a_2[V_4]@HQE>>
                        /
                        <<B->D::a_0[V_4]@HQE>>
                        )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B_s -> V(ector)
    // {{{

    // B_s -> K^*
    // {{{
    ObservableGroup
    make_bs_to_kstar_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \bar{K}^*$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \bar{K}^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B_s->K^*::V(q2)", R"(V^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_0(q2)", R"(A_0^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_1(q2)", R"(A_1^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_2(q2)", R"(A_2^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::A_12(q2)", R"(A_{12}^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::T_1(q2)", R"(T_1^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::T_2(q2)", R"(T_2^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::T_3(q2)", R"(T_3^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->K^*::T_23(q2)", R"(T_{23}^{B_s\to \bar{K}^*}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_expression_observable("B_s->K^*::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->K^*::V(q2)>> / <<B_s->K^*::A_1(q2)>> )"),

                make_expression_observable("B_s->K^*::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->K^*::A_2(q2)>> / <<B_s->K^*::A_1(q2)>> )"),

                make_expression_observable("B_s->K^*::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->K^*::A_12(q2)>> / <<B_s->K^*::A_1(q2)>> )"),

                make_expression_observable("B_s->K^*::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B_s->K^*::T_23(q2)>> / <<B_s->K^*::T_2(q2)>> )"),

                make_observable("B_s->K^*::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B_s->K^*::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> phi
    // {{{
    ObservableGroup
    make_bs_to_phi_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \phi$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \phi$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B_s->phi::V(q2)", R"(V^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_0(q2)", R"(A_0^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_1(q2)", R"(A_1^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_2(q2)", R"(A_2^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::A_12(q2)", R"(A_{12}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::T_1(q2)", R"(T_1^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::T_2(q2)", R"(T_2^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::T_3(q2)", R"(T_3^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::T_23(q2)", R"(T_{23}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_perp(q2)", R"(\mathcal{F}_\perp^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_perp, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_para(q2)", R"(\mathcal{F}_\parallel^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_para, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_long(q2)", R"(\mathcal{F}_0^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_long, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_perp_T(q2)", R"(\mathcal{F}_{\perp,T}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_perp_T, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_para_T(q2)", R"(\mathcal{F}_{\parallel,T}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_para_T, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_long_T(q2)", R"(\mathcal{F}_{0,T}^{B_s\to \phi}(q^2))",
                        &FormFactors<PToV>::f_long_T, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->phi::F_long_T(q2)_Normalized",
                        &FormFactors<PToV>::f_long_T_Normalized, std::make_tuple("q2")),

                make_expression_observable("B_s->phi::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::V(q2)>> / <<B_s->phi::A_1(q2)>> )"),

                make_expression_observable("B_s->phi::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::A_2(q2)>> / <<B_s->phi::A_1(q2)>> )"),

                make_expression_observable("B_s->phi::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::A_12(q2)>> / <<B_s->phi::A_1(q2)>> )"),

                make_expression_observable("B_s->phi::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::T_23(q2)>> / <<B_s->phi::T_2(q2)>> )"),

                make_observable("B_s->phi::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B_s->phi::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),

                make_expression_observable("B_s->phi::F_perp_T(q2)/F_perp(q2)", R"(\mathcal{F}_{\perp,T}(q^2)/\mathcal{F}_\perp(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::F_perp_T(q2)>> / <<B_s->phi::F_perp(q2)>> )"),

                make_expression_observable("B_s->phi::F_para_T(q2)/F_para(q2)", R"(\mathcal{F}_{\parallel,T}(q^2)/\mathcal{F}_\parallel(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::F_para_T(q2)>> / <<B_s->phi::F_para(q2)>> )"),

                make_expression_observable("B_s->phi::F_long_T(q2)/F_long(q2)", R"(\mathcal{F}_{0,T}(q^2)/\mathcal{F}_0(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::F_long_T(q2)>> / <<B_s->phi::F_long(q2)>> )"),

                make_expression_observable("B_s->phi::F_long_T(q2)/F_long(q2)_Normalized", R"(\mathcal{N} \mathcal{F}_{0,T}(q^2)/\mathcal{F}_0(q^2))",
                        Unit::None(),
                        R"( <<B_s->phi::F_long_T(q2)_Normalized>> / <<B_s->phi::F_long(q2)>> )"),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_s -> D_s^*
    // {{{
    ObservableGroup
    make_bs_to_dsstar_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B_s\to \bar{D}_s^*$ transitions)",
            R"(Pseudo observables representing the full basis of $B_s\to \bar{D}_s^*$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B_s->D_s^*::V(q2)", R"(V^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::v, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_0(q2)", R"(A_0^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::a_0, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_1(q2)", R"(A_1^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::a_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_2(q2)", R"(A_2^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::a_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::A_12(q2)", R"(A_{12}^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::a_12, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_1(q2)", R"(T_1^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::t_1, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_2(q2)", R"(T_2^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::t_2, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_3(q2)", R"(T_3^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::t_3, std::make_tuple("q2")),

                make_form_factor_adapter("B_s->D_s^*::T_23(q2)", R"(T_{23}^{B_s\to D_s^*}(q^2))",
                        &FormFactors<PToV>::t_23, std::make_tuple("q2")),

                make_expression_observable("B_s->D_s^*::V(q2)/A_1(q2)", R"(V(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s^*::V(q2)>> / <<B_s->D_s^*::A_1(q2)>> )"),

                make_expression_observable("B_s->D_s^*::A_2(q2)/A_1(q2)", R"(A_2(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s^*::A_2(q2)>> / <<B_s->D_s^*::A_1(q2)>> )"),

                make_expression_observable("B_s->D_s^*::A_12(q2)/A_1(q2)", R"(A_{12}(q^2)/A_1(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s^*::A_12(q2)>> / <<B_s->D_s^*::A_1(q2)>> )"),

                make_expression_observable("B_s->D_s^*::T_23(q2)/T_2(q2)", R"(T_{23}(q^2)/T_2(q^2))",
                        Unit::None(),
                        R"( <<B_s->D_s^*::T_23(q2)>> / <<B_s->D_s^*::T_2(q2)>> )"),

                make_observable("B_s->D_s^*::A_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_a_1,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::A_2[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_a_2,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::A_30[s^1/s^0](q2)",
                        Unit::GeV2(),
                         &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_a_30,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::V[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_v,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::T_1[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_t_1,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::T_23A[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_t_23A,
                        std::make_tuple("q2")),

                make_observable("B_s->D_s^*::T_23B[s^1/s^0](q2)",
                        Unit::GeV2(),
                        &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::normalized_moment_1_t_23B,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // B -> P P
    // {{{

    // B -> pi pi
    // {{{
    ObservableGroup
    make_b_to_pi_pi_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $B \to \pi \pi$ transitions)",
            R"(Pseudo observables representing the $B \to \pi \pi$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("B->pipi::Im{F_perp}(q2,k2,z)", R"(\mathrm{Im}\,F_\perp^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_perp, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{F_para}(q2,k2,z)", R"(\mathrm{Im}\,F_\parallel^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_para, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{F_long}(q2,k2,z)", R"(\mathrm{Im}\,F_0^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_long, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{F_time}(q2,k2,z)", R"(\mathrm{Im}\,F_t^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::im_f_time, std::make_tuple("q2", "k2", "z")),

                make_form_factor_adapter("B->pipi::Im{Res{F_perp}}(q2,k2)", R"(\mathrm{Res}\,\mathrm{Im}\,F_\perp^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::f_perp_im_res_qhat2, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->pipi::Im{Res{F_para}}(q2,k2)", R"(\mathrm{Res}\,\mathrm{Im}\,F_\parallel^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::f_para_im_res_qhat2, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->pipi::Im{Res{F_long}}(q2,k2)", R"(\mathrm{Res}\,\mathrm{Im}\,F_0^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::f_long_im_res_qhat2, std::make_tuple("q2", "k2")),

                make_form_factor_adapter("B->pipi::Im{Res{F_time}}(q2,k2)", R"(\mathrm{Res}\,\mathrm{Im}\,F_t^{B\to \pi\pi}(q^2,k^2,z))",
                        &FormFactors<PToPP>::f_time_im_res_qhat2, std::make_tuple("q2", "k2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // 1/2^+ -> 1/2^+
    // {{{

    // Lambda_b -> Lambda
    // {{{
    ObservableGroup
    make_lambdab_to_lambda_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $\Lambda_b \to \Lambda$ transitions)",
            R"(Pseudo observables representing the full basis of $\Lambda_b \to \Lambda$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("Lambda_b->Lambda::f_time^V(q2)", R"(f_t^{V,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_long^V(q2)", R"(f_0^{V,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_perp^V(q2)", R"(f_\perp^{V,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_time^A(q2)", R"(f_t^{A,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_long^A(q2)", R"(f_0^{A,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_perp^A(q2)", R"(f_\perp^{A,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_long^T(q2)", R"(f_0^{T,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_perp^T(q2)", R"(f_\perp^{T,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_long^T5(q2)", R"(f_0^{T5,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t5, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda::f_perp^T5(q2)", R"(f_\perp^{T5,\Lambda_b\to\Lambda}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t5, std::make_tuple("q2")),

                make_observable("Lambda_b->Lambda::Bound[0^+]@BMRvD2022", R"(\textrm{Bound[$0^+$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_0p),

                make_observable("Lambda_b->Lambda::Bound[1^-]@BMRvD2022", R"(\textrm{Bound[$1^-$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_1m),

                make_observable("Lambda_b->Lambda::Bound[0^-]@BMRvD2022", R"(\textrm{Bound[$0^-$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_0m),

                make_observable("Lambda_b->Lambda::Bound[1^+]@BMRvD2022", R"(\textrm{Bound[$1^+$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_1p),

                make_observable("Lambda_b->Lambda::Bound[T]@BMRvD2022", R"(\textrm{Bound[$T$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_T),

                make_observable("Lambda_b->Lambda::Bound[T5]@BMRvD2022", R"(\textrm{Bound[$T5$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_T5),

                make_observable("Lambda_b->Lambda::Prior[0^+]@BMRvD2022", R"(\textrm{Prior[$0^+$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_0p_prior),

                make_observable("Lambda_b->Lambda::Prior[1^-]@BMRvD2022", R"(\textrm{Prior[$1^-$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_1m_prior),

                make_observable("Lambda_b->Lambda::Prior[0^-]@BMRvD2022", R"(\textrm{Prior[$0^-$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_0m_prior),

                make_observable("Lambda_b->Lambda::Prior[1^+]@BMRvD2022", R"(\textrm{Prior[$1^+$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_1p_prior),

                make_observable("Lambda_b->Lambda::Prior[T]@BMRvD2022", R"(\textrm{Prior[$T$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_T_prior),

                make_observable("Lambda_b->Lambda::Prior[T5]@BMRvD2022", R"(\textrm{Prior[$T5$]})", Unit::None(),
                        &BMRvD2022FormFactors<LambdaBToLambda>::bound_T5_prior),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Lambda_b -> Lambda_c
    // {{{
    ObservableGroup
    make_lambdab_to_lambdac_form_factors_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Form factors for $\Lambda_b \to \Lambda_c$ transitions)",
            R"(Pseudo observables representing the full basis of $\Lambda_b \to \Lambda_c$ form factors. )"
            R"(The specific parametrization can be chosen via the "form-factors" option.)",
            {
                make_form_factor_adapter("Lambda_b->Lambda_c::f_time^V(q2)", R"(f_t^{V,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_long^V(q2)", R"(f_0^{V,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_perp^V(q2)", R"(f_\perp^{V,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_v, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_time^A(q2)", R"(f_t^{A,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_time_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_long^A(q2)", R"(f_0^{A,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_perp^A(q2)", R"(f_\perp^{A,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_a, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_long^T(q2)", R"(f_0^{T,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_perp^T(q2)", R"(f_\perp^{T,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_long^T5(q2)", R"(f_0^{T5,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_t5, std::make_tuple("q2")),

                make_form_factor_adapter("Lambda_b->Lambda_c::f_perp^T5(q2)", R"(f_\perp^{T5,\Lambda_b\to\Lambda_c}(q^2))",
                        &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_t5, std::make_tuple("q2")),

                // Zero-Recoil Sum Rule for the Lambda_b -> Lambda_c Form Factors
                make_observable("Lambda_b->Lambda_c::F(1)",
                        Unit::None(),
                        &ZeroRecoilSumRule<LambdaBToC>::vector_current),

                make_observable("Lambda_b->Lambda_c::G(1)",
                        Unit::None(),
                        &ZeroRecoilSumRule<LambdaBToC>::axialvector_current),

                make_observable("Lambda_b->Lambda_c::F_inel(1)",
                        Unit::None(),
                        &ZeroRecoilSumRule<LambdaBToC>::vector_current_inel),

                make_observable("Lambda_b->Lambda_c::G_inel(1)",
                        Unit::None(),
                        &ZeroRecoilSumRule<LambdaBToC>::axialvector_current_inel),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // }}}

    // unitarity bounds
    // {{{

    ObservableGroup
    make_unitarity_bounds_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Unitarity Bounds)",
            R"(Pseudo observables arising in the various unitarity bounds for $b\to c$ semileptonic form factors.)",
            {
                make_observable("b->c::Bound[0^+]@CLN", R"(B^{b\to c}_{0^+})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_0p),

                make_observable("b->c::Bound[0^-]@CLN", R"(B^{b\to c}_{0^-})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_0m),

                make_observable("b->c::Bound[1^+]@CLN", R"(B^{b\to c}_{1^+})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_1p),

                make_observable("b->c::Bound[1^-]@CLN", R"(B^{b\to c}_{1^-})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_1m),

                make_observable("b->c::Prior[0^+]@CLN", R"(B^{b\to c}_{0^+})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_0p_prior),

                make_observable("b->c::Prior[0^-]@CLN", R"(B^{b\to c}_{0^-})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_0m_prior),

                make_observable("b->c::Prior[1^+]@CLN", R"(B^{b\to c}_{1^+})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_1p_prior),

                make_observable("b->c::Prior[1^-]@CLN", R"(B^{b\to c}_{1^-})",
                        Unit::None(),
                        &HQETUnitarityBounds::bound_1m_prior),

                make_observable("b->c::Bound[0^+]@OPE", R"(B^{b\to c}_{0^+})",
                        Unit::None(),
                        &OPEUnitarityBounds::bound_0p),

                make_observable("b->c::Bound[0^-]@OPE", R"(B^{b\to c}_{0^-})",
                        Unit::None(),
                        &OPEUnitarityBounds::bound_0m),

                make_observable("b->c::Bound[1^+]@OPE", R"(B^{b\to c}_{1^+})",
                        Unit::None(),
                        &OPEUnitarityBounds::bound_1p),

                make_observable("b->c::Bound[1^-]@OPE", R"(B^{b\to c}_{1^-})",
                        Unit::None(),
                        &OPEUnitarityBounds::bound_1m),

                make_observable("b->c::Bound[0^+]@BGL", R"(B^{b\to c}_{0^+})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_0p),

                make_observable("b->c::Bound[0^-]@BGL", R"(B^{b\to c}_{0^-})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_0m),

                make_observable("b->c::Bound[1^+]@BGL", R"(B^{b\to c}_{1^+})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_1p),

                make_observable("b->c::Bound[1^-]@BGL", R"(B^{b\to c}_{1^-})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_1m),

                make_observable("b->c::Prior[0^+]@BGL", R"(B^{b\to c}_{0^+})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_0p_prior),

                make_observable("b->c::Prior[0^-]@BGL", R"(B^{b\to c}_{0^-})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_0m_prior),

                make_observable("b->c::Prior[1^+]@BGL", R"(B^{b\to c}_{1^+})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_1p_prior),

                make_observable("b->c::Prior[1^-]@BGL", R"(B^{b\to c}_{1^-})",
                        Unit::None(),
                        &BGLUnitarityBounds::bound_1m_prior),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_form_factors_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Form factors",
            "",
            {
                // B -> P
                make_b_to_pi_form_factors_group(),
                make_b_to_k_form_factors_group(),
                make_b_to_d_form_factors_group(),

                // B_s -> P
                make_bs_to_k_form_factors_group(),
                make_bs_to_ds_form_factors_group(),

                // B -> V
                make_b_to_omega_form_factors_group(),
                make_b_to_rho_form_factors_group(),
                make_b_to_kstar_form_factors_group(),
                make_b_to_dstar_form_factors_group(),

                // B_s -> V
                make_bs_to_kstar_form_factors_group(),
                make_bs_to_phi_form_factors_group(),
                make_bs_to_dsstar_form_factors_group(),

                // B -> P P
                make_b_to_pi_pi_form_factors_group(),

                // Lb -> 1/2^+
                make_lambdab_to_lambda_form_factors_group(),
                make_lambdab_to_lambdac_form_factors_group(),

                // unitarity bounds
                make_unitarity_bounds_group(),
            }
        );

        return ObservableSection(imp);
    }
}
