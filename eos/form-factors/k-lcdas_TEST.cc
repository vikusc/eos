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

#include <test/test.hh>
#include <eos/form-factors/k-lcdas.hh>

#include <eos/utils/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class KaonLCDAsTest :
    public TestCase
{
    public:
        KaonLCDAsTest() :
            TestCase("k_lcdas_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["QCD::alpha_s(MZ)"]    =  0.1176;
            p["mass::s(2GeV)"]       =  0.095;
            p["K::a1@1GeV"]          = -0.07;
            p["K::a2@1GeV"]          = +0.24;
            p["mass::K_u"]           =  0.49368;
            p["decay-constant::K_u"] =  0.1561;

            /* Diagnostics */
            {
                KaonLCDAs k(p, Options{ });
                Diagnostics diagnostics = k.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(+1.00000, 1e-5), // c_rge(mu = 1.0 GeV)
                    std::make_pair(+0.94850, 1e-5), // c_rge(mu = 2.0 GeV)
                    std::make_pair(+0.92874, 1e-5), // c_rge(mu = 3.0 GeV)
                    std::make_pair(+0.91708, 1e-5), // c_rge(mu = 4.0 GeV)
                    std::make_pair(+0.90893, 1e-5), // c_rge(mu = 5.0 GeV)
                };

                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
            }

            /* Twist 2 */
            {
                KaonLCDAs k(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.17,      k.a2(1.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.126731,  k.a2(2.0),   eps);
                TEST_CHECK_NEARLY_EQUAL( 0.112741,  k.a2(3.0),   eps);

                // scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0,     k.phi(0.0, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.09617, k.phi(0.3, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.28625, k.phi(0.5, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.09617, k.phi(0.7, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,     k.phi(1.0, 1.0), eps);

                // scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0,     k.phi(0.0, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.14718, k.phi(0.3, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.32488, k.phi(0.5, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.14718, k.phi(0.7, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0,     k.phi(1.0, 2.0), eps);
            }

            /* Twist 3 */
            {
                KaonLCDAs k(p, Options{ });

                // coefficients at mu = 1.0 GeV, 2.0 GeV and 3.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.0045,         k.f3(1.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.003257533016, k.f3(2.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.002864248153, k.f3(3.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 1.846675085,    k.mu(1.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 2.434973113,    k.mu(2.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 2.697094626,    k.mu(3.0),        eps);
                TEST_CHECK_NEARLY_EQUAL( 0.01868720856,  k.eta3(1.0),      eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0102592843,   k.eta3(2.0),      eps);
                TEST_CHECK_NEARLY_EQUAL( 0.008143983155, k.eta3(3.0),      eps);
                TEST_CHECK_NEARLY_EQUAL(-1.5,            k.omega3(1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-1.124800783,    k.omega3(2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL(-1.002982205,    k.omega3(3.0),    eps);

                // phi3p, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.23828994,   k.phi3p(0.1, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.9881149354, k.phi3p(0.2, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.8447373277, k.phi3p(0.3, 1.0),    eps);

                // phi3p, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 1.133511907,  k.phi3p(0.1, 2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.9981866083, k.phi3p(0.2, 2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 0.9160656408, k.phi3p(0.3, 2.0),    eps);

                // phi3s, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.7314784825, k.phi3s(0.1, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.083784069,  k.phi3s(0.2, 1.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.219383352,  k.phi3s(0.3, 1.0),    eps);

                // phi3s, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 0.6416920521, k.phi3s(0.1, 2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.025740317,  k.phi3s(0.2, 2.0),    eps);
                TEST_CHECK_NEARLY_EQUAL( 1.238428959,  k.phi3s(0.3, 2.0),    eps);

                // phi3s first derivative, scale mu = 1.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 5.109460174,  k.phi3s_d1(0.1, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.207429218,  k.phi3s_d1(0.2, 1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.6979690447, k.phi3s_d1(0.3, 1.0), eps);

                // phi3s first derivative, scale mu = 2.0 GeV
                TEST_CHECK_NEARLY_EQUAL( 4.964350791,  k.phi3s_d1(0.1, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.860421439,  k.phi3s_d1(0.2, 2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.496070648,  k.phi3s_d1(0.3, 2.0), eps);
            }
        }
} k_lcdas_test;
