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
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <test/test.hh>
#include <eos/rare-b-decays/charm-loops.hh>

#include <cmath>

using namespace test;
using namespace eos;
using namespace agv_2019a;

class LoopParameterTest :
    public TestCase
{
    public:
        LoopParameterTest() :
            TestCase("loop_parameter_test")
        {
        }

        virtual void run() const
        {
            {
                static const double eps = 1e-12;
                agv_2019a::CharmLoopsParameters testclp = {/*muhat =*/ 1.0, /*s =*/ -4.0, /*z =*/ 0.15, /*feynepsilonhat*/ 1e-10};

                /* Check, that the square root in C++ of complex arguments is handled correctly */
                static const complex<double> m1(-1.0, 0.0);
                TEST_CHECK_NEARLY_EQUAL(0.0,     std::sqrt(m1).real(), eps);
                TEST_CHECK_NEARLY_EQUAL(1.0,     std::sqrt(m1).imag(), eps);

                /* Comparison with Mathematica results */

                TEST_CHECK_NEARLY_EQUAL(1.5811388300841895,     testclp.x_a.real(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.7745966689726898,     testclp.x_b.real(), eps);
                TEST_CHECK_NEARLY_EQUAL(1.5811388300841895,     testclp.x_c.real(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.7745966689726898,     testclp.x_d.real(), eps);
                TEST_CHECK_NEARLY_EQUAL(1.5811388300841895,     testclp.x_e.real(), eps);

                TEST_CHECK_NEARLY_EQUAL(-6.719840027857805e-10, testclp.x_a.imag(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.6324555318142068,     testclp.x_b.imag(), eps);
                TEST_CHECK_NEARLY_EQUAL(-6.719840027857805e-10, testclp.x_c.imag(), eps);
                TEST_CHECK_NEARLY_EQUAL(0.6324555318142068,     testclp.x_d.imag(), eps);
                TEST_CHECK_NEARLY_EQUAL(-6.719840027857805e-10, testclp.x_e.imag(), eps);
            }
        }
} loop_parameter_test ;
