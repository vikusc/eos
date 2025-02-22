#include <test/test.hh>
#include <eos/rare-b-decays/b-to-k-charmonium.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

class BToKCharmoniumGvDV2020Test :
    public TestCase
{
    public:
    BToKCharmoniumGvDV2020Test() :
            TestCase("b_to_k_charmonium_GvDV2020_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();
            p["CKM::abs(V_cb)"]                          =  0.041996951916414726;
            p["CKM::arg(V_cb)"]                          =  0.0;
            p["CKM::abs(V_cs)"]                          =  0.9734061815416853;
            p["CKM::arg(V_cs)"]                          = -3.304199362533668e-05;
            p["mass::B_d"]                               =  5.279;
            p["mass::K_d"]                               =  0.492;
            p["mass::J/psi"]                             =  3.0969;
            p["mass::psi(2S)"]                           =  3.6860;
            p["mass::D^0"]                               =  1.86723;
            p["b->sccbar::t_0"]                          =  4.0;
            p["b->sccbar::t_s"]                          =  -17.4724;
            p["b->sccbar::chiOPE@GvDV2020"]              =  1.81e-4;

            p["B->Kccbar::Re{alpha_0^plus}@GvDV2020"]  = 2.0;
            p["B->Kccbar::Im{alpha_0^plus}@GvDV2020"]  = 3.0;
            p["B->Kccbar::Re{alpha_1^plus}@GvDV2020"]  = 4.0;
            p["B->Kccbar::Im{alpha_1^plus}@GvDV2020"]  = 5.0;
            p["B->Kccbar::Re{alpha_2^plus}@GvDV2020"]  = 6.0;
            p["B->Kccbar::Im{alpha_2^plus}@GvDV2020"]  = 7.0;

            Options oo
            {
                {"model",               "WET"},
                {"q",                   "d"},
                {"nonlocal-formfactor", "GvDV2020"},
                {"psi",                 "J/psi"}
            };

            BToKCharmonium c(p, oo);

            TEST_CHECK_RELATIVE_ERROR(c.branching_ratio(),  109982.6197301, 1e-5);

        }
} b_to_k_charmonium_GvDV2020_test;


class BToKCharmoniumGRvDV2021Test :
    public TestCase
{
    public:
    BToKCharmoniumGRvDV2021Test() :
            TestCase("b_to_k_charmonium_GRvDV2021_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();
            p["CKM::abs(V_cb)"]                           =  0.041996951916414726;
            p["CKM::arg(V_cb)"]                           =  0.0;
            p["CKM::abs(V_cs)"]                           =  0.9734061815416853;
            p["CKM::arg(V_cs)"]                           = -3.304199362533668e-05;
            p["mass::B_d"]                                = 5.279;
            p["mass::K_d"]                                = 0.492;
            p["mass::J/psi"]                              = 3.0969;
            p["mass::psi(2S)"]                            = 3.6860;
            p["mass::D^0"]                                = 1.86723;
            p["b->sccbar::t_0"]                           = 4.0;
            p["b->sccbar::t_s"]                           = -17.4724;
            p["b->sccbar::chiOPE@GRvDV2021"]              = 1.81e-4;

            p["B->Kccbar::Re{alpha_0^plus}@GRvDV2021"]  = 2.0;
            p["B->Kccbar::Im{alpha_0^plus}@GRvDV2021"]  = 3.0;
            p["B->Kccbar::Re{alpha_1^plus}@GRvDV2021"]  = 4.0;
            p["B->Kccbar::Im{alpha_1^plus}@GRvDV2021"]  = 5.0;
            p["B->Kccbar::Re{alpha_2^plus}@GRvDV2021"]  = 6.0;
            p["B->Kccbar::Im{alpha_2^plus}@GRvDV2021"]  = 7.0;

            Options oo
            {
                {"model",               "WET"},
                {"q",                   "d"},
                {"nonlocal-formfactor", "GRvDV2021"},
                {"psi",                 "J/psi"}
            };

            BToKCharmonium c(p, oo);

            TEST_CHECK_RELATIVE_ERROR(c.branching_ratio(),  756.3214, 1e-5);

        }
} b_to_k_charmonium_GRvDV2021_test;