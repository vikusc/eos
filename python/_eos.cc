/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2016, 2019, 2020 Danny van Dyk
 * Copyright (c) 2021 Philip Lüghausen
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

#include "config.h"

#include "eos/constraint.hh"
#include "eos/observable.hh"
#include "eos/reference.hh"
#include "eos/signal-pdf.hh"
#include "eos/models/model.hh"
#include "eos/utils/kinematic.hh"
#include "eos/utils/log.hh"
#include "eos/utils/parameters.hh"
#include "eos/utils/options.hh"
#include "eos/utils/qualified-name.hh"
#include "eos/utils/reference-name.hh"
#include "eos/utils/units.hh"
#include "eos/statistics/goodness-of-fit.hh"
#include "eos/statistics/log-likelihood.hh"
#include "eos/statistics/log-posterior.hh"
#include "eos/statistics/log-prior.hh"
#include "eos/statistics/test-statistic-impl.hh"

#include "eos/rare-b-decays/charm-loops-impl.hh"

#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>

using namespace boost::python;
using namespace eos;

namespace impl
{
    // raw constructor for class Kinematics
    object
    Kinematics_ctor(tuple args, dict kwargs)
    {
        // strip off self
        object self = args[0];
        args = tuple(args.slice(1,_));

        self.attr("__init__")();

        if (1 < len(args))
        {
            PyErr_SetString(PyExc_TypeError, "eos.Kinematics expects exactly one argument, or keyword arguments, but not both");
            return object();
        }

        dict kinematics;

        if (1 == len(args))
        {
            kinematics = dict(args[0]);
            args = tuple(args.slice(1,_));
        }
        else
        {
            kinematics = kwargs;
        }

        list items = kinematics.items();
        for (unsigned i = 0 ; i < len(items) ; ++i)
        {
            object name = items[i][0];
            object value = items[i][1];
            self.attr("declare")(name, value);
        }

        return object();
    }

    // raw constructor for class Options
    object
    Options_ctor(tuple args, dict kwargs)
    {
        // strip off self
        object self = args[0];
        args = tuple(args.slice(1,_));

        self.attr("__init__")();

        if (1 < len(args))
        {
            PyErr_SetString(PyExc_TypeError, "eos.Options expects exactly one argument, or keyword arguments, but not both");
            return object();
        }

        dict options;

        if (1 == len(args))
        {
            options = dict(args[0]);
            args = tuple(args.slice(1,_));
        }
        else
        {
            options = kwargs;
        }

        list items = options.items();
        for (unsigned i = 0 ; i < len(items) ; ++i)
        {
            object name = items[i][0];
            object value = items[i][1];
            self.attr("declare")(name, value);
        }

        return object();
    }

    // converter for std::pair
    // converts a std::pair instance to a Python tuple, from Boost Python example
    template <typename T1, typename T2>
    struct std_pair_to_tuple
    {
        static PyObject* convert(std::pair<T1, T2> const& p)
        {
            return boost::python::incref(
                    boost::python::make_tuple(p.first, p.second).ptr());
        }
        static PyTypeObject const *get_pytype () {return &PyTuple_Type; }
    };

    // Helper for convenience.
    template <typename T1, typename T2>
    struct std_pair_to_python_converter
    {
        std_pair_to_python_converter()
        {
          boost::python::to_python_converter<
              std::pair<T1, T2>,
              std_pair_to_tuple<T1, T2>,
              true //std_pair_to_tuple has get_pytype
              >();
        }
    };

    // converter for std::vector
    // converts a std::vector instance to a Python list
    template <typename T>
    struct std_vector_to_list
    {
        static PyObject* convert(const std::vector<T> & v)
        {
            boost::python::list ret;
            for (const T & e : v)
            {
                ret.append(e);
            }

		    return incref(ret.ptr());
        }
        static PyTypeObject const * get_pytype () { return &PyList_Type; }
    };

    // Helper for convenience.
    template <typename T>
    struct std_vector_to_python_converter
    {
        std_vector_to_python_converter()
        {
          boost::python::to_python_converter<
              std::vector<T>,
              std_vector_to_list<T>,
              true //std_vector_to_list has get_pytype
              >();
        }
    };

    static const char version[] = PACKAGE_VERSION;

    void translate_exception(const Exception & e)
    {
        PyErr_SetString(PyExc_RuntimeError, e.what());
    }

    void logging_callback(PyObject * c, const std::string & id, const LogLevel & l, const std::string & m)
    {
        call<void>(c, id, l, m);
    }

    void register_log_callback(PyObject * c)
    {
        Log::instance()->register_callback(std::bind(&logging_callback, c, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
    }

    void set_native_log_level(const LogLevel & log_level)
    {
        Log::instance()->set_log_level(log_level);
    }

    // for testing
    void emit_native_log(const std::string & id, const LogLevel & log_level, const std::string & m)
    {
        Log::instance()->message(id, log_level) << m;
    }

    // wrapper to avoid issues with virtual inheritance and overloading
    double m_b_pole_wrapper_noargs(Model& m) { return m.m_b_pole(); }
}

BOOST_PYTHON_MODULE(_eos)
{
    using namespace boost::python;
    using namespace eos;

    // enable manually defined docstrings and python signatures, but disable
    // automatically generated C++ signatures.
    docstring_options local_docstring_options(true, true, false);

    // eos::Exception
    register_exception_translator<Exception>(&::impl::translate_exception);

    // native eos logging: provide functions and enum type
    def("_register_log_callback", &::impl::register_log_callback);
    def("_emit_native_log", &::impl::emit_native_log);
    def("_set_native_log_level", &::impl::set_native_log_level);
    enum_<LogLevel>("_NativeLogLevel")
        .value("SILENT",  ll_silent)
        .value("ERROR",   ll_error)
        .value("WARNING", ll_warning)
        .value("INFO",    ll_informational)
        .value("DEBUG",   ll_debug)
        ;

    // {{{ eos/utils
    // qnp::Prefix
    class_<qnp::Prefix>("qnpPrefix", init<std::string>())
        .def("__repr__", &qnp::Prefix::str, return_value_policy<copy_const_reference>())
        .def("__str__", &qnp::Prefix::str, return_value_policy<copy_const_reference>())
        .def("__lt__", &qnp::Prefix::operator<)
        ;

    // qnp::Name
    class_<qnp::Name>("qnpName", init<std::string>())
        .def("__repr__", &qnp::Name::str, return_value_policy<copy_const_reference>())
        .def("__str__", &qnp::Name::str, return_value_policy<copy_const_reference>())
        .def("__lt__", &qnp::Name::operator<)
        ;

    // qnp::Suffix
    class_<qnp::Suffix>("qnpSuffix", init<std::string>())
        .def("__repr__", &qnp::Suffix::str, return_value_policy<copy_const_reference>())
        .def("__str__", &qnp::Suffix::str, return_value_policy<copy_const_reference>())
        .def("__lt__", &qnp::Suffix::operator<)
        ;

    // QualifiedName
    class_<QualifiedName>("QualifiedName", R"(
            Represent a qualified (i.e. complete and syntactically correct) name.

            EOS uses qualified names when naming any observable or constraint. The
            composition is approximately::

                PREFIX::NAME@SUFFIX;OPTIONS

        )", init<std::string>())
        .def("__repr__", &QualifiedName::full, return_value_policy<copy_const_reference>())
        .def("__str__", &QualifiedName::str, return_value_policy<copy_const_reference>())
        .def("__eq__", &QualifiedName::operator==)
        .def("__ne__", &QualifiedName::operator!=)
        .def("__lt__", &QualifiedName::operator<)
        .def("prefix_part", &QualifiedName::prefix_part, return_value_policy<copy_const_reference>(), R"(
            Returns the prefix part of the name, i.e., the part preceeding the '::'.
        )")
        .def("name_part", &QualifiedName::name_part, return_value_policy<copy_const_reference>(), R"(
            Returns the name part of the name, i.e., the part following the '::' and preceeding any
            optional '@'.
        )")
        .def("suffix_part", &QualifiedName::suffix_part, return_value_policy<copy_const_reference>(), R"(
            Returns the optional suffix part of the name, i.e., the part following the optional '@'.
        )")
        ;
    implicitly_convertible<std::string, QualifiedName>();


    // ParameterSection
    class_<ParameterSection>("ParameterSection", no_init)
        .def("__iter__", range(&ParameterSection::begin, &ParameterSection::end))
        .def("name", &ParameterSection::name, return_value_policy<copy_const_reference>())
        .def("description", &ParameterSection::description, return_value_policy<copy_const_reference>())
        ;

    // ParameterGroup
    class_<ParameterGroup>("ParameterGroup", no_init)
        .def("__iter__", range(&ParameterGroup::begin, &ParameterGroup::end))
        .def("name", &ParameterGroup::name, return_value_policy<copy_const_reference>())
        .def("description", &ParameterGroup::description, return_value_policy<copy_const_reference>())
        ;

    // Parameters
    class_<Parameters>("_Parameters", no_init)
        .def("Defaults", &Parameters::Defaults)
        .staticmethod("Defaults")
        .def("__getitem__", (Parameter (Parameters::*)(const QualifiedName &) const) &Parameters::operator[])
        .def("by_id", (Parameter (Parameters::*)(const Parameter::Id &) const) &Parameters::operator[])
        .def("__iter__", range(&Parameters::begin, &Parameters::end))
        .def("declare", &Parameters::declare, return_value_policy<return_by_value>())
        .def("sections", range(&Parameters::begin_sections, &Parameters::end_sections))
        .def("set", &Parameters::set,
            R"(
            Set the value of a parameter.

            :param name: The name of the parameter to set.
            :type name: str
            :param value: The value to set the parameter to.
            :type value: float
            )")
        .def("has", &Parameters::has)
        .def("override_from_file", &Parameters::override_from_file)
        ;

    // Parameter
    class_<Parameter>("Parameter", R"(
            Represents a single real-valued scalar parameter in EOS.

            Users cannot directly create new objects of this class. However, new named sets of parameters can be created,
            via the :class:`eos.Parameters <eos.Parameters>` class, from which the parameter of interest can be extracted, inspected, and altered.
        )", no_init)
        .def(float_(self))
        .def("central", &Parameter::central, return_value_policy<copy_const_reference>())
        .def("max", &Parameter::max, return_value_policy<copy_const_reference>())
        .def("min", &Parameter::min, return_value_policy<copy_const_reference>())
        .def("name", &Parameter::name, return_value_policy<copy_const_reference>(),
            R"(
            Returns the name of the parameter.
            )")
        .def("latex", &Parameter::latex, return_value_policy<copy_const_reference>(),
            R"(
            Returns the LaTeX representation of the parameter.
            )")
        .def("unit", &Parameter::unit)
        .def("set", &Parameter::set,
            R"(
            Set the value of a parameter.

            :param value: The value to set the parameter to.
            :type value: float
            )")
        .def("set_max", &Parameter::set_max)
        .def("set_min", &Parameter::set_min)
        .def("evaluate", &Parameter::evaluate,
            R"(
            Return the current value of a parameter.
            )")
        ;

    // ParameterUser
    class_<ParameterUser>("ParameterUser", no_init)
        .def("used_parameter_ids", range(&ParameterUser::begin, &ParameterUser::end))
        ;

    // ParameterRange
    class_<ParameterRange>("ParameterRange", init<double, double>())
        ;

    // Kinematics
    class_<Kinematics>("Kinematics", R"(
            Represents the set of kinematic variables relevant to an :class:`observable <eos.Obserable>`.

            Initialize a new set of kinematic variables. The inital set of variables and their initial
            set of values can be provided through keyword arguments, e.g. using

            .. code-block::

               k = eos.Kinematics(q2=0.4, k2=0.0)                      # default keyword arguments
               k = eos.Kinematics({'q2': 0.4, 'cos(theta_l)': -1.0})   # use a dictionary if variable names are not
                                                                       # valid python identifiers
        )", no_init)
        .def("__init__", raw_function(&::impl::Kinematics_ctor))
        .def(init<>())
        .def("__add__", &Kinematics::operator+)
        .def("__iter__", range(&Kinematics::begin, &Kinematics::end))
        .def("__getitem__", (KinematicVariable (Kinematics::*)(const std::string &) const) &Kinematics::operator[])
        .def("declare", &Kinematics::declare, return_value_policy<return_by_value>(), R"(
            Declares a new kinematic variable.

            :param name: The name of the new kinematic variable.
            :type name: str
            :param value: The initial value for the new kinematic variable.
            :type value: float
        )", args("self", "name", "value"))
        .def("__str__", &Kinematics::as_string)
        ;

    // KinematicVariable
    class_<KinematicVariable>("KinematicVariable", no_init)
        .def(float_(self))
        .def("name", &KinematicVariable::name, return_value_policy<copy_const_reference>())
        .def("set", &KinematicVariable::set)
        .def("evaluate", &KinematicVariable::evaluate)
        ;

    // Options
    ::impl::std_pair_to_python_converter<const std::string, std::string> converter_options_iter;
    class_<Options>("Options", R"(
            Represents the set of options provided to an observable.

            Options are pairs of (key, value) pairs. The list of valid keys and their
            respective valid options are specific to each observable. The initialization
            accepts keyword arguments, e.g.:

            .. code-block::

               o = eos.Options(model='WET')                   # default keyword arguments
               o = eos.Options({'form-factors': 'BSZ2015'})   # use a dictionary if option keys are not
                                                              # valid python identifiers
        )", no_init)
        .def("__init__", raw_function(&::impl::Options_ctor))
        .def("__iter__", range(&Options::begin, &Options::end))
        .def(init<>())
        .def("declare", &Options::declare)
        .def("__str__", &Options::as_string)
        ;

    // OptionSpecification
    ::impl::std_vector_to_python_converter<std::string> converter_option_specifications;
    class_<OptionSpecification>("OptionSpecification")
        .def_readonly("key", &OptionSpecification::key)
        .add_property("allowed_values", make_getter(&OptionSpecification::allowed_values, return_value_policy<return_by_value>()))
        .def_readonly("default_value", &OptionSpecification::default_value)
        ;

    // Units
    class_<Unit>("Unit", R"(
            Represents the unit of the observables.

            Seven possible entries are currently implemented in EOS:
            - Undefined
            - None
            - GeV
            - GeV2
            - InverseGeV2
            - InverseGeV4
            - InversePicoSecond
        )", no_init)
        .def("Undefined", &Unit::Undefined)
        .staticmethod("Undefined")
        .def("Unity", &Unit::None)
        .staticmethod("Unity")
        .def("GeV", &Unit::GeV)
        .staticmethod("GeV")
        .def("GeV2", &Unit::GeV2)
        .staticmethod("GeV2")
        .def("InverseGeV2", &Unit::InverseGeV2)
        .staticmethod("InverseGeV2")
        .def("InverseGeV4", &Unit::InverseGeV4)
        .staticmethod("InverseGeV4")
        .def("InversePicoSecond", &Unit::InversePicoSecond)
        .staticmethod("InversePicoSecond")
        .def("latex", &Unit::latex, return_value_policy<copy_const_reference>())
        .def("__eq__", &Unit::operator==)
        ;

    // WilsonCoefficients
    class_<WilsonCoefficients<eos::BToS>>("BToSWilsonCoefficients", no_init)
        .def("c1", &WilsonCoefficients<eos::BToS>::c1)
        .def("c2", &WilsonCoefficients<eos::BToS>::c2)
        ;

    // Model
    register_ptr_to_python<std::shared_ptr<Model>>();
    class_<Model, boost::noncopyable>("Model", no_init)
        .def("make", &Model::make, return_value_policy<return_by_value>())
        .staticmethod("make")
        // CKM component
        .def("ckm_cd", &Model::ckm_cd)
        .def("ckm_cs", &Model::ckm_cs)
        .def("ckm_cb", &Model::ckm_cb)
        .def("ckm_ud", &Model::ckm_ud)
        .def("ckm_us", &Model::ckm_us)
        .def("ckm_ub", &Model::ckm_ub)
        .def("ckm_td", &Model::ckm_td)
        .def("ckm_ts", &Model::ckm_ts)
        .def("ckm_tb", &Model::ckm_tb)
        // QCD component
        .def("m_t_msbar",  &Model::m_t_msbar)
        .def("m_t_pole",   &Model::m_t_pole)
        .def("m_b_kin",    &Model::m_b_kin)
        .def("m_b_msbar",  &Model::m_b_msbar)
        .def("m_b_pole",   &Model::m_b_pole)
        .def("m_b_pole",   &::impl::m_b_pole_wrapper_noargs)
        .def("m_c_kin",    &Model::m_c_kin)
        .def("m_c_msbar",  &Model::m_c_msbar)
        .def("m_c_pole",   &Model::m_c_pole)
        .def("m_s_msbar",  &Model::m_s_msbar)
        .def("m_ud_msbar", &Model::m_ud_msbar)
        // WilsonCoefficients
        .def("wilson_coefficients_b_to_s", &Model::wilson_coefficients_b_to_s)
        // alpha_s
        .def("alpha_s", &Model::alpha_s)
        ;

    // ObservableCache
    class_<ObservableCache>("ObservableCache", no_init)
        .def("__iter__", range(&ObservableCache::begin, &ObservableCache::end))
        .def("__getitem__", &ObservableCache::operator[])
        .def("add", &ObservableCache::add)
        .def("update", &ObservableCache::update)
        ;

    // ReferenceName
    class_<ReferenceName>("ReferenceName", init<std::string>())
        .def("__str__", &ReferenceName::str, return_value_policy<copy_const_reference>())
        .def("__eq__", &ReferenceName::operator==)
        .def("__ne__", &ReferenceName::operator!=)
        .def("__lt__", &ReferenceName::operator<)
        ;
    implicitly_convertible<std::string, ReferenceName>();

    // }}}

    // {{{ eos/statistics
    // LogLikelihoodBlock
    register_ptr_to_python<std::shared_ptr<LogLikelihoodBlock>>();
    class_<LogLikelihoodBlock, boost::noncopyable>("LogLikelihoodBlock", no_init)
        .def("__str__", &LogLikelihoodBlock::as_string)
        ;

    // LogLikelihood
    class_<LogLikelihood>("LogLikelihood", R"(
            Represents the log(likelihood) of a Bayesian analysis undertaken with the :class:`Analysis <eos.Analysis>` class.
        )", init<Parameters>())
        .def("add", (void (LogLikelihood::*)(const Constraint &)) &LogLikelihood::add)
        .def("__iter__", range(&LogLikelihood::begin, &LogLikelihood::end))
        .def("observable_cache", &LogLikelihood::observable_cache)
        .def("evaluate", &LogLikelihood::operator())
        ;

    // Constraint
    class_<Constraint>("Constraint", no_init)
        .def("make", &Constraint::make, return_value_policy<return_by_value>())
        .staticmethod("make")
        .def("name", &Constraint::name, return_value_policy<copy_const_reference>())
        .def("blocks", range(&Constraint::begin_blocks, &Constraint::end_blocks))
        .def("observables", range(&Constraint::begin_observables, &Constraint::end_observables))
        ;

    // ConstraintEntry
    register_ptr_to_python<std::shared_ptr<const ConstraintEntry>>();
    class_<eos::ConstraintEntry, boost::noncopyable>("ConstraintEntry", no_init)
        .def("make", &ConstraintEntry::make, return_value_policy<return_by_value>())
        .def("name", &ConstraintEntry::name, return_value_policy<copy_const_reference>())
        .def("type", &ConstraintEntry::type, return_value_policy<copy_const_reference>())
        .def("observables", range(&ConstraintEntry::begin_observable_names, &ConstraintEntry::end_observable_names))
        .def("serialize", (std::string (ConstraintEntry::*)(void) const) &ConstraintEntry::serialize, return_value_policy<return_by_value>())
        .def("deserialize", (ConstraintEntry * (*)(const QualifiedName &, const std::string &)) &ConstraintEntry::FromYAML, return_value_policy<manage_new_object>())
        .staticmethod("deserialize")
        ;

    // Constraints
    ::impl::std_pair_to_python_converter<const QualifiedName, std::shared_ptr<const ConstraintEntry>> converter_constraints_iter;
    class_<Constraints>("_Constraints")
        .def("__getitem__", (std::shared_ptr<const ConstraintEntry> (Constraints::*)(const QualifiedName &) const) &Constraints::operator[])
        .def("__iter__", range(&Constraints::begin, &Constraints::end))
        .def("insert", &Constraints::insert)
        ;

    // LogPrior
    register_ptr_to_python<std::shared_ptr<LogPrior>>();
    class_<LogPrior, boost::noncopyable>("LogPrior", R"(
            Represents a Bayesian prior on the log scale.

            New LogPrior objects can only be created using the capitalized static methods:
            :meth:`LogPrior.Uniform`, :meth:`LogPrior.Gaussian`, and :meth:`LogPrior.Scale`.
        )", no_init)
        .def("Uniform", &LogPrior::Flat, return_value_policy<return_by_value>(), R"(
            Returns a new uniform prior as a LogPrior.

            The prior's support is provided by the `range` parameter.

            :param parameters: The parameters to which this LogPrior is bound.
            :type parameters: eos.Parameters
            :param name: The name of the parameter for which the LogPrior is defined.
            :type name: str
            :param range: The range [min, max] for the values that the parameter is allowed to take.
            :type range: tuple of two floating point numbers
        )", args("parameters", "name", "range"))
        .staticmethod("Uniform")
        .def("Flat", &LogPrior::Flat, return_value_policy<return_by_value>(), "Alias for :meth:`LogPrior.Uniform`.",
            args("parameters", "name", "range"))
        .staticmethod("Flat")
        .def("Gauss", &LogPrior::Gauss, return_value_policy<return_by_value>(), R"(
            Returns a new Gaussian prior as a LogPrior.

            The prior's support is provided by the `range` parameter, with the
            68% probability interval [`lower`, `upper`] and the mode provided
            by the parameter `central`.

            :param parameters: The parameters to which this LogPrior is bound.
            :type parameters: eos.Parameters
            :param name: The name of the parameter for which the LogPrior is defined.
            :type name: str
            :param range: The range [min, max] for the values that the parameter is allowed to take.
            :type range: tuple of two floating point numbers
            :param lower: The lower boundary of the 68% probability interval.
            :type lower: float
            :param central: The mode and median of the prior.
            :type central: float
            :param upper: The upper boundary of the 68% probability interval.
            :type upper: float
        )", args("parameters", "name", "range", "lower", "central", "upper"))
        .staticmethod("Gauss")
        .def("Scale", &LogPrior::Scale, return_value_policy<return_by_value>(), R"(
            Returns a new Scale prior as a LogPrior.

            The prior's support is provided by the `range` parameter, which should
            coincide with [`mu_0 / lambda`, `mu_0 * lambda`]. The PDF is chosen
            such that a renormalization scale is varied in this range and with
            central value `mu_0` such that :math:`\ln x / \mu_0` is uniformly
            distributed in the interval :math:`[-\ln \lambda, +\ln \lambda]`.

            :param parameters: The parameters to which this LogPrior is bound.
            :type parameters: eos.Parameters
            :param name: The name of the parameter for which the LogPrior is defined.
            :type name: str
            :param range: The range [min, max] for the values that the parameter is allowed to take.
            :type range: tuple of two floating point numbers
            :param mu_0: The central value of the parameter.
            :type mu_0: float, strictly positive
            :param lambda: The scale factor.
            :type lambda: float, strictly positive
        )", args("parameters", "name", "range", "mu_0", "scale"))
        .staticmethod("Scale")
        .def("inverse_cdf", &LogPrior::inverse_cdf, R"(
            Returns the parameter value corresponding to the cumulative propability :math:`p`.

            :param p: The cumulative propability.
            :type p: float, [0.0, 1.0]
        )", args("p"))
        ;

    // LogPosterior
    class_<LogPosterior>("LogPosterior", init<LogLikelihood>())
        .def("add", &LogPosterior::add)
        .def("log_likelihood", &LogPosterior::log_likelihood)
        .def("log_priors", range(&LogPosterior::begin_priors, &LogPosterior::end_priors))
        .def("evaluate", &LogPosterior::evaluate)
        ;

    // test_statistics::ChiSquare
    class_<test_statistics::ChiSquare>("test_statisticsChiSquare", no_init)
        .def_readonly("chi2", &test_statistics::ChiSquare::chi2)
        .def_readonly("dof", &test_statistics::ChiSquare::dof)
        ;

    // GoodnessOfFit
    ::impl::std_pair_to_python_converter<const QualifiedName, test_statistics::ChiSquare> converter_goodnessoffit_chi_square_iter;
    class_<GoodnessOfFit>("GoodnessOfFit", R"(
            Represents the goodness of fit characteristics of the log(posterior).
        )", init<LogPosterior>())
        .def("__iter__", range(&GoodnessOfFit::begin_chi_square, &GoodnessOfFit::end_chi_square))
        .def("total_chi_square", &GoodnessOfFit::total_chi_square, R"(
            Returns the total :math:`\chi^2` value of the log(likelihood). Only (multivariate) gaussian
            likelihoods are considered for this result.
        )")
        .def("total_degrees_of_freedom", &GoodnessOfFit::total_degrees_of_freedom, R"(
            Returns the total number of degrees of freedom in the log(posterior).
        )")
        ;

    // }}}

    // {{{ eos/
    // Reference
    register_ptr_to_python<ReferencePtr>();
    class_<Reference>("Reference", no_init)
        .def("name", &Reference::name, return_value_policy<copy_const_reference>())
        .def("authors", &Reference::authors, return_value_policy<copy_const_reference>())
        .def("eprint_archive", &Reference::eprint_archive, return_value_policy<copy_const_reference>())
        .def("eprint_id", &Reference::eprint_id, return_value_policy<copy_const_reference>())
        .def("title", &Reference::title, return_value_policy<copy_const_reference>())
        .def("inspire_id", &Reference::inspire_id, return_value_policy<copy_const_reference>())
        ;

    // References
    ::impl::std_pair_to_python_converter<const ReferenceName, ReferencePtr> converter_references_iter;
    class_<References>("_References")
        .def("__getitem__", &References::operator[])
        .def("__iter__", range(&References::begin, &References::end))
        ;

    // ReferenceUser
    class_<ReferenceUser>("ReferenceUser", no_init)
        .def("references", range(&ReferenceUser::begin_references, &ReferenceUser::end_references))
        ;


    // Observable
    register_ptr_to_python<std::shared_ptr<Observable>>();
    class_<Observable, bases<ParameterUser, ReferenceUser>, boost::noncopyable>("Observable", R"(
            Represents an observable or pseudo observable known to EOS.

            New observable objects are created using the :meth:`make <eos.Observable.make>` static method.
            See also `the complete list of observables <../observables.html>`_.
        )", no_init)
        .def("make", &Observable::make, return_value_policy<return_by_value>(), R"(
            Makes a new :class:`Observable` object.

            :param name: The name of the observable. See `the complete list of observables <../observables.html>`_.
            :type name: eos.QualifiedName
            :param parameters: The set of parameters to which this observable is bound.
            :type parameters: eos.Parameters
            :param kinematics: The set of kinematic variables to which this observable is bound.
            :type kinematics: eos.Kinematics
            :param options: The set of options relevant to this observable.
            :type options: eos.Options

            :return: The new observable object.
            :rtype: eos.Observable
        )", args("name", "parameters", "kinematics", "options"))
        .staticmethod("make")
        .def("evaluate", &Observable::evaluate, R"(
            Evaluates the observable for the present values of its bound set of parameters and set of kinematic variables.

            :return: The value of the observable.
            :rtype: float
        )", args("self"))
        .def("name", &Observable::name, return_value_policy<copy_const_reference>(), R"(
            Returns the name of the observable.
        )")
        .def("parameters", &Observable::parameters, R"(
            Returns the set of parameters bound to this observable.
        )")
        .def("kinematics", &Observable::kinematics, R"(
            Returns the set of kinematic variables bound to this observable.
        )")
        .def("options", &Observable::options, R"(
            Returns the set of options used when creating the observable.
        )")
        ;

    // ObservableEntry
    register_ptr_to_python<std::shared_ptr<const ObservableEntry>>();
    class_<ObservableEntry, boost::noncopyable>("ObservableEntry", no_init)
        .def("name", &ObservableEntry::name, return_value_policy<copy_const_reference>())
        .def("latex", &ObservableEntry::latex, return_value_policy<copy_const_reference>())
        .def("unit", &ObservableEntry::unit, return_internal_reference<>())
        .def("kinematic_variables", range(&ObservableEntry::begin_kinematic_variables, &ObservableEntry::end_kinematic_variables))
        .def("options", range(&ObservableEntry::begin_options, &ObservableEntry::end_options))
        ;

    // ObservableGroup
    register_ptr_to_python<std::shared_ptr<ObservableGroup>>();
    class_<ObservableGroup>("ObservableGroup", no_init)
        .def("__iter__", range(&ObservableGroup::begin, &ObservableGroup::end))
        .def("name", &ObservableGroup::name, return_value_policy<copy_const_reference>())
        .def("description", &ObservableGroup::description, return_value_policy<copy_const_reference>())
        ;

    // ObservableSection
    register_ptr_to_python<std::shared_ptr<ObservableSection>>();
    class_<ObservableSection>("ObservableSection", no_init)
        .def("__iter__", range(&ObservableSection::begin, &ObservableSection::end))
        .def("name", &ObservableSection::name, return_value_policy<copy_const_reference>())
        .def("description", &ObservableSection::description, return_value_policy<copy_const_reference>())
        ;

    // Observables
    ::impl::std_pair_to_python_converter<const QualifiedName, ObservableEntryPtr> converter_observables_iter;
    class_<Observables>("_Observables")
        .def("__getitem__", &Observables::operator[])
        .def("__iter__", range(&Observables::begin, &Observables::end))
        .def("insert", &Observables::insert, R"(
            Insert a new observable to EOS by parsing the input string.

            :param name: The name of the new observable.
            :type name: eos.QualifiedName
            :param latex: The latex representation of the new observable.
            :type latex: std::string
            :param unit: The unit of the new observable.
            :type unit: eos.Unit
            :param options: The set of options relevant to this new observable. These options apply to **all** the observables in the expression.
            :type options: eos.Options
            :param expression: The expression to be parsed.
            :type expression: std::string
        )", args("name", "latex", "unit", "options", "expression"))
        .def("sections", range(&Observables::begin_sections, &Observables::end_sections))
        ;


    // SignalPDF
    register_ptr_to_python<std::shared_ptr<SignalPDF>>();
    class_<SignalPDF, boost::noncopyable>("_SignalPDF", R"(
            Represents a probability density function (PDF) for any of the physics signals known to EOS.

            New PDF objects are created using the :meth:`make <eos.SignalPDF.make>` static method.
            See also `the complete list of PDFs <../signal-pdfs.html>`_.
    )", no_init)
        .def("make", &SignalPDF::make, return_value_policy<return_by_value>()) // docstring is maintained in python/eos/signal_pdf.py
        .staticmethod("make")
        .def("evaluate", &SignalPDF::evaluate, R"(
            Evaluates the (unnormalized) PDF for the present values of the sets of parameters and kinematic variables that it is bound to.

            :return: The value of the PDF.
            :rtype: float
        )", args("self"))
        .def("normalization", &SignalPDF::normalization, R"(
            Evaluates the normalization of the PDF.

            To speed up sampling from the PDF, the :meth:`evaluate <eos.SignalPDF.evaluate>` returns values of an
            unnormalized function proportional to the actual PDF. To ensure that the integral over the PDF is
            normalized to 1, the values returned by evaluate need to be divided by the return value of this method.
        )")
        .def("name", &SignalPDF::name, return_value_policy<copy_const_reference>(), R"(
            Returns the name of the PDF.
        )")
        .def("parameters", &SignalPDF::parameters, R"(
            Returns the set of parameters bound to this PDF.
        )")
        .def("options", &SignalPDF::options, R"(
            Returns the set of options used when creating the PDF.
        )")
        .def("kinematics", &SignalPDF::kinematics, R"(
            Returns the set of kinematic variables bound to this PDF.
        )")
        ;

    // SignalPDFEntry
    register_ptr_to_python<std::shared_ptr<SignalPDFEntry>>();
    class_<SignalPDFEntry, boost::noncopyable>("SignalPDFEntry", no_init)
        .def("name", &SignalPDFEntry::name, return_value_policy<copy_const_reference>())
        .def("description", &SignalPDFEntry::description, return_value_policy<copy_const_reference>())
        ;

    // SignalPDFs
    ::impl::std_pair_to_python_converter<const QualifiedName, std::shared_ptr<SignalPDFEntry>> converter_signalpdfs_iter;
    class_<SignalPDFs>("_SignalPDFs")
        .def("__iter__", range(&SignalPDFs::begin, &SignalPDFs::end))
        ;

    // }}}

    // EOS version
    scope().attr("__version__") = PACKAGE_VERSION;

    // Analytic modification of C_7 and C_9 at NLO
    def("delta_c7", &agv_2019a::delta_c7);
    def("delta_c7_Qc", &agv_2019a::delta_c7_Qc);
    def("delta_c9", &agv_2019a::delta_c9);
    def("delta_c9_Qc", &agv_2019a::delta_c9_Qc);
}
