#include <iostream>
#include <filesystem>

////////////////////////////////////////////////////////////////////////////////
#include <esl/computation/environment.hpp>
#include <esl/simulation/model.hpp>
#include <esl/economics/markets/quote.hpp>
#include <esl/law/jurisdictions.hpp>
#include <esl/economics/markets/walras/tatonnement.hpp>
#include <esl/economics/markets/walras/price_setter.hpp>
#include <esl/economics/markets/walras/quote_message.hpp>
#include <esl/economics/finance/bond.hpp>
#include <esl/data/file.hpp>
#include <esl/computation/thread_pool.hpp>

using namespace esl;
using namespace esl::computation;
using namespace esl::economics;
using namespace esl::economics::markets::walras;
using namespace esl::economics::finance;
using namespace esl::law;
using namespace esl::economics::currencies;
using namespace esl::simulation;
using namespace esl::simulation::parameter;

#include "traded_company.hpp"

////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::to_string;
using std::pair;
using std::vector;
using std::map;
using std::shared_ptr;
using std::make_shared;
using std::tuple;
using std::get;
using std::make_tuple;
using std::dynamic_pointer_cast;
////////////////////////////////////////////////////////////////////////////////

#include "strategies/constant_demand.hpp"
#include "strategies/fundamental_value/dividend_discount.hpp"
#include "strategies/fundamental_value/mean_reverting_noise.hpp"
#include "strategies/technical/momentum.hpp"
#include <experiment/experiment_1_population_fluctuations.hpp>
#include <experiment/experiment_3_simplex.hpp>
#include <experiment/experiment_4_trajectories.hpp>
#include <experiment/experiment_5_flows.hpp>
#include <experiment/experiment_6_kelly.hpp>

////////////////////////////////////////////////////////////////////////////////

#include <boost/program_options.hpp>
using namespace boost::program_options;

////////////////////////////////////////////////////////////////////////////////

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#endif

///
/// \param argument_count
/// \param arguments
/// \return
int main(int argument_count, char *arguments[])
{
    // on windows hosts, prefer to run simulations with below normal priority
    // because otherwise the host may become difficult to interact with
#if defined(WIN32) || defined(_WIN32) || (defined(__WIN32) && !defined(__CYGWIN__))
    SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
    // set preferred floating point format
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.setf(std::ios::showpoint);

    for(auto a = 0; a < argument_count; ++a) {
        LOG(notice) << arguments[a] << std::endl;
    }

    options_description description_("Allowed options");
    description_.add_options()
        ("help", "produce help message")
        ("experiment", value<std::string>()->default_value("experiment_6"), "choose experiment")
        ("reinvestment", value<double>()->default_value(1.), "reinvestment rate")
        ("stocks", value<unsigned int>()->default_value(1), "set the number of stocks")
        ;

    variables_map arguments_;
    store(parse_command_line(argument_count, arguments, description_), arguments_);

    notify(arguments_);

    if (arguments_.count("help")) {
        std::cout << description_ << std::endl;
        return 0;
    }

    unsigned int stocks_count_;
    if (arguments_.count("stocks")) {
        stocks_count_ = arguments_["stocks"].as<unsigned int>();
    }

    LOG(notice) << "stocks_count_ " << stocks_count_ << std::endl;

    std::string experiment_ = arguments_["experiment"].as<std::string>();
    std::transform(experiment_.begin(), experiment_.end(), experiment_.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if("experiment_3" == experiment_) {
        experiment_3();

    }else if("experiment_4" == experiment_){
        experiment_4();

    }else if("experiment_5" == experiment_){
        double reinvestment_rate = arguments_["reinvestment"].as<double>();
        experiment_5(reinvestment_rate);

    }else if("experiment_3_best" == experiment_){
        experiment_3_best(32);

    }else if("experiment_6" == experiment_){
        experiment_6_best();

    }else if("volatility_illustration" == experiment_){
        volatility_illustration(1);
    }

    return 0;
}