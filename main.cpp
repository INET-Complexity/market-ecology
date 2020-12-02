#include <iostream>
#include <filesystem>

////////////////////////////////////////////////////////////////////////////////
#include <esl/computation/environment.hpp>
#include <esl/simulation/model.hpp>
using esl::simulation::parameter::parametrization;
#include <esl/economics/markets/quote.hpp>
#include <esl/law/jurisdictions.hpp>

using esl::economics::currencies::USD;

#include <esl/economics/markets/walras/tatonnement.hpp>
#include <esl/economics/markets/walras/price_setter.hpp>
using esl::economics::markets::walras::price_setter;
#include <esl/economics/markets/walras/quote_message.hpp>

#include <esl/economics/finance/bond.hpp>
#include <esl/data/file.hpp>

#include <esl/computation/thread_pool.hpp>

using namespace esl;
using namespace esl::economics;
using namespace esl::economics::finance;
using namespace esl::law;

#include "traded_company.hpp"

using esl::simulation::time_duration;
using esl::economics::rate;
using esl::computation::environment;
using esl::simulation::model;

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

#include "strategies/fundamental_value/mean_reverting_noise.hpp"
#include "strategies/fundamental_value/dividend_discount.hpp"
#include "strategies/technical/momentum.hpp"
#include "strategies/constant_demand.hpp"

////////////////////////////////////////////////////////////////////////////////

#include <boost/program_options.hpp>
#include <experiment/e_3_grid.hpp>
using namespace boost::program_options;

////////////////////////////////////////////////////////////////////////////////

#include "experiment/e_6.hpp"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#endif

///
/// \brief main market ecology model
///
struct main_model
: public model
{
    using model::model;

    ///
    /// Maximum number of steps for the simulation
    ///
    constexpr static time_point end_time = 252 * 300;

    ///
    /// \param step
    /// \return
    time_point step(time_interval step) override
    {
        auto ps_ =  dynamic_pointer_cast<price_setter>(agents.local_agents_[identity<agent>{1}]);
        auto nt = dynamic_pointer_cast<fund>(agents.local_agents_[identity<agent>{2}]);
        auto fv = dynamic_pointer_cast<fund>(agents.local_agents_[identity<agent>{3}]);
        auto tf = dynamic_pointer_cast<fund>(agents.local_agents_[identity<agent>{4}]);

        // if time is larger than this, test whether the simulation is stuck possibly
        // this does not happen normally, but is used to warn the modeler on making changes to the code
        // that could potentially cause this
        constexpr size_t early_exit_after_ = 10000;
        if(is_stuck(early_exit_after_)){
            return end_time;
        }

        // let the agents act
        auto result_ = model::step(step);

        // multiple cases of resetting wealth or re-investment
        if(!dynamic_pointer_cast<fund>(agents.local_agents_[identity<agent>{2}])->output_net_asset_value->values.empty()){
            auto nt_ = double(get<1>(nt->output_net_asset_value->values.back()));
            auto fv_ = double(get<1>(fv->output_net_asset_value->values.back()));
            auto tf_ = double(get<1>(tf->output_net_asset_value->values.back()));

            auto total_ = nt_ + fv_ + tf_;
            // because total wealth changes and we want to reset to a percentage
            // of total rather than an absolute amount, we compute the reset
            // amount every time step
            auto reset_amount_ = total_ * 0.001;
            double nw[3] = {nt_, fv_, tf_};
            std::shared_ptr<fund> agents_[3] = {nt, fv, tf};

            for(size_t i = 0; i < 3; ++i){
                // if one agent's wealth exceeds 99%
                if(nw[i] / total_ > 0.99){
                    for(auto &[p, q]: agents_[i]->inventory) {
                        auto cast_ = dynamic_pointer_cast<cash>(p);
                        if(cast_) {
                            agents_[i]->inventory[p].amount -= std::uint64_t(total_ * 0.01 * 100);
                        }
                    }
                }
                // if one agent's wealth is less than 1%
                if(nw[i] / total_ < 0.01 && agents_[i]->target_net_asset_value && double(agents_[i]->target_net_asset_value.value())> 0.) {
                    LOG(notice) << step.lower << " " << i << " falls below 1%" << std::endl;
                    for(auto &[p, q] :  agents_[i]->inventory) {
                        auto cast_ = dynamic_pointer_cast<cash>(p);
                        if(cast_) {
                            agents_[i]->inventory[p].amount +=
                                std::uint64_t(total_ * 0.01 * 100);
                        }
                    }
                }
            }
            nt->reset_amount = price::approximate(reset_amount_, currencies::USD);
            fv->reset_amount = price::approximate(reset_amount_, currencies::USD);
            tf->reset_amount = price::approximate(reset_amount_, currencies::USD);
        }else{
            nt->reset_amount = price::approximate(100'000.00, currencies::USD);
            fv->reset_amount = price::approximate(100'000.00, currencies::USD);
            tf->reset_amount = price::approximate(100'000.00, currencies::USD);
        }
        return result_;
    }

    ///
    /// \brief  This is to warn the modeler that the model is stuck on something
    ///         unexpected for the model:
    ///         prices not changing, or prices becoming negative
    ///
    /// \param early_exit_after If the model is stuck after this period, stop the model automatically
    /// \return
    bool is_stuck(size_t early_exit_after = 10000)
    {
        auto ps_ =  dynamic_pointer_cast<price_setter>(agents.local_agents_[identity<agent>{1}]);

        if(ps_->output_clearing_prices_->values.size() > early_exit_after){
            if(0. > double(get<1>(ps_->output_clearing_prices_->values.back())[0])){
                LOG(errorlog) << "negative prices!" << std::endl;
                return true;
            }
            auto j = ps_->output_clearing_prices_->values.rbegin();
            auto k = ps_->output_clearing_prices_->values.rbegin();

            bool equal_ = true;
            for(auto i = 0; i < early_exit_after; ++i, ++k){
                if(get<1>(*j)[0] != get<1>(*k)[0]){
                    equal_ = false;
                }
            }
            if(equal_){
                LOG(errorlog) << " stuck on same price!" << std::endl;
                return true;
            }
        }
        return false;
    }
};

///
/// \brief  Figure in paper with a single trace
/// \param stocks_count
/// \return
int volatility_illustration( unsigned int stocks_count )
{
    size_t sample = 1;

    environment e;
    main_model model_(e, parametrization(0, 0, 252 * 200));

    vector<tuple<shared_ptr<traded_company>, share_class>> shares_;
    map<tuple<identity<traded_company>, share_class>, identity<property>> stocks_;
    property_map<markets::quote> traded_assets_;
    property_map<size_t> shares_outstanding_;

    double nt = 0.42;
    double fv = 0.33;
    double tf = 0.25;

    std::string prefix_ = "output/volatility_illustration";

    for (size_t a = 0; a < stocks_count; ++a) {
        auto traded_company_ = model_.template create<traded_company>(model_.world, jurisdictions::US, sample);

        std::cout << (std::to_string(a) + "_dividend.txt", prefix_) << std::endl;

        traded_company_->outputs["dividend_per_share"]->streams.push_back(
            make_shared<data::file>(std::to_string(a) + "_dividend.txt", prefix_));

        auto main_issue_ = share_class();
        size_t main_issue_amount_ = 3 * 1'000'000u;
        traded_company_->shares_outstanding[main_issue_] = main_issue_amount_;

        for (const auto &[share_, quantity]: traded_company_->shares_outstanding) {
            auto stock_ = std::make_shared<stock>(*traded_company_, share_);

            shares_outstanding_[stock_] = main_issue_amount_;
            traded_assets_.emplace(stock_, markets::quote(price::approximate(100.00, USD)));
            shares_.emplace_back(make_tuple(traded_company_, share_));
            auto key_ = make_tuple<identity<traded_company>, share_class>(*traded_company_,
                                                                          share_class(share_));
            stocks_.insert({key_, stock_->identifier});
        }
    }

    auto market_ = model_.template create<price_setter>();
    market_->outputs["clearing_prices"]->buffered = false;
#if !defined(ESL_RELEASE) || (ESL_RELEASE == 0)
    market_->outputs["clearing_prices"]->streams.push_back(make_shared<data::terminal>(data::terminal::out));
#endif
    market_->outputs["clearing_prices"]->streams.push_back(make_shared<data::file>("prices.txt", prefix_));
    market_->outputs["volumes"]->streams.push_back(make_shared<data::file>("volumes.txt", prefix_));

    LOG(notice) << market_->describe() << std::endl;
    market_->traded_properties = traded_assets_;

    auto set_outputs = [&](shared_ptr<fund> f) {
        auto repr_ = f->identifier.representation();
        std::replace(repr_.begin(), repr_.end(), '"', '_');

        f->outputs["signal"]->streams.push_back(make_shared<data::file>(repr_ + "_signal.txt", prefix_));
        f->outputs["net_asset_value"]->streams.push_back(make_shared<data::file>(repr_ + "_net_asset_value.txt", prefix_));
        f->outputs["cash"]->streams.push_back(make_shared<data::file>(repr_ + "_cash.txt", prefix_));
        f->outputs["stocks"]->streams.push_back(make_shared<data::file>(repr_ + "_stocks.txt", prefix_));
        f->outputs["loans"]->streams.push_back(make_shared<data::file>(repr_ + "_loans.txt", prefix_));
        f->outputs["lending"]->streams.push_back(make_shared<data::file>(repr_ + "_lending.txt", prefix_));
        f->outputs["pnl"]->streams.push_back(make_shared<data::file>(repr_ + "_pnl.txt", prefix_));
    };

    auto target_nav_nt_ = price::approximate(nt * 200'000'000.00, USD);
    auto target_nav_fv_ = price::approximate(fv * 200'000'000.00, USD);
    auto target_nav_tf_ = price::approximate(tf * 200'000'000.00, USD);

    auto participants_ = vector<shared_ptr<fund>>();

    size_t target_date_     = 252;
    size_t noise_traders    = 1;
    size_t value_investors  = 1;
    size_t trend_followers  = 1;

    for(size_t i = 0; i < noise_traders; ++i) {
        auto mr_ = model_.template create<mean_reverting_noise_trader>(model_.world);
        mr_->target_net_asset_value.emplace(target_nav_nt_);
        mr_->target_date = target_date_;
        mr_->aggression = 5.0;
        mr_->maximum_leverage = 1.0;
        participants_.push_back(mr_);
        set_outputs(mr_);
    }

    for(size_t i = 0; i < value_investors; ++i) {
        auto fv_ = model_.template create<dividend_discount>(model_.world);
        fv_->target_net_asset_value.emplace(target_nav_fv_);
        fv_->target_date = target_date_;
        fv_->aggression = 6.0;
        fv_->maximum_leverage = 4.0;
        participants_.push_back(fv_);
        set_outputs(fv_);
    }

    for(size_t i = 0; i < trend_followers; ++i) {
        auto tf_ = model_.template create<momentum>(model_.world,esl::law::jurisdictions::US);
        tf_->target_net_asset_value.emplace(target_nav_tf_);
        tf_->target_date = target_date_;
        tf_->aggression = 1.5;
        tf_->maximum_leverage = 0.8;
        participants_.push_back(tf_);
        set_outputs(tf_);
    }

    for (const auto &fund: participants_) {
        LOG(trace) << fund->identifier << " " << fund->describe() << std::endl;
    }

    auto cash_ = make_shared<cash>(USD);

    //
    // Set cash for each agent
    //
    quantity cash_amounts_[noise_traders + value_investors + trend_followers];
    for(size_t i = 0; i < noise_traders; ++i) {
        cash_amounts_[i] = cash_->amount(nt > 0 ? 10'000'000'000.00 : 0.00);
    }

    for(size_t i = noise_traders; i < noise_traders + value_investors; ++i) {
        cash_amounts_[i] = cash_->amount(fv > 0 ? 10'000'000'000.00 : 0.00);
    }

    for(size_t i = noise_traders + value_investors; i < noise_traders + value_investors + trend_followers; ++i) {
        cash_amounts_[i] = cash_->amount(tf > 0 ? 10'000'000'000.00 : 0.00);
    }

    //
    //  Set the amount of loans for each agent
    //
    quantity loan_amounts[noise_traders + value_investors + trend_followers];
    for(size_t i = 0; i < noise_traders; ++i) {
        loan_amounts[i] = quantity(nt > 0 ? 9'900'000'000 + (1. - nt) * 100'000'000 : 0);
    }

    for(size_t i = noise_traders; i < noise_traders + value_investors; ++i) {
        loan_amounts[i] = quantity(fv > 0 ? 9'900'000'000 + (1. - fv) * 100'000'000 : 0);
    }

    for(size_t i = noise_traders + value_investors; i < noise_traders + value_investors + trend_followers; ++i) {
        loan_amounts[i] = quantity(tf > 0 ? 9'900'000'000 + (1. - tf) * 100'000'000 : 0);
    }

    //
    //  Set the starting amounts of stocks for the traders
    //
    quantity stock_amounts[noise_traders + value_investors + trend_followers];
    for(size_t i = 0; i < noise_traders; ++i) {
        stock_amounts[i] = quantity(nt > 0 ? nt * 1000'000 : 0);
    }
    for(size_t i = noise_traders; i < noise_traders + value_investors; ++i) {
        stock_amounts[i] = quantity(fv > 0 ? fv * 1000'000 : 0);
    }
    for(size_t i = noise_traders + value_investors; i < noise_traders + value_investors + trend_followers; ++i) {
        stock_amounts[i] = quantity(0);
    }

    for (auto[i, p] : enumerate(participants_)) {
        // 1. Add cash to participants
        (*p).shareholder::owner<cash>::take(cash_, cash_amounts_[i]);

        auto l_ = make_shared<loan>(p->identifier, p->identifier);
        (*p).shareholder::owner<loan>::take(l_, loan_amounts[i]);

        market_->participants.insert(*p);   // 2. Add stocks to participant

        for (const auto &[k, v] : stocks_) {
            p->stocks.insert(std::make_pair(k, v));
        }
        for (const auto &[k, v] : traded_assets_) {
            if (dynamic_pointer_cast<stock>(k)) {
                (*p).shareholder::owner<stock>::take(dynamic_pointer_cast<stock>(k), stock_amounts[i]);
            }
        }
        LOG(notice) << p->describe() << " " << p->identifier << " " << p->inventory << std::endl;
        for (auto &share: shares_) {
            get<0>(share)->shareholders[*p] = {{get<1>(share), stock_amounts[i].amount}};
        }
    }
    e.run(model_);

    return 0;
}

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

    if("experiment_3" == experiment_){
        experiment_3();
    }else if("experiment_3_best" == experiment_){
        experiment_3_best(64);
    }else if("experiment_6" == experiment_){
        experiment_6();
    }

    return 0;
}