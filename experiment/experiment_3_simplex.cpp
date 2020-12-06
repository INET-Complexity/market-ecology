//
// Created by Maarten on 11/05/2020.
//

#include "experiment_3_simplex.hpp"

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

#include "../traded_company.hpp"

using esl::simulation::time_duration;
using esl::economics::rate;
using esl::computation::environment;
using esl::simulation::model;
#include <esl/data/log.hpp>
#include <esl/data/file.hpp>
#include <esl/algorithms.hpp>

////////////////////////////////////////////////////////////////////////////////


using std::to_string;
using std::tuple;
using std::string;
using std::pair;
using std::vector;
using std::map;
using std::make_shared;
using std::shared_ptr;
using std::make_tuple;


////////////////////////////////////////////////////////////////////////////////

#include "../strategies/fundamental_value/mean_reverting_noise.hpp"
#include "../strategies/fundamental_value/dividend_discount.hpp"
#include "../strategies/technical/momentum.hpp"
#include "../strategies/constant_demand.hpp"

////////////////////////////////////////////////////////////////////////////////

#include <boost/program_options.hpp>
using namespace boost::program_options;

int experiment_3_task(std::uint64_t sample, std::uint64_t assets, double nt, double fv, double tf

                      , double nt_agg, double nt_lev
                      , double fv_agg, double fv_lev
                      , double tf_agg, double tf_lev)
{
    unsigned int stocks_count               = assets;
    environment e;
    model model_(e, parametrization(0, 0, 11*252));

    std::dynamic_pointer_cast<simulation::parameter::constant<std::uint64_t>>(model_.parameters.values["sample"])->choice = sample;

    vector<tuple<shared_ptr<traded_company>, share_class>> shares_;
    map<tuple<identity<traded_company>, share_class>, identity<property>> stocks_;
    property_map<markets::quote> traded_assets_;
    property_map<size_t> shares_outstanding_;

    const std::string prefix_ = std::string("output/experiment3_alternative1_")
                                + std::to_string(nt_agg) + "_"
                                + std::to_string(nt_lev) + "_"
                                + std::to_string(fv_agg) + "_"
                                + std::to_string(fv_lev) + "_"
                                + std::to_string(tf_agg) + "_"
                                + std::to_string(tf_lev)
                                + "/run_"
                                + std::to_string(nt) + "_"
                                + std::to_string(fv) + "_"
                                + std::to_string(tf) + "_sample_"
                                + std::to_string(sample) +"/" ;

    if(std::filesystem::exists(prefix_)){
        return -1;
    }

    std::cout << "starting experiment: " << prefix_ << std::endl;

    std::filesystem::create_directories(prefix_);

    for(size_t a = 0; a < stocks_count; ++a){
        auto traded_company_ = model_.template create<traded_company>(model_.world, jurisdictions::US);

        traded_company_->outputs["dividend_per_share"]->streams.push_back(make_shared<data::file>(std::to_string(a) + "_dividend.txt", prefix_));

        auto main_issue_ = share_class();
        size_t main_issue_amount_ = 3 * 1'0'000u;
        traded_company_->shares_outstanding[main_issue_] = main_issue_amount_;

        for(const auto &[share_, quantity]: traded_company_->shares_outstanding){
            auto stock_ = std::make_shared<stock>(*traded_company_, share_);

            shares_outstanding_[stock_] = main_issue_amount_;
            traded_assets_.emplace(stock_, markets::quote(price::approximate(100.0, USD)));
            shares_.emplace_back(make_tuple(traded_company_, share_));
            auto key_ = make_tuple<identity<traded_company>, share_class>(*traded_company_, share_class(share_));
            stocks_.insert({key_, stock_->identifier});
        }
    }

    auto market_ = model_.template create<price_setter>();
    market_->outputs["clearing_prices"]->buffered = false;
#if !defined(ESL_RELEASE) || (ESL_RELEASE == 0)
    market_->outputs["clearing_prices"]->streams.push_back(make_shared<data::terminal>(data::terminal::out));
#endif
    market_->outputs["clearing_prices"]->streams.push_back(make_shared<data::file>("prices.txt", prefix_));

    LOG(notice) << market_->describe() << std::endl;
    market_->traded_properties = traded_assets_;
    auto set_outputs = [&](shared_ptr<fund> f){
        auto repr_ = f->identifier.representation();
        std::replace(repr_.begin(), repr_.end(), '"', '_');

        //f->outputs["signal"]->streams.push_back(make_shared<data::file>(repr_ + "_signal.txt", prefix_));
        f->outputs["net_asset_value"]->streams.push_back(make_shared<data::file>(repr_ + "_net_asset_value.txt", prefix_));
        f->outputs["pnl"]->streams.push_back(make_shared<data::file>(repr_ + "_pnl.txt", prefix_));
    };

    auto target_nav_nt_ = price::approximate(nt * 2'000'000.00, USD);
    auto target_nav_fv_ = price::approximate(fv * 2'000'000.00, USD);
    auto target_nav_tf_ = price::approximate(tf * 2'000'000.00, USD);

    auto participants_ = vector<shared_ptr<fund>>();

    auto mr_ = model_.template create<mean_reverting_noise_trader>();
    mr_->seed = sample;
    mr_->target_net_asset_value.emplace(target_nav_nt_);
    mr_->aggression = nt_agg;
    mr_->maximum_leverage = nt_lev;
    participants_.push_back(mr_);
    set_outputs(mr_);

    auto fv_ = model_.template create<dividend_discount>();
    fv_->target_net_asset_value.emplace(target_nav_fv_);
    fv_->aggression = fv_agg;
    fv_->maximum_leverage = fv_lev;
    participants_.push_back(fv_);
    set_outputs(fv_);

    auto tf_ = model_.template create<momentum>();
    tf_->target_net_asset_value.emplace(target_nav_tf_);
    tf_->aggression = tf_agg;
    tf_->maximum_leverage = tf_lev;
    participants_.push_back(tf_);
    set_outputs(tf_);

    for(const auto &fund: participants_){
        LOG(trace) << fund->identifier << " " << fund->describe() << std::endl;
    }
    auto cash_ = std::make_shared<cash>(USD);
    quantity cash_amounts_[3] =  {
            cash_->amount(nt > 0 ?  1000'000'000.00 : 0.00),
            cash_->amount(fv > 0 ?  1000'000'000.00 : 0.00),
            cash_->amount(tf > 0 ?  1000'000'000.00 : 0.00),
    };
    quantity loan_amounts[3] =  {
            quantity(nt > 0 ? 999'000'000 + (1. - nt) * 1'000'000 : 0),
            quantity(fv > 0 ? 999'000'000 + (1. - fv) * 1'000'000 : 0),
            quantity(tf > 0 ? 999'000'000 + (1. - tf) * 1'000'000 : 0),
    };
    quantity stock_amounts[3] =  {
            quantity(nt > 0 ? nt * 10'000 : 0),
            quantity(fv > 0 ? fv * (nt == 0? 10'000 : 0) : 0),
            quantity(tf > 0 ? tf * (nt == 0 && fv == 0? 10'000 : 0) : 0),
    };

    for(auto [i, p] : enumerate(participants_)){           // 1. Add cash to participants
        (*p).shareholder::owner<cash>::take(cash_, cash_amounts_[i]);
        auto l_ = std::make_shared<loan>(p->identifier, p->identifier);
        (*p).shareholder::owner<loan>::take(l_, loan_amounts[i]);

        market_->participants.insert(*p);   // 2. Add stocks to participant

        for(const auto &[k, v] : stocks_){
            p->stocks.insert(std::make_pair(k, v));
        }
        for(const auto &[k, v] : traded_assets_){
            if(std::dynamic_pointer_cast<stock>(k)){
                (*p).shareholder::owner<stock>::take( std::dynamic_pointer_cast<stock>(k), stock_amounts[i]);
            }
        }
        LOG(notice) << p->describe() << " " << p->identifier << " " << p->inventory << std::endl;
        for(auto &share: shares_) {
            std::get<0>(share)->shareholders[*p] = {{std::get<1>(share), stock_amounts[i].amount}};
        }
    }

    try {
        e.run(model_);
    }catch(...){
        std::cout << "ERROR " << prefix_ << std::endl;
    }

    return 0;
}


///
/// \param sample
/// \return
int experiment_3_(uint64_t sample, double nt_agg, double nt_lev, double fv_agg, double fv_lev, double tf_agg, double tf_lev
                 , size_t subsamples = 64 )
{
    uint64_t subsamples_ = subsamples; // 64
    esl::computation::thread_pool pool_;
    std::vector<std::vector<double>> combinations_ = esl::sample_barycentric(3, subsamples_);
    LOG(warning) << combinations_ << std::endl;
    std::cout << "starting " << combinations_.size() << " tasks" << std::endl;
    std::vector<std::future<int>> results_;
    for (auto c: combinations_) {
        results_.emplace_back(pool_.enqueue_task(experiment_3_task
                                                 , sample
                                                 , 1
                                                 , c[0]
                                                 , c[1]
                                                 , c[2]
                                                 , nt_agg
                                                 , nt_lev
                                                 , fv_agg
                                                 , fv_lev
                                                 , tf_agg
                                                 , tf_lev) );
        if(results_.size() >= 8){
            for(auto &r: results_){
                r.wait();
            }
            results_.clear();
        }
    }
    return 0;
}

///
/// \return
int experiment_3()
{
    for(double nt_agg: {5.0, 2.0, 10.0, 1.0}){
        for(double nt_lev: {1.0, 2.0, 4.0}){
            for(double fv_agg: {10., 12., 8.0, 4.0, 2.0, 1.0, 16.0}){
                for(double fv_lev: {2.0, 4.0, 1.0,  8.0}){
                    for(double tf_agg: {1.0, 0.5, 2.0, 3.0, 4.0}){
                        for(double tf_lev: {1.0, 0.8,  1.5, 2.0}){
                            for(size_t s = 0; s < 2; ++s){
                                experiment_3_(s , nt_agg
                                    , nt_lev
                                    , fv_agg
                                    , fv_lev
                                    , tf_agg
                                    ,  tf_lev);
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}



int experiment_3_best(unsigned int precision)
{
    for(size_t s = 0; s < 8; ++s){
        experiment_3_(s
            , 5.
            , 1.
            , 10.
            , 8.
            , 4.
            , 0.8
            , precision);
    }
    return 0;
}