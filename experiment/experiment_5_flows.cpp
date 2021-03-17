//
// Created by Maarten on 02/06/2020.
//

#include "experiment_5_flows.hpp"


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

#include "../strategies/constant_demand.hpp"
#include "../strategies/fundamental_value/dividend_discount.hpp"
#include "../strategies/fundamental_value/mean_reverting_noise.hpp"
#include "../strategies/technical/trend_follower.hpp"

////////////////////////////////////////////////////////////////////////////////

struct me_model5
    : public model
{
    using model::model;

    constexpr static time_point end_time = 99999999999;

    time_point step(time_interval step) override
    {
        auto ps_ =  std::dynamic_pointer_cast<price_setter>(agents.local_agents_[identity<agent>{1}]);

        auto nt = std::dynamic_pointer_cast<fund>(agents.local_agents_[identity<agent>{2}]);
        auto fv = std::dynamic_pointer_cast<fund>(agents.local_agents_[identity<agent>{3}]);
        auto tf = std::dynamic_pointer_cast<fund>(agents.local_agents_[identity<agent>{4}]);

        if(ps_->output_clearing_prices_->values.size() > 100){

            // the model is not robust against negative prices..
            if(0. > double(std::get<1>(ps_->output_clearing_prices_->values.back())[0])){
                return end_time;
            }


            auto j = ps_->output_clearing_prices_->values.rbegin();
            auto k = ps_->output_clearing_prices_->values.rbegin();
            bool equal_ = true;
            for(auto i = 0; i < 25; ++i, ++k){
                if(std::get<1>(*j)[0] != std::get<1>(*k)[0]){
                    equal_ = false;
                }
            }

            if(equal_){
                return end_time;
            }
        }

        auto result_ = model::step(step);

        if(!std::dynamic_pointer_cast<fund>(agents.local_agents_[identity<agent>{2}])->output_net_asset_value->values.empty()){
            auto nt_ = double(std::get<1>(nt->output_net_asset_value->values.back()));
            auto fv_ = double(std::get<1>(fv->output_net_asset_value->values.back()));
            auto tf_ = double(std::get<1>(tf->output_net_asset_value->values.back()));

            auto total_ = nt_ + fv_ + tf_;

            //std::cout << std::setprecision(3) <<  nt_ << " ";
            //std::cout << std::setprecision(3) <<  fv_ << " ";
            //std::cout << std::setprecision(3) <<  tf_ << " ";
            //std::cout << std::endl;

            auto reset_amount_ = total_ * 0.001;

            double nw[3] = {nt_, fv_, tf_};
            std::shared_ptr<fund> agents_[3] = {nt, fv, tf};

            for(size_t i = 0; i < 3; ++i){
                if(nw[i] / total_ > 0.99){
                    return end_time;
                    for(auto &[p, q]: agents_[i]->inventory) {
                        auto cast_ = std::dynamic_pointer_cast<cash>(p);
                        if(cast_) {
                            agents_[i]->inventory[p].amount -= total_ * 0.01 * 100;
                        }
                    }
                }

                if(nw[i] / total_ < 0.01) {
                    return end_time;
                    //std::cout << i << " falls below 1%" << std::endl;
                    for(auto &[p, q] :  agents_[i]->inventory) {
                        auto cast_ = std::dynamic_pointer_cast<cash>(p);
                        if(cast_) {
                            agents_[i]->inventory[p].amount +=
                                std::uint64_t(total_ * 0.01 * 100);
                            //(nt->inventory[p].amount * 0.1);
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

};

////////////////////////////////////////////////////////////////////////////////

#include <boost/program_options.hpp>
using namespace boost::program_options;

std::vector<double> experiment_5_task(double ntv, double dpv, std::uint64_t sample, std::uint64_t assets, double nt, double fv, double tf, double reinvestment_rate = 1.0)
{
    unsigned int stocks_count               = assets;
    environment e;
    me_model5 model_(e, parametrization(0, 0, 252 * 100));

    std::dynamic_pointer_cast<simulation::parameter::constant<std::uint64_t>>(model_.parameters.values["sample"])->choice = sample;
    vector<tuple<shared_ptr<traded_company>, share_class>> shares_;
    map<tuple<identity<traded_company>, share_class>, identity<property>> stocks_;
    property_map<markets::quote> traded_assets_;
    property_map<size_t> shares_outstanding_;

    const std::string prefix_ = "output/experiment5_" + std::to_string(ntv) + "_"
                                + std::to_string(dpv) + "_"
                                + std::to_string(reinvestment_rate) + "/run_"
                                + std::to_string(nt) + "_"
                                + std::to_string(fv) + "_"
                                + std::to_string(tf) + "_"
                                  + "_sample_" + std::to_string(sample) +"/";

    if(std::filesystem::exists(prefix_)){
        return {nt, fv, tf};
    }
    std::filesystem::create_directories(prefix_);
    for(size_t a = 0; a < stocks_count; ++a){
        auto traded_company_ = model_.template create<traded_company>(model_.world, jurisdictions::US, sample);
        traded_company_->sigma = dpv / sqrt(double(traded_company_->day_count));
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

        f->outputs["net_asset_value"]->streams.push_back(make_shared<data::file>(repr_ + "_net_asset_value.txt", prefix_));
        f->outputs["pnl"]->streams.push_back(make_shared<data::file>(repr_ + "_pnl.txt", prefix_));
    };

    auto target_nav_nt_ = price::approximate(nt * 2'000'000.00, USD);
    auto target_nav_fv_ = price::approximate(fv * 2'000'000.00, USD);
    auto target_nav_tf_ = price::approximate(tf * 2'000'000.00, USD);

    auto participants_ = vector<shared_ptr<fund>>();

    double nt_agg = 5.;
    double nt_lev = 1.;
    double fv_agg = 10.;
    double fv_lev = 8.;
    double tf_agg = 4.;
    double tf_lev = 0.8;

    auto mr_ = model_.template create<mean_reverting_noise_trader>();
    mr_->seed = sample;
    mr_->target_net_asset_value.emplace(target_nav_nt_);
    mr_->aggression = nt_agg;
    mr_->maximum_leverage = nt_lev;
    mr_->target_date = 252;
    mr_->reinvestment_rate = reinvestment_rate;
    mr_->sigma_ = ntv * std::sqrt(1./252);
    participants_.push_back(mr_);
    set_outputs(mr_);

    auto fv_ = model_.template create<dividend_discount>();
    fv_->target_net_asset_value.emplace(target_nav_fv_);
    fv_->aggression = fv_agg;
    fv_->maximum_leverage = fv_lev;
    fv_->target_date = 252;
    fv_->reinvestment_rate = reinvestment_rate;
    participants_.push_back(fv_);
    set_outputs(fv_);

    auto tf_ = model_.template create<trend_follower>();
    tf_->target_net_asset_value.emplace(target_nav_tf_);
    tf_->aggression = tf_agg;
    tf_->maximum_leverage = tf_lev;
    tf_->target_date = 252;
    tf_->reinvestment_rate = reinvestment_rate;
    participants_.push_back(tf_);
    set_outputs(tf_);

    for(const auto &fund: participants_){
        LOG(trace) << fund->identifier << " " << fund->describe() << std::endl;
    }

    double large_precision_ = 1'000'000'000.00;

    auto cash_ = std::make_shared<cash>(USD);
    quantity cash_amounts_[3] =  {
        cash_->amount(nt > 0 ?  10'000'000'000.00 : 0.00),
        cash_->amount(fv > 0 ?  10'000'000'000.00 : 0.00),
        cash_->amount(tf > 0 ?  10'000'000'000.00 : 0.00),
    };
    quantity loan_amounts[3] =  {
        quantity(nt > 0 ? 9'999'000'000 + (1. - nt) * 1'000'000 : 0),
        quantity(fv > 0 ? 9'999'000'000 + (1. - fv) * 1'000'000 : 0),
        quantity(tf > 0 ? 9'999'000'000 + (1. - tf) * 1'000'000 : 0),
    };
    quantity stock_amounts[3] =  {
        quantity(nt > 0 ? nt * 10'000 : 0),
        quantity(fv > 0 ? fv * 10'000 : 0),
        quantity(tf > 0 ? tf * 10'000 : 0),
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
    }catch(std::exception &e){
        std::cout << "ERROR " << prefix_ << "  message: " << e.what() << std::endl;
    }catch(...){
        std::cout << "ERROR " << prefix_ << std::endl;
    }

    std::vector<double> result_;
    result_.push_back(double(mr_->net_asset_value({model_.time, model_.end})));
    result_.push_back(double(fv_->net_asset_value({model_.time, model_.end})));
    result_.push_back(double(tf_->net_asset_value({model_.time, model_.end})));

    auto norm_ = result_[0] + result_[1] + result_[2];

    result_[0] /= norm_;
    result_[1] /= norm_;
    result_[2] /= norm_;

    return result_;
}


///
/// \param sample
/// \return
int experiment_5(double ntv, double dpv, double reinvestment_rate, std::size_t threads)
{
    size_t progress_ = 0;
    uint64_t subsamples_ = 16;
    esl::computation::thread_pool pool_;
    std::vector<std::future<std::vector<double>>> results_;
    std::vector<std::vector<double>> combinations2_;

    for(size_t s = 0; s < 64; ++s) {
        std::vector<std::vector<double>> combinations_;
        for(auto c:  esl::sample_barycentric(3, subsamples_)){
            if (c[0] <= 0. || c[1] <= 0. || c[2] <= 0.) {
                continue;
            }
            combinations_.push_back(c);
        }

        sort( combinations_.begin()
            , combinations_.end()
            , [](const auto & a, const auto & b) -> bool
            {
                auto ga = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
                auto gb = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
                return ga < gb;
            });

        LOG(warning) << combinations_ << std::endl;
        std::cout << "starting " << combinations_.size() << " tasks" << std::endl;

        for (auto c: combinations_) {
            results_.emplace_back(pool_.enqueue_task(experiment_5_task, ntv, dpv, progress_, 1, c[0], c[1], c[2], reinvestment_rate));
            progress_ += 1;

            if (results_.size() >= std::thread::hardware_concurrency()){
                combinations2_.clear();
                for (auto &r: results_) {
                    r.wait();
                    combinations2_.push_back(r.get());
                }

                results_.clear();
                //combinations_ = combinations2_;
            }
        }
    }

    for (auto &r: results_) {
        r.wait();
        combinations2_.push_back(r.get());
    }
    return 0;
}
