//
// Created by Maarten on 03/05/2020.
//

#include "mean_reverting_noise.hpp"

using namespace esl;
using namespace esl::economics;
using namespace esl::simulation;
using namespace esl::economics::finance;
using esl::economics::markets::walras::quote_message;
using esl::simulation::time_interval;
using esl::identity;
using esl::law::property;


mean_reverting_noise_trader::mean_reverting_noise_trader(const identity<fund> &i, const jurisdiction &j, size_t seed)
: agent(i)
, owner<cash>(i)
, owner<stock>(i)
, fund(i, j)
, seed(seed)
{

}

time_point mean_reverting_noise_trader::invest(
        std::shared_ptr<quote_message> message,
        time_interval interval, std::seed_seq &seed)
{
    if(message->received < interval.lower){
        return interval.lower;
    }

    if(this->target_net_asset_value.has_value() && double(target_net_asset_value.value()) <= 1.){
        return interval.upper;
    }

    auto nav_ = net_asset_value(interval);
    if(nav_.value <= 0){
        return interval.upper;
    }

    //LOG(trace) << describe() << " " << identifier <<  " inventory " <<  inventory << std::endl;

    std::default_random_engine generator(seed );



    if(variates.empty()) {
        std::normal_distribution<double> distribution_(1.0, 0.0000001);
        std::map<identity<law::property>, double> variates_;
        for(auto [k, v]: message->proposed) {
            auto noise_ = distribution_(generator);
            LOG(errorlog) << "noise " << noise_ << std::endl;
            variates_.insert({k->identifier, noise_});
        }
        variates.insert(std::make_pair(interval.lower, variates_));
    } else {

        auto iterator_ = variates.rbegin();

        std::normal_distribution<double> distribution_(0.0, 1.0);

        double sigma_   = 0.2 * std::sqrt(1./252);
        double rho_     = 0.00045832561;//0.00045832561;
        double mu_      = 1.0;

        while(variates.size() > 1) {
            auto t_j = variates.begin();
            variates.erase(t_j);
        }

        if(iterator_->first != interval.lower) {
            std::map< identity< law::property>, double> variates_;
            for(auto [k, v] : message->proposed) {
                auto j  = iterator_->second.find(k->identifier);
                auto Yt = j->second;
                auto zi = distribution_(generator);
                auto dt = std::sqrt(interval.lower - iterator_->first);
                auto E1 = rho_ * (std::log2(mu_) - std::log2(Yt))   * dt;
                auto E2 = sigma_ * Yt * zi;

                if(iterator_->second.end() == j) {
                    throw std::logic_error("stochastic process was started with no history");
                }
                variates_.insert({k->identifier, Yt + E1 + E2 });
            }
            variates.insert(std::make_pair(interval.lower, variates_));
        }

    }

    std::map<identity<property>, price> valuations_;
    for(auto [property_, quote_]: message->proposed) {

        auto price_ = std::get<price>(quote_.type);
        //auto [dd, b_] = drawdown.emplace(property_->identifier, std::map<simulation::time_point, price>());
        //dd->second.emplace(interval.lower, price_);
        auto i = market_data.dividend_per_share.find(property_->identifier);
        // it is assumed all assets are quoted with a price
        auto val_ = std::get<price>(quote_.type);
        if(i != market_data.dividend_per_share.end()){
            // in this simplified model, we assume there is one dividend payment
            // per period, and therefore the dividend yield is in USD/day
            assert(i->second.valuation == price_.valuation);

            auto payment_ = double(i->second);
            auto shares_outstanding_ = market_data.shares_outstanding.find(property_->identifier)->second;
            auto dividend_rate_ = (payment_ / shares_outstanding_);

            auto growth_ = 0.000'1;
            auto lower_limit_ = 0.000'000'1;
            growth_ = std::max(growth_, lower_limit_);

            auto compounded_rate_return_ = 0.000'2; // FVI expects company to appreciate at RFR

            auto gordon_ = dividend_rate_ / (std::max(0.000'1, compounded_rate_return_ - growth_));
            gordon_ = std::max(std::min(gordon_, 1'000.00), lower_limit_);

            //LOG(trace) << "Gordon growth model " << gordon_<< "(mu=" << std::fixed << std::setprecision(12)
            //           << growth_ << ",d_1=" << std::setprecision(5)<< dividend_rate_ << ")" <<  std::endl;

            val_.value =  int64_t(gordon_ * variates[interval.lower][property_->identifier] * 100) ;
            //output_signal->put(interval.lower,gordon_ *  variates[interval.lower][property_->identifier]);
        }

        if(valuations_.end() != valuations_.find(property_->identifier)){
            valuations_.find(property_->identifier)->second.value = val_.value;
        }else{
            valuations_.emplace(property_->identifier, val_);
        }
    }

    //LOG(trace) << "create_message<dividend_discount_ddsf> with "  << valuations_.size() << " valuations: " << valuations_ << std::endl;
    auto message_ = this->template create_message<mean_reverting_noise_ddsf>(
            message->sender, interval.lower, (*this), message->sender,
            interval.lower, interval.lower, nav_, valuations_);


    message_->agression = this->aggression;
    message_->leverage = this->maximum_leverage;

    for(auto [p,q]: owner<stock>::properties.items){
        if(0 == q.amount){
            continue;
        }
        message_->supply.emplace(p->identifier, std::make_tuple(q, quantity(0)));
    }
    for(auto [p,q]: owner<securities_lending_contract>::properties.items){
        if(0 == q.amount){
            continue;
        }
        if(message_->supply.end() != message_->supply.find(p->security)){
            std::get<1>( message_->supply.find(p->security)->second ) = q;
        }else{
            message_->supply.emplace(p->security, std::make_tuple(quantity(0), q));
        }
    }

    return interval.upper;
}

///
///
///
mean_reverting_noise_ddsf::mean_reverting_noise_ddsf
        ( const esl::identity<esl::agent> &sender
                , const esl::identity<esl::agent> &recipient
                , time_point sent
                , time_point received
                , price nav
                , std::map<identity<law::property>, economics::price> valuations
        )
        : differentiable_order_message(sender, recipient, sent, received)
        , valuations(/*std::move*/(valuations))
        , net_asset_value(nav)
{

}

std::map<identity<law::property>, variable>
mean_reverting_noise_ddsf::excess_demand(
        const std::map<identity<law::property>,
                std::tuple<markets::quote, esl::variable>> &quotes) const
{
    std::map<esl::identity<property>, esl::variable> result_;
    auto scale_ = double(net_asset_value) / quotes.size();
    //LOG(warning) << "net asset value" << scale_ << std::endl;
    for(auto &[k, v] : quotes) {
        const auto &[quote_, variable_] = v;
        const auto quoted_price_ = double(std::get<price>(quote_.type));
        auto i = valuations.find(k);
        if(valuations.end() != i) {
            auto value_ = double(i->second);
            auto j = this->supply.find(k);
            if(supply.end() == j){
                result_.emplace(k,  scale_ / (quoted_price_ * variable_) * (value_ -  (quoted_price_ * variable_)) - 0);
            }else{
                auto supply_long_ = double(std::get<0>(j->second));
                auto supply_short_ = double(std::get<1>(j->second));
                auto lambda_ = leverage * 2;
                result_.emplace(k, scale_ *
                                   (lambda_ / (1. + adept::exp( - agression * (std::log2(value_) - adept::log2(quoted_price_ * variable_)))) - lambda_/2 + 0.5)
                                   - (supply_long_ - supply_short_) * (quoted_price_ * variable_)
                );
            }
        }
    }
    //LOG(trace) << "excess demand " << result_ << std::endl;
    return result_;
}