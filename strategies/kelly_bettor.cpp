/// \file   kelly_bettor.cpp
///
/// \brief
///
/// \authors    maarten
/// \date       2020-11-06
/// \copyright  Copyright 2017-2020 The Institute for New Economic Thinking,
///             Oxford Martin School, University of Oxford
///
///             Licensed under the Apache License, Version 2.0 (the "License");
///             you may not use this file except in compliance with the License.
///             You may obtain a copy of the License at
///
///                 http://www.apache.org/licenses/LICENSE-2.0
///
///             Unless required by applicable law or agreed to in writing,
///             software distributed under the License is distributed on an "AS
///             IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
///             express or implied. See the License for the specific language
///             governing permissions and limitations under the License.
///
///             You may obtain instructions to fulfill the attribution
///             requirements in CITATION.cff
///
#include "kelly_bettor.hpp"

#include <utility>

#include <boost/accumulators/statistics/mean.hpp>

#include <esl/data/log.hpp>
using economics::finance::securities_lending_contract;

using namespace esl;
using namespace esl::economics;
using namespace esl::economics::markets;


kelly_bettor::kelly_bettor(const identity<fund> &i, const jurisdiction &j)
: agent(i)
, owner<cash>(i)
, owner<stock>(i)
, fund(i, j)
, past_price(100'00, currencies::USD)
{
    estimates(0.);
}

time_point kelly_bettor::invest(std::shared_ptr<walras::quote_message> message, simulation::time_interval interval, std::seed_seq &seed)
{
    if(message->received < interval.lower){
        return interval.lower;
    }

    demand.clear();
    for(auto [p,q]: message->proposed){
        demand.emplace(p->identifier, 1.0);
    }

    auto nav_ =  net_asset_value(interval);
    LOG(trace) << describe() << " " << identifier << " inventory " <<  inventory << std::endl;
    double update_value_ = 0.;

    if(interval.lower <= 100){
        return interval.lower;
    }

    auto m = this->template create_message<kelly_bettor_ddsf>(
        message->sender, interval.lower, (*this), message->sender,
        interval.lower, interval.lower,nav_,demand);


    price price_ = price::approximate(100.00, currencies::USD);
    double dividend_rate = 0.;
    for(auto [property_, quote_] : message->proposed) {

        price_ = std::get<price>(quote_.type);

        auto i = market_data.dividend_per_share.find(*property_);
        // it is assumed all assets are quoted with a price
        if(i != market_data.dividend_per_share.end()) {
            // in this simplified model, we assume there is one dividend payment per period, and therefore the dividend yield is in USD/day

            auto payment_ = double(i->second);
            auto shares_outstanding_ = market_data.shares_outstanding.find(property_->identifier)->second;
            dividend_rate = (payment_ / shares_outstanding_) ;
            update_value_ = 0.01 * std::min(interval.lower, 100ul) * (
                                (double(price_) / double(past_price)) - 1 + dividend_rate/ double(price_));
            past_price = price_;
        }
    }
    m->estimates = estimates;
    m->time = interval.lower;

    for(auto [p,q]: inventory){
        auto cast_ = std::dynamic_pointer_cast<stock>(p);
        if(cast_){
            if(0 == q.amount){
                continue;
            }
            m->supply.emplace(p->identifier, std::make_tuple(q, quantity(0)));
        }else{
            auto cast2_ = std::dynamic_pointer_cast<securities_lending_contract>(p);
            if(cast2_){
                if(0 == q.amount){
                    continue;
                }
                if(m->supply.end() != m->supply.find(cast2_->security)){
                    std::get<1>( m->supply.find(cast2_->security)->second ) = q;
                }else{
                    m->supply.emplace(cast2_->security, std::make_tuple(quantity(0), q));
                }
            }
        }
    }

    for(auto [p,q]: owner<stock>::properties.items){


    }
    for(auto [p,q]: owner<securities_lending_contract>::properties.items){

    }



    m->past_price = past_price;
    m->dividend_rate = dividend_rate;
    estimates(update_value_);
    return interval.lower;
}

std::string kelly_bettor::describe() const
{
    return "kelly bettor";
}

///
///
///
kelly_bettor_ddsf::kelly_bettor_ddsf(
    const identity<agent> &sender,
    const identity<agent> &recipient,
    simulation::time_point sent, simulation::time_point received,
    price net_asset_value,
    std::map<identity<law::property>, double> variates
)
    : differentiable_order_message(sender, recipient, sent, received)
    , variates(std::move(variates))
    , net_asset_value(net_asset_value)
{

}

std::map<identity<law::property>, esl::variable>
kelly_bettor_ddsf::excess_demand(
    const std::map<identity<law::property>,
        std::tuple<economics::markets::quote, variable>> &quotes)
const
{
    std::map<esl::identity<property>, esl::variable> result_;
    auto scale_ = double(net_asset_value) / quotes.size();

    for(auto &[k, v] : quotes) {
        const auto &[quote_, variable_] = v;
        const auto quoted_price_ = double(std::get<price>(quote_.type));
        auto i = variates.find(k);
        if(variates.end() != i) {
            auto value_ = i->second;
            auto j = supply.find(k);
            if(supply.end() == j){
                result_.emplace(k,  scale_ * value_ * value_);
            }else{
                auto supply_long_  = double(std::get<0>(j->second));
                auto supply_short_ = double(std::get<1>(j->second));

                double ease_in_ = 0.01 * std::min(time, 100ul);
                ease_in_ *= ease_in_;

                auto m = boost::accumulators::mean(estimates);
                m = std::max(-0.5, std::min(0.5, std::pow(1. + ease_in_ * m, 252.)-1));

                double x = quoted_price_ * adept::value(variable_);
                auto update_value_ = ease_in_ * ((x / double(past_price)) - 1 + dividend_rate/ x);


                const_cast<decltype(estimates) &>( estimates)(update_value_);

                auto sigma = std::sqrt(boost::accumulators::variance(estimates));
                sigma = std::max(0.001, sigma);
                sigma *= std::sqrt(252);

                // fractional kelly
                double c = ease_in_;
                double r = 0.01;  //   annual
//
//                std::cout << std::setprecision(6) << "price=" << x
//                          << " signal " << signal
//                          << " mu " << mu
//                          << " sigma " << sigma
//                          << " d " << dividend_rate
//                          << " t " << time << std::endl;

                auto lambda_ = 0.9 * 2;
                result_.emplace(k, scale_ * (lambda_ / (1. + adept::exp( - (c * ((m*time + (ease_in_ * ((quoted_price_ * variable_) / double(past_price) - 1 + dividend_rate/ (quoted_price_ * variable_))))/(time+1) - r) / (sigma*sigma)) )) - lambda_/2 + 0.5)
                                       - (supply_long_ - supply_short_) * (quoted_price_ * variable_)
                );
            }
        }
    }
    return result_;
}
