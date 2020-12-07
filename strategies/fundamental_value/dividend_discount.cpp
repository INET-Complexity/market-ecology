/// \file   dividend_discount.cpp
///
/// \brief
///
/// \authors    Maarten P. Scholl
/// \date       2020-01-30
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
#include "dividend_discount.hpp"

using namespace esl;
using namespace esl::economics;
using namespace esl::simulation;
using namespace esl::economics::finance;
using esl::economics::markets::walras::quote_message;
using esl::simulation::time_interval;
using esl::identity;
using esl::law::property;


dividend_discount::dividend_discount(const identity<fund> &i, const jurisdiction &j)
: agent(i)
, owner<cash>(i)
, owner<stock>(i)
, fund(i, j)
{

}

time_point dividend_discount::invest(
    std::shared_ptr<quote_message> message,
    time_interval interval, std::seed_seq &seed)
{
    if(message->received + 1 < interval.lower){
        return interval.lower;
    }
    auto nav_ = net_asset_value(interval);
    if(nav_.value <= 0){
        return interval.upper;
    }


    if(target_net_asset_value.has_value() && double(target_net_asset_value.value()) <= 1.){
        return interval.upper;
    }

    std::map<identity<property>, price> valuations_;
    for(auto [property_, quote_]: message->proposed) {
        auto price_ = std::get<price>(quote_.type);

        auto i = market_data.dividend_per_share.find(*property_);
        // it is assumed all assets are quoted with a price
        auto fundamental_value_ = std::get<price>(quote_.type);
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
            auto compounded_rate_return_ = 0.000'2; // FVI expects company to appreciate at this rate
            auto gordon_ = dividend_rate_ / (std::max(0.000'1, compounded_rate_return_ - growth_));
            gordon_ = std::max(std::min(gordon_, 1'000.00), lower_limit_);
            fundamental_value_.value =  int64_t(gordon_* 100) ;
            output_signal->put(interval.lower, gordon_);
        }

        if(valuations_.end() != valuations_.find(property_->identifier)){
            valuations_.find(property_->identifier)->second.value = fundamental_value_.value;
        }else{
            valuations_.emplace(property_->identifier, fundamental_value_);
        }
    }

    auto message_ = this->template create_message<dividend_discount_ddsf>(
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
dividend_discount_ddsf::dividend_discount_ddsf
( const identity<agent> &sender
, const identity<agent> &recipient
, time_point sent
, time_point received
, price nav
, std::map<identity<law::property>, price> valuations
)
: differentiable_order_message(sender, recipient, sent, received)
, valuations(valuations)
, net_asset_value(nav)
{

}

std::map<identity<law::property>, variable>
dividend_discount_ddsf::excess_demand(
    const std::map<identity<law::property>,
    std::tuple<markets::quote, variable>> &quotes) const
{
    std::map<identity<property>, variable> result_;
    auto scale_ = double(net_asset_value) / quotes.size();
    for(auto &[k, v] : quotes) {
        const auto &[quote_, variable_] = v;
        const auto quoted_price_ = double(std::get<price>(quote_.type));
        auto i = valuations.find(k);
        if(valuations.end() != i) {
            auto value_ = double(i->second);
            auto j = supply.find(k);

            // Scale by two,
            auto lambda_ = leverage * 2;
            // if we have a neutral position
            if(supply.end() == j){
                result_.emplace(k,  scale_ *(lambda_ / (1. + adept::exp(- agression * (std::log2(value_) - adept::log2(quoted_price_ * variable_)))) - lambda_/2 +0.5  )
                                    * (quoted_price_ * variable_));
            }else{
                auto supply_long_ = double(std::get<0>(j->second));
                auto supply_short_ = double(std::get<1>(j->second));
                result_.emplace(k, scale_ *(lambda_ / (1. + adept::exp(- agression * (std::log2(value_) - adept::log2(quoted_price_ * variable_)))) - lambda_/2 +0.5  )
                - (supply_long_ - supply_short_) * (quoted_price_ * variable_));
            }
        }
    }
    //LOG(trace) << "excess demand " << result_ << std::endl;
    return result_;
}