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




    std::map<identity<property>, price> valuations_;
    for(auto [property_, quote_]: message->proposed) {
        auto price_ = std::get<price>(quote_.type);

        auto i = market_data.dividend_per_share.find(*property_);
        // it is assumed all assets are quoted with a price
        auto fundamental_value_ = std::get<price>(quote_.type);
        if(i != market_data.dividend_per_share.end()){
            // in this simplified model, we assume there is one dividend payment
            // per period, and therefore the dividend yield is in USD/day


            estimates(std::log(double(price_)) -  std::log(double(past_price)));

            auto m = std::exp(boost::accumulators::mean(estimates)) - 1;
            auto s = std::exp(boost::accumulators::moment<2>(estimates)) - 1;

            double sig = m/s;

            output_signal->put(interval.lower, sig);
        }

        if(valuations_.end() != valuations_.find(property_->identifier)){
            valuations_.find(property_->identifier)->second.value = fundamental_value_.value;
        }else{
            valuations_.emplace(property_->identifier, fundamental_value_);
        }
    }


    auto m = this->template create_message<kelly_bettor_ddsf>(
        message->sender, interval.lower, (*this), message->sender,
        interval.lower, interval.lower,nav_,demand);

    output_signal->put(interval.lower, 1.0);


    for(auto [p,q]: owner<stock>::properties.items){
        if(0 == q.amount){
            continue;
        }
        m->supply.emplace(p->identifier,std::make_tuple(q, quantity(0)));
    }

    for(auto [p,q]: owner<securities_lending_contract>::properties.items){
        if(0 == q.amount){
            continue;
        }
        m->supply.emplace(p->security,std::make_tuple(quantity(0), q));
    }

    return interval.lower;
}

std::string kelly_bettor::describe() const
{
    return "constant demand";
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
                result_.emplace(k,  scale_ * 0.5 * value_);
            }else{
                auto supply_long_  = double(std::get<0>(j->second));
                auto supply_short_ = double(std::get<1>(j->second));

                auto signal_ = 1.0;

                result_.emplace(k, scale_ * value_ - (supply_long_ - supply_short_) * (quoted_price_ * variable_)
                );
            }
        }
    }
    return result_;
}
