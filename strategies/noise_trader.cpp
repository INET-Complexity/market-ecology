/// \file   brownian_motion.cpp
///
/// \brief
///
/// \authors    maarten
/// \date       2020-01-21
/// \copyright  Copyright 2017-2020 The Institute for New Economic Thinking,
/// Oxford Martin School, University of Oxford
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
#include "noise_trader.hpp"

#include <utility>

#include <esl/data/log.hpp>
using economics::finance::securities_lending_contract;

using namespace esl;
using namespace esl::economics;
using namespace esl::economics::markets;

#include "allocation.hpp"


noise_trader::noise_trader(const identity<fund> &i, const jurisdiction &j)
: agent(i)
, owner<cash>(i)
, owner<stock>(i)
, fund(i, j)
{

}

time_point noise_trader::invest(
    std::shared_ptr<walras::quote_message> message,
    simulation::time_interval interval, std::seed_seq &seed)
{
    if(message->received < interval.lower){
        return interval.lower;
    }

    auto nav_ = net_asset_value(interval);
    if(nav_.value <= 0){
        //return interval.upper;
    }

    //auto st = accounting::standard(currencies::USD);
    std::default_random_engine generator(seed);

    if(variates.empty()) {
        std::normal_distribution<double> distribution_(1.0, 0.0000001);
        std::map<identity<law::property>, double> variates_;
        for(auto [k, v]: message->proposed) {
            auto noise_ = distribution_(generator);
            variates_.insert({k->identifier, noise_});
        }
        variates.insert(std::make_pair(interval.lower, variates_));
    } else {
        auto iterator_ = variates.rbegin();
        std::normal_distribution<double> distribution_(0.0, 1.0);

        double sigma_ = 0.03;
        double rho_ = 0.0003;
        double mu_ = 0.01;

        double intercept_ = std::pow(1.0002, interval.lower) - 1;


        if(iterator_->first != interval.lower) {
            std::map< identity< law::property>, double> variates_;
            for(auto [k, v] : message->proposed) {
                auto j  = iterator_->second.find(k->identifier);
                auto Yt = j->second;

                auto zi = distribution_(generator);

                auto dt = std::sqrt(interval.lower - iterator_->first);

                auto E1 = rho_ * ( std::log2(mu_ + intercept_) - std::log2(Yt)) * Yt * dt;
                auto E2 = sigma_ * Yt * zi;

                output_signal->put(interval.lower,  Yt + E1 + E2);

                if(iterator_->second.end() == j) {
                    throw std::logic_error("stochastic process was started with no history");
                }
                variates_.insert({k->identifier, Yt + E1 + E2 });
            }
            variates.insert(std::make_pair(interval.lower, variates_));
        }
    }
    //LOG(trace) << describe() << " " << identifier << " inventory " <<  inventory << " @" << nav_<< std::endl;
    auto m = this->template create_message<brownian_motion_ddsf>(
        message->sender, interval.lower, (*this),
        message->sender,
        interval.lower,
        interval.lower,
        nav_,
        variates.rbegin()->second);

    for(const auto &[p,q]: owner<stock>::properties.items){
        if(0 == q.amount){
            continue;
        }
        //LOG(warning) << describe() << this->identifier << " is reporting long equal to " << q << " for " << p->identifier <<  std::endl;
        m->supply.insert(std::make_pair(identity<property>(p->identifier),std::make_tuple(q, 0)));
    }

    for(const auto &[p,q]: owner<securities_lending_contract>::properties.items){
        if(0 == q.amount){
            continue;
        }
        //LOG(warning) << describe() << this->identifier << " is reporting short equal to " << q << " for " << p->identifier <<  std::endl;
        m->supply.emplace(p->security,std::make_tuple(0, q));
    }

    return interval.lower;
}

std::string noise_trader::describe() const
{
    return "noise trader";
}

///
///
///
brownian_motion_ddsf::brownian_motion_ddsf(
    const identity<agent> &sender,
    const identity<agent> &recipient,
    simulation::time_point sent, simulation::time_point received,
    price net_asset_value,
    std::map<identity<law::property>, double> variates
    )
: differentiable_order_message(sender, recipient, sent, received)
, variates(std::move(variates))
, net_asset_value(net_asset_value)
{ }

std::map<identity<law::property>, esl::variable>
brownian_motion_ddsf::excess_demand(
    const std::map<identity<law::property>,
                   std::tuple<quote, variable>> &quotes)
    const {
    std::map<identity<property>, esl::variable> result_;
    auto scale_ = double(net_asset_value) / quotes.size();

    for (auto &[k, v] : quotes) {
        const auto &[quote_, variable_] = v;
        const auto quoted_price_ = double(std::get<price>(quote_.type));
        auto i = variates.find(k);
        if (variates.end() == i) {
            continue;
        }


        if (variates.end() != i) {
            auto value_ = i->second;
            auto j = this->supply.find(k);
            if (supply.end() == j) {
                result_.emplace(k, scale_ / (quoted_price_ * variable_) * (value_ - (quoted_price_ * variable_)) - 0);
            } else {
                auto supply_long_ = double(std::get<0>(j->second));
                auto supply_short_ = double(std::get<1>(j->second));

                result_.emplace(k, (scale_ * 0.1 * value_)   -
                                   (supply_long_ + supply_short_) * quoted_price_ * variable_
                );
            }
        }
    }
    return result_;
}