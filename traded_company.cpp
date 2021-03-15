/// \file   traded_company.cpp
///
/// \brief  
///
/// \authors    Maarten P. Scholl
/// \date       2020-01-20
/// \copyright  Copyright 2017-2020 The Institute for New Economic Thinking, Oxford Martin School, University of Oxford
///
///             Licensed under the Apache License, Version 2.0 (the "License");
///             you may not use this file except in compliance with the License.
///             You may obtain a copy of the License at
///
///                 http://www.apache.org/licenses/LICENSE-2.0
///
///             Unless required by applicable law or agreed to in writing, software
///             distributed under the License is distributed on an "AS IS" BASIS,
///             WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
///             See the License for the specific language governing permissions and
///             limitations under the License.
///
///             You may obtain instructions to fulfill the attribution requirements in CITATION.cff
///
#include "traded_company.hpp"
using esl::economics::finance::dividend_announcement_message;

#include <random>

using esl::identity;
using esl::economics::cash;
using esl::economics::price;
using esl::economics::accounting::inventory_filter;
using esl::economics::finance::stock;
using esl::economics::finance::share_class;
using esl::law::owner;
using esl::economics::currencies::USD;
using esl::law::property;
using esl::interaction::transfer;
using esl::economics::accounting::inventory_filter;
using esl::economics::finance::dividend_policy;

///
/// \param i
/// \param j
traded_company::traded_company( const identity<traded_company> &i
                              , const jurisdiction &j
                              , uint64_t sample
                              , time_point end_point
                              )
: agent(i)
, owner<cash>(i)
, owner<stock>(i)
, identifiable_as<traded_company>()
, company(i,j)
, dividend_per_share(0.01)
{
    create_output<double>("dividend_per_share");

    std::seed_seq seed_ {
            std::uint64_t(std::hash<identity<agent>>()(i))
            , sample + 2
    };

    std::minstd_rand0 generator_(seed_);

    unsigned int dividend_receiving_shares_ = 0;
    for(const auto &[s, q] : shares_outstanding) {
        if(s.dividend) {
            dividend_receiving_shares_ += q;
        }
    }

    constexpr double dt = 1.; // discretization step

    for(size_t t = 0; t < end_point; ++t){
        std::normal_distribution<double> standard_normal_(0,1);
        auto z = standard_normal_(generator_);

        if(wiener.size() > tau){
            // zero based index, so subtract another 1 with tau
            z = (1.0 - omega) * z + omega * wiener[wiener.size() - 1 - tau];
        }
        wiener.push_back(z);

        // reflect so that dividend is always positive
        dividend_per_share = abs(dividend_per_share
                                 + mu * dividend_per_share * dt
                                 + sigma * dividend_per_share * wiener.back());

        dividends_per_share.push_back(dividend_per_share);
    }
}

///
/// \param interval
/// \param s
/// \return
time_point traded_company::act(time_interval interval, std::seed_seq &seed)
{
    auto possible_policy_ = upcoming_dividend(interval, seed);
    auto next_event_      = interval.lower + 1;

    if(!possible_policy_.has_value()) {
        return next_event_;
    }

    auto policy_ = possible_policy_.value();

    if(interval.lower >= policy_.announcement_date) {
        if(last_announced_ < policy_.announcement_date) {
            last_announced_ = policy_.announcement_date;

            for(const auto &s : unique_shareholders()) {
                this->template create_message<dividend_announcement_message>(
                        s, interval.lower, this->identifier, s, policy_);
            }
        }
    } else {
        next_event_ = std::min<time_point>(
                policy_.announcement_date, next_event_);
    }

    if(interval.lower >= policy_.payable_date) {
        if(last_payment_ < policy_.announcement_date) {
            last_payment_ = policy_.announcement_date;
        }
    } else {
        next_event_ = std::min<time_point>(
                policy_.announcement_date, next_event_);
    }
    return next_event_;
}

///
/// \param interval
/// \return
std::optional<dividend_policy> traded_company::upcoming_dividend(time_interval interval, std::seed_seq &seed)
{
    std::minstd_rand0 generator_(seed);

    unsigned int dividend_receiving_shares_ = 0;
    for(const auto &[s, q] : shares_outstanding) {
        if(s.dividend) {
            dividend_receiving_shares_ += q;
        }
    }


    outputs["dividend_per_share"]->write(std::make_tuple(interval.lower,
                                                         dividends_per_share[interval.lower]));
    if(interval.lower > 1 && dividends_per_share[interval.lower] != dividends_per_share[interval.lower-1]){
        //std::cout << interval << " DIV!" << std::endl;
    }
    auto unappropriated_profit_ =
        esl::economics::cash(USD).price( dividends_per_share[interval.lower] *
                                        dividend_receiving_shares_);
    //LOG(trace) << "company " << identifier << " pays out "
    //           << unappropriated_profit_
    //           << " in profits, distributed over "
    //           << dividend_receiving_shares_ << " shares, "
    //           << dividend_per_share
    //           << " (" << (dividend_per_share / (0.02 - mu)) << ")"
    //           << std::endl;

    auto dividends_ =
        compute_dividend_per_share(unappropriated_profit_);

    dividend_policy policy_(
            interval.lower
            , interval.lower
            , interval
            , interval.lower
            , primary_jurisdiction.tender
            , dividends_
            );

    return policy_;
}
