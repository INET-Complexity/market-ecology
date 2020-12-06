/// \file   dividend_discount.hpp
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
#ifndef ME_DIVIDEND_DISCOUNT_HPP
#define ME_DIVIDEND_DISCOUNT_HPP


#include <esl/economics/interest_rate.hpp>

namespace esl::economics::markets::walras {
    struct quote_message;
}

using esl::economics::nominal_interest_rate;
using esl::economics::price;
using esl::economics::markets::walras::quote_message;
using esl::simulation::time_point;
using esl::simulation::time_interval;
using esl::economics::nominal_interest_rate;

#include "../../fund.hpp"

///
/// \brief  A value investor that uses discounting to value a stream of future
///         dividends.
///
class dividend_discount
: public fund
{
public:
    explicit dividend_discount(const esl::identity<fund> &i, const jurisdiction &j = esl::law::jurisdictions::US);

    law::property_map<std::vector<nominal_interest_rate>> dividend_yields;

    time_point invest(std::shared_ptr<quote_message> message,
                      time_interval interval, std::seed_seq &seed) override;

    [[nodiscard]] std::string describe() const override
    {
        return "value investor (dividend)";
    }
};


#include <esl/economics/markets/walras/differentiable_order_message.hpp>
#include <esl/economics/price.hpp>
using namespace esl;
using esl::economics::markets::walras::differentiable_order_message;
///
///
///
struct dividend_discount_ddsf
: public differentiable_order_message
{
public:

    esl::economics::price net_asset_value;

    ///
    /// \brief
    /// \param sender
    /// \param recipient
    /// \param sent
    /// \param received
    std::map<identity<law::property>, economics::price> valuations;

    double agression;
    double leverage;

    dividend_discount_ddsf
            ( const identity<esl::agent> &sender
            , const identity<esl::agent> &recipient
            , simulation::time_point sent     = simulation::time_point()
            , simulation::time_point received = simulation::time_point()
            , esl::economics::price nav = esl::economics::price(0.00, esl::economics::currencies::USD)
            , std::map<esl::identity<esl::law::property>, economics::price> valuations = {}
            );

    virtual std::map<esl::identity<esl::law::property>, esl::variable>
    excess_demand(
        const std::map<esl::identity<esl::law::property>,
                       std::tuple<esl::economics::markets::quote, esl::variable>>
            &quotes) const override;

    template<class archive_t>
    void serialize(archive_t &archive, const unsigned int version)
    {
        (void)version;
        archive &BOOST_SERIALIZATION_BASE_OBJECT_NVP(differentiable_order_message);
    }
};

#endif //ME_DIVIDEND_DISCOUNT_HPP
