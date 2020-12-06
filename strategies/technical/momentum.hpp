/// \file   momentum_investor.hpp
///
/// \brief
///
/// \authors    maarten
/// \date       2020-02-13
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
#ifndef ME_momentum_HPP
#define ME_momentum_HPP


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

using namespace esl;

#include "../../fund.hpp"


class momentum
: public fund
{
public:
    explicit momentum( const esl::identity<fund> &i
                     , const jurisdiction &j = esl::law::jurisdictions::US
                     , size_t window = 1); 
    ///
    /// theta - 1
    ///
    size_t window;

    ///
    /// 
    ///
    esl::law::property_map<std::map<time_point, price>> historic_prices;

    /// 
    /// \param message 
    /// \param interval 
    /// \param seed 
    /// \return 
    time_point invest(std::shared_ptr<quote_message> message,
                      time_interval interval, std::seed_seq &seed) override;

    /// 
    ///
    ///
    /// \return 
    [[nodiscard]] std::string describe() const override
    {
        return "momentum trader";
    }
};



#include <esl/economics/markets/walras/differentiable_order_message.hpp>
#include <esl/economics/price.hpp>
using namespace esl;
using economics::markets::walras::differentiable_order_message;
///
///
///
struct momentum_ddsf
    : public differentiable_order_message
{
public:

    economics::price net_asset_value;

    ///
    /// \brief
    /// \param sender
    /// \param recipient
    /// \param sent
    /// \param received
    std::map<identity<law::property>, double> valuations;

    double agression;
    double leverage;

    momentum_ddsf
        ( const identity<agent> &sender
            , const identity<agent> &recipient
            , simulation::time_point sent     = simulation::time_point()
            , simulation::time_point received = simulation::time_point()
            , economics::price nav = economics::price(0, economics::currencies::USD)
            , std::map<identity<law::property>, double> valuations = {}
        );

    virtual std::map<identity<law::property>, variable>
    excess_demand(
        const std::map<identity<law::property>,
            std::tuple<economics::markets::quote, variable>>
        &quotes) const override;

    template<class archive_t>
    void serialize(archive_t &archive, const unsigned int version)
    {
        (void)version;
        archive &BOOST_SERIALIZATION_BASE_OBJECT_NVP(differentiable_order_message);
    }
};


#endif  // ME_momentum_HPP
