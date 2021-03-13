/// \file   momentum_investor.cpp
///
/// \brief
///
/// \authors    Maarten P. Scholl
/// \date       2020-03-05
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
#ifndef ME_MEAN_REVERTING_NOISE_HPP
#define ME_MEAN_REVERTING_NOISE_HPP

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




class mean_reverting_noise_trader
: public fund
{
public:
    explicit mean_reverting_noise_trader( const identity<fund> &i
                                        , const jurisdiction &j = law::jurisdictions::US
                                        , size_t seed = 0
                                         );

    law::property_map<std::vector<nominal_interest_rate>> dividend_yields;

    std::map<time_point, std::map<identity<law::property>, double>> variates;

    size_t seed;

    double sigma_   = 0.2 * std::sqrt(1./252);

    time_point invest(std::shared_ptr<quote_message> message,
                      time_interval interval, std::seed_seq &seed) override;

    [[nodiscard]] std::string describe() const override
    {
        return "mean reverting noise trader";
    }
};


#include <esl/economics/markets/walras/differentiable_order_message.hpp>
#include <esl/economics/price.hpp>
using namespace esl;
using economics::markets::walras::differentiable_order_message;
///
///
///
struct mean_reverting_noise_ddsf
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
    std::map<identity<law::property>, economics::price> valuations;

    double agression;
    double leverage;

    mean_reverting_noise_ddsf
            ( const identity<agent> &sender
                    , const identity<agent> &recipient
                    , simulation::time_point sent     = simulation::time_point()
                    , simulation::time_point received = simulation::time_point()
                    , economics::price nav = economics::price(0, economics::currencies::USD)
                    , std::map<identity<law::property>, economics::price> valuations = {}
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



#endif //ME_MEAN_REVERTING_NOISE_HPP
