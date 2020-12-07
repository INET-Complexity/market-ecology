/// \file   fund.hpp
///
/// \brief  
///
/// \authors    maarten
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
#ifndef ME_FUND_HPP
#define ME_FUND_HPP

#include <esl/economics/finance/company.hpp>
using namespace esl;
using namespace esl::economics;
using esl::economics::company;
using esl::economics::cash;
using esl::economics::price;
using esl::economics::finance::bond;
using esl::economics::finance::stock;
using esl::identity;
using esl::law::property;
using esl::law::jurisdiction;
using esl::simulation::time_point;
using esl::simulation::time_interval;
#include <esl/law/jurisdictions.hpp>
using esl::law::jurisdictions::US;
#include <esl/economics/markets/walras/quote_message.hpp>
using esl::economics::markets::walras::quote_message;
#include <esl/economics/owner.hpp>
#include <esl/economics/finance/securities_lending_contract.hpp>
#include <esl/economics/price.hpp>

#include "financial_ratios.hpp"


///
/// \brief  Defines an investment fund specific to the market ecology model
///
class fund
: public company
, public identifiable_as<fund>
, public law::owner<finance::loan>  // for leverage
, public law::owner<finance::securities_lending_contract>   // for short selling
{
public:
    ///
    /// \brief
    ///
    /// \param i    identity
    /// \param j    jurisdiction
    explicit fund(const identity<fund> &i, const jurisdiction &j = US);

    ///
    /// \brief  Computes the net-asset value for the fund
    ///
    /// \param ti
    /// \return
    price net_asset_value(simulation::time_interval ti);

    ///
    /// \brief  Until which date wealth is reset. Used during the transient
    ///
    simulation::time_point target_date = std::numeric_limits<simulation::time_point >::max();

    ///
    /// \param interval
    /// \param s
    /// \return
    time_point act(time_interval interval, std::seed_seq &s) override;

    ///
    /// \param message
    /// \param seed
    /// \return
    virtual time_point invest( std::shared_ptr<quote_message> message
        , time_interval, std::seed_seq &seed) = 0;

    ///
    /// \brief  During transient, the value to reset the NAV to
    ///
    std::optional<price> target_net_asset_value;

    ///
    /// \brief When the fund goes bankrupt, it is reset to this amount
    ///
    price reset_amount;

    ///
    /// \brief  For feedback on paper:
    ///         funds get reinvestment_rate * PNL for the next time step
    ///
    double reinvestment_rate = 1.0;


    std::optional<price> previous_net_asset_value;


    ///
    /// \brief  Used to arbitrarily scale up/down trade signal
    ///
    double aggression = 1.;

    ///
    /// \brief  Maximum leverage allowed (long or short)
    ///
    double maximum_leverage = 1.;

    ///
    /// \brief  Data structure to store information about the investments
    ///
    financial_ratios market_data;


    // outputs:

    std::shared_ptr<data::output<price>>    output_net_asset_value;
    std::shared_ptr<data::output<double>>   output_signal;
    std::shared_ptr<data::output<price>>    output_cash;
    std::shared_ptr<data::output<price>>    output_stocks;
    std::shared_ptr<data::output<price>>    output_loans;
    std::shared_ptr<data::output<price>>    output_lending;
    std::shared_ptr<data::output<price>>    output_pnl;

    ///
    /// \brief  Used to describe the agent when debugging
    /// \return
    [[nodiscard]] std::string describe() const override
    {
        std::stringstream stream_;
        stream_ << "fund " << identifier;
        return stream_.str();
    }
};


#endif //ME_FUND_HPP
