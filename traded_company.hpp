/// \file   traded_company.hpp
///
/// \brief  In the market ecology model traded companies generate fictitious revenues.
///
/// \authors    Maarten
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
#ifndef ME_TRADED_COMPANY_HPP
#define ME_TRADED_COMPANY_HPP

#include <esl/economics/company.hpp>
using esl::economics::company;
using esl::law::jurisdiction;
using esl::simulation::time_point;
using esl::simulation::time_interval;
using esl::economics::finance::dividend_policy;
using esl::economics::price;


///
/// \brief  A traded company issues shares and pays dividends
///
class traded_company
: public company
, public esl::identifiable_as<traded_company>
{
public:
    ///
    /// \param i    company identity
    /// \param j    jurisdiction (for accounting)
    /// \param sample   seed for random dividend process
    traded_company( const esl::identity<traded_company> &i
                  , const jurisdiction &j
                  , uint64_t sample = 0
                  , time_point end_point = 201 * 252
                  );

    ///
    /// \brief  Generates new profits and decides on next dividend
    ///
    /// \param interval
    /// \param s
    /// \return
    time_point act(time_interval interval, std::seed_seq &s) override;

    ///
    /// \brief  Estimates the first upcoming dividend payment
    ///
    /// \param interval
    /// \param seed
    /// \return
    std::optional<dividend_policy> upcoming_dividend( time_interval interval
                                                    , std::seed_seq &seed
                                                    ) override;

    ///
    /// \brief  Number of days in a year. Used to annualize economic quantities
    ///
    constexpr static size_t day_count  = 252;

    ///
    /// \brief  2% growth rate annually, converted to daily
    ///
    double mu = pow(1.02, 1./day_count) - 1;

    ///
    /// \brief  10% volatility annually, converted to daily
    ///
    double sigma = 0.10 / sqrt(double(day_count));

    ///
    /// \brief  autocorrelation lag of 1 Month = 252 / 12
    ///
    constexpr static size_t tau        = day_count / 12;

    ///
    /// \brief  autocorrelation factor
    ///
    constexpr static double omega      = 0.1;

    ///
    /// \brief  Wiener process realization
    ///
    std::vector<double> wiener;

    ///
    /// \brief  Raw dividends per share
    /// TODO: convert to price instead of double
    ///
    std::vector<double> dividends_per_share;

    ///
    /// \brief  current dividend per share
    /// TODO: convert to price instead of double
    ///
    double dividend_per_share;

    ///
    /// \brief  company income
    ///
    [[nodiscard]] esl::economics::price income() const
    {
        return esl::economics::price(0, primary_jurisdiction.tender);
    }

    ///
    /// \brief  company book value
    ///
    [[nodiscard]] esl::economics::price book_value() const
    {
        return esl::economics::price(0, primary_jurisdiction.tender);
    }

    ///
    /// \brief  debt
    ///
    [[nodiscard]] esl::economics::price debt() const
    {
        return esl::economics::price(0, primary_jurisdiction.tender);
    }

    ///
    /// \brief sales for current time
    /// \return
    [[nodiscard]] esl::economics::price sales() const
    {
        return esl::economics::price(0, primary_jurisdiction.tender);
    }

    ///
    /// \brief costs for current time
    /// \return
    [[nodiscard]] esl::economics::price costs() const
    {
        return esl::economics::price(0, primary_jurisdiction.tender);
    }

    ///
    /// \brief earnings at current time
    /// \return
    [[nodiscard]] esl::economics::price historic_earnings(/*time span*/) const
    {
        return esl::economics::price(0, primary_jurisdiction.tender);
    }
};


#endif //ME_TRADED_COMPANY_HPP
