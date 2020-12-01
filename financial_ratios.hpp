/// \file   market_data.hpp
///
/// \brief
///
/// \authors    Maarten P. Scholl
/// \date       2020-01-31
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
#ifndef ME_FINANCIAL_RATIOS_HPP
#define ME_FINANCIAL_RATIOS_HPP

#include <unordered_map>

#include <esl/economics/finance/company.hpp>
using namespace esl;

///
/// \brief  A structure to hold information about stocks
///
struct financial_ratios
{
    ///
    /// \brief  For each stock, how many shares are outstanding in the market
    ///
    std::unordered_map<identity<law::property>, std::uint64_t> shares_outstanding;

    ///
    /// \brief  Latest dividend per share certificate underlying stock
    ///
    std::unordered_map<identity<law::property>, economics::price> dividend_per_share;

    ///
    /// \brief  Book value of the company underlying the stock
    ///
    std::unordered_map<identity<law::property>, economics::price> book_value;
};


#endif  // ME_FINANCIAL_RATIOS_HPP
