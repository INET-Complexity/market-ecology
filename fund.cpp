/// \file   fund.cpp
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
#include "fund.hpp"

#include <esl/economics/markets/walras/price_setter.hpp>

#include <esl/economics/markets/walras/quote_message.hpp>
using esl::economics::markets::walras::quote_message;

using esl::law::owner;
using namespace esl::economics;
using namespace esl::economics::accounting;
using namespace esl::economics::finance;



fund::fund(const identity<fund> &i, const jurisdiction &j)
    : agent(i)
    , law::owner<cash>(i)
    , owner<stock>(i)
    , identifiable_as<fund>()
    , company(i,j)
    , lookup_(primary_jurisdiction.tender)
    , reset_amount(1.00, currencies::USD)
{
    output_net_asset_value  = create_output<price>("net_asset_value");
    output_signal           = create_output<double>("signal");
    output_cash             = create_output<price>("cash");
    output_stocks           = create_output<price>("stocks");
    output_loans            = create_output<price>("loans");
    output_lending          = create_output<price>("lending");
    output_pnl              = create_output<price>("pnl");

    auto invest_ = [this](auto msg, simulation::time_interval ti, std::seed_seq &seed) {
        return invest(msg, ti, seed);
    };

    ESL_REGISTER_CALLBACK(quote_message, -100, invest_, "make investment decisions");

    auto process_dividends_ = [this](std::shared_ptr<dividend_announcement_message> m,
                                     simulation::time_interval step,
                                     std::seed_seq &) {

        for(const auto &[share_, sharedetails_] : m->policy.dividend_per_share){
            std::tuple<identity<company>, share_class> key_ =
                {esl::reinterpret_identity_cast<company>(m->sender), share_};

            auto iterator_ = stocks.find(key_);
            if (stocks.end() == iterator_){
                LOG(warning) << "fund.cpp(" << __LINE__ << ") fund receives info about stock it is not tracking" << std::endl;
            }else{
                identity<property> stock_ = iterator_->second;
                market_data.shares_outstanding[stock_] = std::get<0>(sharedetails_);

                uint64_t stocks_ = 0;
                for(const auto& s: owner<stock>::inventory){
                    if(s.first->identifier == stock_){
                        stocks_ += s.second.amount;
                    }
                }

                auto dps_ = std::get<1>(sharedetails_);
                auto usd_ = std::make_shared<cash>(currencies::USD);
                auto [i, b] = owner<cash>::inventory.emplace(usd_, quantity(0));

                /// TODO: this needs to be moved someplace else
                auto dividend_received_ = uint64_t((stocks_ * dps_.value * (1.+risk_free_rate)) /  std::get<0>(sharedetails_));
                i->second.amount += dividend_received_;

                for(auto [p, q]: owner<securities_lending_contract>::properties.items){
                    if(p->security == stock_){
                        auto pay_on_short_ = uint64_t((q.amount * dps_.value * (1.+risk_free_rate)) / std::get<0>(sharedetails_));

                        if(i->second.amount < pay_on_short_){
                            LOG(trace) << "pay " << pay_on_short_ << " to fund short position " << std::endl;
                            throw std::invalid_argument("out of cash " + identifier.representation());
                        }
                        i->second.amount -= pay_on_short_;
                    }
                }

                if(market_data.dividend_per_share.end() == market_data.dividend_per_share.find(stock_)){
                    market_data.dividend_per_share.insert(std::make_pair(stock_,std::get<1>(sharedetails_)));
                }else{
                    market_data.dividend_per_share[stock_].value = std::get<1>(sharedetails_).value;
                }
            }
        }
        return step.upper;
    };

    ESL_REGISTER_CALLBACK(quote_message, 16, [&](std::shared_ptr<quote_message> m, simulation::time_interval ti, std::seed_seq &s){

        for(auto [property_, quote_]: m->proposed) {
            auto price_ = std::get<price>(quote_.type);
            lookup_.mark_to_market.erase(property_->identifier);
            auto i = lookup_.mark_to_market.emplace(property_->identifier, price_);
        }

            return ti.upper;
        }, "update prices");

    ESL_REGISTER_CALLBACK(dividend_announcement_message, 0, process_dividends_, "process dividend announcement and store dividends")
}

template<typename T>
price valuate(const esl::law::property_map<quantity> &inventory, const standard &s)
{
    price value_ = price::approximate(0.00, s.reporting_currency);
    for(const auto &[p,q]: inventory){
        auto cast_ = std::dynamic_pointer_cast<T>(p);
        if(cast_){
            value_ += s.value(*cast_, q);
        }
    }
    return value_;
}


///
/// If we are re-setting wealth in an experiment then
/// `target_net_asset_value` is set to the value that we reset wealth to.
/// `target_date` determines the last date on which we reset, for example
/// in experiments where we always reset it is set to 200 years
///
void fund::reset_wealth(price &net_asset_value_, simulation::time_interval ti)
{
    // if the target NAV rule does not apply, exit
    // if we are no longer resetting, exit
    if(!target_net_asset_value.has_value() || double(target_net_asset_value.value()) <= 0. || ti.lower > target_date){
        if(! output_net_asset_value->values.empty() ){
            auto pnl = net_asset_value_ - std::get<1>( output_net_asset_value->values.back() );
            output_pnl->put(ti.lower, pnl);
        }
        return;
    }

    // bailout as a fraction of assets and liabilities
    double bailout_ratio = double(target_net_asset_value.value()) / double(net_asset_value_);

    //std::cout << this->describe() << " bailout_ratio " << std::setprecision(5) << bailout_ratio << " (" << target_net_asset_value.value() << " / " << net_asset_value_ << ")" <<std::endl;
    // bailout monetary amount
    auto bailout_ = price::approximate((double(target_net_asset_value.value()) - double(net_asset_value_)), currencies::USD);
    auto remainder_ = bailout_.value;

    output_pnl->put(ti.lower, -bailout_);



    for(auto &[p, q]: inventory) {
        auto cast_ = std::dynamic_pointer_cast<loan>(p);
        if(cast_) {
            if(remainder_ > 0) {
                auto dec = std::min<std::uint64_t>(inventory[p].amount, remainder_ / 100);
                inventory[p].amount -= dec;  // std::uint64_t(inventory[p].amount * bailout_ratio);
                remainder_ -= dec;
            } else if(remainder_ < 0) {
                inventory[p].amount += (-remainder_/ 100);  // std::uint64_t(inventory[p].amount * bailout_ratio);
                remainder_ = 0;
            }
        }
    }

    if(remainder_){

    for(auto &[p, q]: inventory) {
        auto cast_ = std::dynamic_pointer_cast<cash>(p);
        if(cast_){
            if(remainder_ < 0){
                auto dec = std::min<std::uint64_t>(inventory[p].amount, -remainder_);
                inventory[p].amount -= dec;//std::uint64_t(inventory[p].amount * bailout_ratio);
                remainder_ += dec;
            }else if(remainder_ > 0){
                inventory[p].amount += remainder_;//std::uint64_t(inventory[p].amount * bailout_ratio);
                remainder_ = 0;
            }
        }
    }



    }



    if(remainder_){
        for(auto &[p, q]: inventory) {
            auto cast_ = std::dynamic_pointer_cast<stock>(p);
            if(cast_){
                auto i = lookup_.mark_to_market.find(p->identifier);
                if(lookup_.mark_to_market.end() == i){
                    continue;
                }
                auto change_ = remainder_ / std::get<price>(i->second.type).value;
                if(change_ < 0){
                    auto dec = std::min<std::uint64_t>(inventory[p].amount, -change_);
                    inventory[p].amount -= dec;//std::uint64_t(inventory[p].amount * bailout_ratio);
                    remainder_ += change_ * std::get<price>(i->second.type).value;
                }else if(change_ > 0){
                    inventory[p].amount += change_;//std::uint64_t(inventory[p].amount * bailout_ratio);
                    remainder_ = 0;
                }
            }


        }
    }
}


void fund::apply_reinvestment(price &net_asset_value_, simulation::time_interval ti)
{
    ///
    /// \brief  This section relates to the reinvestment rate experiments,
    ///         where we introduce (unseen) retail investors that
    ///         invest based on the performance of the fund
    ///
    if(previous_net_asset_value.has_value() && reinvestment_rate != 1.){
        double returns_ =  double(net_asset_value_) / previous_net_asset_value.value() - 1;
        double compound_returns_f_ = std::pow(1 + returns_, reinvestment_rate - 1) - 1;

        //std::cout <<"returns: " << returns_ << " compounded " << compound_returns_f_ << std::endl;

        // OLD rule matching papers
        //int64_t change_ = static_cast<int64_t>(round(  (reinvestment_rate -1.) * (double(net_asset_value_) - previous_net_asset_value.value()) ));
        // new rule

        double change_ = compound_returns_f_ * double(net_asset_value_);

        for(auto &[p, q]: inventory){
            auto cast_ = std::dynamic_pointer_cast<cash>(p);
            if(cast_){
                net_asset_value_ += price::approximate(change_, net_asset_value_.valuation);
                auto qchange_ = std::max<int64_t>(-((int64_t)inventory[p].amount), int64_t(change_ * 100) );
                inventory[p].amount += qchange_;
            }
        }
    }

}


price fund::net_asset_value(esl::simulation::time_interval ti)
{
    for(auto p: shareholder::prices){
        lookup_.mark_to_market.emplace(p.first->identifier, p.second);
    }

    for(auto p: bond_prices){
        lookup_.mark_to_market.emplace(p.first->identifier, p.second);
    }

    for(auto p: owner<stock>::properties.items){
        if(lookup_.mark_to_market.end() == lookup_.mark_to_market.find(p.first->identifier)){
            lookup_.mark_to_market.emplace(p.first->identifier, price::approximate(100.00, currencies::USD));
        }
    }

    auto cash1_ = valuate<cash>(inventory, lookup_);
    auto stocks_ = valuate<stock>(inventory, lookup_);
    output_stocks->put(ti.lower, stocks_);

    auto loans_value_ = valuate<loan>(inventory, lookup_);
    auto loans_ = price::approximate(double(loans_value_) * (1.+risk_free_rate), loans_value_.valuation);
    output_loans->put(ti.lower, loans_);

    auto lending_ = valuate<securities_lending_contract>(inventory, lookup_);
    output_lending->put(ti.lower, lending_);
    auto cash_ = price::approximate(double(cash1_) * (1.+risk_free_rate), cash1_.valuation);
    auto net_asset_value_ = cash_ + stocks_ + lending_ + loans_;

    output_cash->put(ti.lower, cash_);
    net_asset_value_ = cash_ + stocks_ + lending_ + loans_;

    if(ti.lower > target_date) {
        apply_reinvestment(net_asset_value_, ti);
    }

    reset_wealth(net_asset_value_, ti);
    previous_net_asset_value = double(net_asset_value_);

    output_net_asset_value->put(ti.lower, net_asset_value_);
    return net_asset_value_;
}

time_point fund::act(time_interval interval, std::seed_seq &s)
{
    (void)s;
    return interval.upper;
}

#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT(esl::data::output<esl::economics::price>);