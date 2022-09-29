#pragma once
#include "marketData.hpp"
#include <ql/models/equity/hestonmodel.hpp>
#include <ql/models/equity/hestonmodelhelper.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/processes/hestonprocess.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/volatility/equityfx/hestonblackvolsurface.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/period.hpp>
#include <cmath>


#include <iostream>
#include <fstream>


using namespace std;
using namespace QuantLib;

namespace data {
    CalibrationSetup getCalibrationSetup(vector<Real> maturities, vector<Rate> iRates, vector<Real> divs, vector<vector<Real>> strikes, vector<vector<Volatility>> vols, Real spot) {

        Date settlementDate(Settings::instance().evaluationDate());
        DayCounter dayCounter = Actual365Fixed();
        Calendar calendar = TARGET();
        
        Handle<Quote> s0(ext::make_shared<SimpleQuote>(spot));



        auto const kNExpirations = maturities.size();

        vector<Date> dates;
        dates.push_back(settlementDate);
        for (Real maturity : maturities)
            dates.push_back(settlementDate + (int)maturity);


        iRates.insert(iRates.begin(), 0.0);
        divs.insert(divs.begin(), 0.0);

        Handle<YieldTermStructure> riskFreeTS(
            ext::make_shared<ZeroCurve>(dates, iRates, dayCounter));
        riskFreeTS->enableExtrapolation();
       
        Handle<YieldTermStructure> dividendYield(
            ext::make_shared<ZeroCurve>(dates, divs, dayCounter));
        dividendYield->enableExtrapolation();
       
        vector<ext::shared_ptr<CalibrationHelper> > options;
    

        cout << "Data to helpers" << endl;
        for (Size m = 0; m < kNExpirations; ++m) {
            for (Size s = 0; s < strikes[m].size(); ++s) {
                Handle<Quote> vol(ext::make_shared<SimpleQuote>(vols[m][s]));
                cout << strikes[m][s] << " ";
                Period maturity(maturities[m], Days);
                maturity.normalize();
                cout << maturity << " ";
                cout << vols[m][s] << endl;
                options.push_back(ext::make_shared<HestonModelHelper>(maturity, calendar,
                    spot, strikes[m][s], vol,
                    riskFreeTS, dividendYield,
                    BlackCalibrationHelper::RelativePriceError));
            }
        }

        CalibrationSetup marketData = { s0, riskFreeTS, dividendYield, options };

        return marketData;

    }

}