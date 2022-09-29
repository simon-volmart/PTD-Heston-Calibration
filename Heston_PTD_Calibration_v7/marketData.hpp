#pragma once

#include <ql/handle.hpp>
#include <ql/quote.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/models/calibrationhelper.hpp>

#include <iostream>


using namespace std;
using namespace QuantLib;

namespace data {
    struct DataToCalibrate {
        Handle<Quote> s0;
        Handle<YieldTermStructure> riskFreeTS, dividendYield;
        std::vector<ext::shared_ptr<CalibrationHelper> > options;
    };
    DataToCalibrate getDataToCalibrate(vector<Real> maturities, vector<Rate> iRates, vector<Real> divs, vector<vector<Real>> strikes, vector<vector<Volatility>> vols, Real spot);
}



