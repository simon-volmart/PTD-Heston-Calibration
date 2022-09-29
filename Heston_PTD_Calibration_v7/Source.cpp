#include <ql/qldefines.hpp>
#if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
#  include <ql/auto_link.hpp>
#endif

#include <ql/instruments/floatfloatswap.hpp>
#include <ql/instruments/floatfloatswaption.hpp>
#include <ql/instruments/nonstandardswaption.hpp>
#include <ql/instruments/payoffs.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/pricingengines/swaption/gaussian1dswaptionengine.hpp>
#include <ql/pricingengines/swaption/gaussian1dnonstandardswaptionengine.hpp>
#include <ql/pricingengines/swaption/gaussian1dfloatfloatswaptionengine.hpp>
#include <ql/models/shortrate/onefactormodels/gsr.hpp>
#include <ql/models/shortrate/onefactormodels/markovfunctional.hpp>
#include <ql/models/shortrate/calibrationhelpers/swaptionhelper.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/cashflows/lineartsrpricer.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/indexes/swap/euriborswap.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/swaption/swaptionconstantvol.hpp>
#include <ql/rebatedexercise.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/handle.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/time/calendars/unitedstates.hpp>
#include <ql/termstructures/volatility/equityfx/blackvariancesurface.hpp>
#include <ql/models/equity/hestonmodelhelper.hpp>
#include <ql/models/equity/hestonmodel.hpp>
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/math/optimization/differentialevolution.hpp>
#include <ql/math/interpolations/bicubicsplineinterpolation.hpp>
#include <ql/math/optimization/leastsquare.hpp>
#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/projectedcostfunction.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/experimental/math/fireflyalgorithm.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/models/equity/piecewisetimedependenthestonmodel.hpp>
#include <ql/pricingengines/vanilla/analyticptdhestonengine.hpp>
#include <ql/pricingengine.hpp>
#include <ql/models/equity/hestonmodel.hpp>
#include "marketData.hpp"

#include <iostream>
#include <memory>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <string>
#include<cmath>

using namespace std;
using namespace QuantLib;
using namespace data;

struct Engine_model {
    ext::shared_ptr<PricingEngine> engine;
    ext::shared_ptr<CalibratedModel> model;
};
struct Market_data {
    std::vector<Rate> rates;
    std::vector<Rate> divs;
    std::vector<vector<Real>> strikes;
    std::vector<vector<Real>> weights;
    std::vector<vector<Real>> vols;
    std::vector<Real> maturities;
};
struct HestonParameters {
    Real kappa;
    Real theta;
    Real sigma;
    Real v0;
    Real rho;
};

vector<Real>  calibrateWith(vector<ext::shared_ptr<CalibrationHelper> > options, Engine_model& engine_model, int scope, vector<bool> fix_params, vector<Real> weights);


Market_data fetchData(string path) {
    fstream newfile;
    Real maturity, strike, vol, weight, rate, div;

    std::vector<Rate> rates;
    std::vector<Rate> divs;
    std::vector<vector<Real>> strikes;
    std::vector<vector<Real>> weights;
    std::vector<vector<Real>> vols;
    std::vector<Real> maturities;

    newfile.open(path, ios::in);
    if (newfile.is_open()) {
        while (newfile >> maturity >> strike >> vol >> weight >> rate >> div) {
            if (maturities.empty() || maturities.back() != maturity) {
                maturities.push_back(maturity);
                rates.push_back(rate);
                divs.push_back(div);
                strikes.push_back(vector<Real>());
                weights.push_back(vector<Real>());
                vols.push_back(vector<Real>());
            }
            strikes.back().push_back(strike);
            weights.back().push_back(weight);
            vols.back().push_back(vol);
        }
        newfile.close();
    }
    Market_data result = { rates, divs, strikes, weights, vols, maturities };
    return result;
}

Engine_model get_PTD_engine(vector<HestonParameters>& params,
    Handle<YieldTermStructure> riskFreeTS,
    Handle<YieldTermStructure> dividendYield,
    Handle<Quote> s0,
    vector<Real> maturities);

Engine_model get_Analytic_engine(HestonParameters params,
    Handle<YieldTermStructure> riskFreeTS,
    Handle<YieldTermStructure> dividendYield,
    Handle<Quote> s0);

vector<Real> slice(vector<Real>& arr, int mat_pos);
vector<vector<Real>> slice(vector<vector<Real>>& arr, int mat_pos);




int main() {

    Real spot = 23208.25;
    Date settlementDate(01, August, 2022);

    Handle<Quote> s0(ext::make_shared<SimpleQuote>(spot));
    Settings::instance().evaluationDate() = settlementDate;

    string data_path = "C:\\Users\\Syimyk\\Desktop\\data_for_model\\2022-08-01_all_market.txt";

    std::vector<Real> weights_;

    Market_data data = fetchData(data_path);
    for (auto v : data.weights) {
        for (auto w : v) {
            weights_.push_back(w);
        }
    }
    HestonParameters t = { 
        2.98319,
        0.87133,
        4.84217,
        0.616182,
        - 0.260703 };
    vector<HestonParameters> temp_input_params = { t };
    vector<Real> temp_output_params;
    CalibrationSetup globalData = getCalibrationSetup(data.maturities, data.rates, data.divs, data.strikes, data.vols, spot);

    auto m_size = data.maturities.size();


    Market_data slice_of_data;

    for (int i = 0; i < m_size; i++) { //calibrate for each expiry date 
        cout << "================================================>" << endl << endl;
        slice_of_data.divs = slice(data.divs, i);
        slice_of_data.maturities = slice(data.maturities, i);
        slice_of_data.rates = slice(data.rates, i);
        slice_of_data.strikes = slice(data.strikes, i);
        slice_of_data.vols = slice(data.vols, i);

        CalibrationSetup temp_helpers = getCalibrationSetup(slice_of_data.maturities, slice_of_data.rates, slice_of_data.divs, slice_of_data.strikes, slice_of_data.vols, spot);

        std::cout << "\nInput Params: ";
        cout << "V0: " << temp_input_params[i].v0 << endl;
        for (const auto& p : temp_input_params) {
            cout << "Kappa: " << p.kappa << endl;
            cout << "theta: " << p.theta << endl;
            cout << "sigma: " << p.sigma << endl;
            cout << "rho: " << p.rho << endl;
        }

        std::vector<Real> temp_weights;
        for (const auto& v : slice(data.weights, i)) {
            for (auto w : v)
                temp_weights.push_back(w);
        }

        vector<bool> fix(temp_input_params.size() * 4 + 1, true);  //fix parameters for the previous expiries 
        for (int j = 1; j <= 4; j++) {
            fix[i * j + j - 1] = false;
        }
        if (i == 0)
            fix.back() = false; //fix v0 after the first expiry

        vector<Real> mat_;
        for (int j = 0; j <= i; j++) {
            mat_.push_back(data.maturities[j]);
        }
        Engine_model temp_PTD_engine;
        if (true)
            temp_PTD_engine = get_PTD_engine(temp_input_params, globalData.riskFreeTS, globalData.dividendYield, globalData.s0, mat_);
        else
            temp_PTD_engine = get_Analytic_engine(t, globalData.riskFreeTS, globalData.dividendYield, globalData.s0);


        temp_output_params = calibrateWith(temp_helpers.options, temp_PTD_engine, 1, fix, temp_weights);

        temp_input_params[i].kappa = temp_output_params[i * 2 + 1];
        temp_input_params[i].theta = temp_output_params[i];
        temp_input_params[i].sigma = temp_output_params[i * 3 + 2];
        temp_input_params[i].rho = temp_output_params[i * 4 + 3];
        temp_input_params[i].v0 = temp_output_params.back();

        temp_input_params.push_back(temp_input_params[i]); //last parameters as initial parameters for the next expiry
    }
    //global search with obtained parameters. (It may prioritize some expiries since weights had been assigned globally)  
    //requires lots of time --> minimal impact
    temp_input_params.pop_back();
    auto temp_PTD_engine = get_PTD_engine(temp_input_params, globalData.riskFreeTS, globalData.dividendYield, globalData.s0, data.maturities);
    temp_output_params = calibrateWith(globalData.options, temp_PTD_engine, 1, vector<bool>(), weights_); 
    for (auto p : temp_output_params) {
        cout << p << " ";
    }
    return 0;
}





vector<Real> calibrateWith(vector<ext::shared_ptr<CalibrationHelper>> options, Engine_model& engine_model,
    int scope, vector<bool> fix_params, vector<Real> weights) {

    for (const auto& option : options)
        ext::dynamic_pointer_cast<BlackCalibrationHelper>(option)->setPricingEngine(engine_model.engine);
    std::cout << "Calibration started..." << endl;


    switch (scope) {

    case 0: 
        {
            std::cout << "LevenbergMarquardt algorithm." << endl;
            EndCriteria end_criteria{ 500, 200, 1e-8, 1e-8, 1e-8 };
            LevenbergMarquardt lm(1e-12, 1e-12, 1e-12);
            try {
                engine_model.model->calibrate(options, lm, end_criteria, Constraint(), weights, fix_params);
            }
            catch (Error e) {
                std::cout << e.what() << endl;
            }
            break;
        }
    case 1:
        {
        std::cout << "Differential Evolution." << endl;
        EndCriteria end_criteria{ 300, 100, 1e-8, 1e-8, 1e-8 };
        DifferentialEvolution::Configuration de_conf;
        de_conf.withBounds(true)
            .withCrossoverProbability(0.8)
            .withPopulationMembers(200)
            .withStepsizeWeight(0.3)
            .withStrategy(DifferentialEvolution::BestMemberWithJitter)
            .withSeed(100);
        DifferentialEvolution de_optimizer(de_conf);

        HestonModelHelper w = (HestonModelHelper&)*options[0];
          
        vector<Real> arr1(weights.size(), -1);
        vector<Real> arr2(weights.size(), 1000);
        arr2.back() = 1;
        Array a1(arr1.begin(), arr1.end());
        Array a2(arr2.begin(), arr2.end());

        NonhomogeneousBoundaryConstraint c(a1, a2);  //can be deleted 
        try {
            engine_model.model->calibrate(options, de_optimizer, end_criteria, Constraint(), weights, fix_params);
        }
        catch (Error e) {
            std::cout << e.what() << endl;
        }
        break;
        }
    case 2: {
        Size agents = 150;
        Real vola = 1.5;
        Real intense = 1.0;
        ext::shared_ptr<FireflyAlgorithm::Intensity> intensity =
            ext::make_shared<ExponentialIntensity>(10.0, 1e-8, intense);
        ext::shared_ptr<FireflyAlgorithm::RandomWalk> randomWalk =
            ext::make_shared<LevyFlightWalk>(vola, 0.5, 1.0);

        std::cout << "Function eggholder, Agents: " << agents
            << ", Vola: " << vola << ", Intensity: " << intense << std::endl;
        FireflyAlgorithm fa(agents, intensity, randomWalk, 40);
        EndCriteria ec(5000, 1000, 1.0e-8, 1.0e-8, 1.0e-8);
        try {
            engine_model.model->calibrate(options, fa, ec, Constraint(), weights, fix_params);
        }
        catch (Error e) {
            std::cout << e.what() << endl;
        }

        break;
    }
    }






    std::cout << "Helpers' report: " << endl;
    Real avr_err = 0;
    Real max = 0;
    std::cout << "marketIV | modelIV | marketPrice | modelPrice | calibrationError | priceDifference" << endl;
    for (auto option : options) {
        Real difference = 0;
        HestonModelHelper w = (HestonModelHelper&)*option;
        try {
            std::cout << w.impliedVolatility(w.marketValue(), 1e-8, 500, 1e-8, 2.5) << " ";
            std::cout << w.impliedVolatility(w.modelValue(), 1e-8, 500, 1e-8, 2.5) << " ";
        }
        catch (Error e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        std::cout << w.marketValue() << " ";
        std::cout << w.modelValue() << " ";
        std::cout << w.calibrationError() << " ";
        difference = w.marketValue() - w.modelValue();
        std::cout << difference << endl;
        avr_err += w.calibrationError();
        if (max < abs(difference))
            max = abs(difference);
    }
    std::cout << "Max error: " << max << endl;
    std::cout << "Average rel error: " << avr_err/options.size() << endl;
    Array result = engine_model.model->params();
    vector<Real> output;
    for (auto t : result) {
        output.push_back(t);
    }

    return output;
}

Engine_model get_PTD_engine(
    vector<HestonParameters>& params,
    Handle<YieldTermStructure> riskFreeTS,
    Handle<YieldTermStructure> dividendYield,
    Handle<Quote> s0,
    vector<Real> maturities) {

    auto dc = riskFreeTS->dayCounter();
    Date settlementDate = Settings::instance().evaluationDate();

    vector<Time> times;
    for (auto m : maturities)
        times.push_back((double)dc.yearFraction(settlementDate, settlementDate + m));
    times.back() = times.back() + 10.0;

    const TimeGrid modelGrid(times.begin(), times.end());

    const std::vector<Time> pTimes(times.begin(), times.end() - 1);

    Real v0 = params[0].v0;
    PiecewiseConstantParameter kappa(pTimes, BoundaryConstraint(0, 50));
    PiecewiseConstantParameter sigma(pTimes, BoundaryConstraint(0, 10));
    PiecewiseConstantParameter theta(pTimes, BoundaryConstraint(0, 10));
    PiecewiseConstantParameter rho(pTimes, BoundaryConstraint(-1, 1));

    for (Size i = 0; i < kappa.size(); ++i) {
        kappa.setParam(i, params[i].kappa);
        theta.setParam(i, params[i].theta);
        sigma.setParam(i, params[i].sigma);
        rho.setParam(i, params[i].rho);
    }

    const ext::shared_ptr<PiecewiseTimeDependentHestonModel> model(
        ext::make_shared<PiecewiseTimeDependentHestonModel>(
            riskFreeTS, dividendYield, s0, v0, theta, kappa, sigma, rho, modelGrid));

    const ext::shared_ptr<AnalyticPTDHestonEngine> engine_(
        ext::make_shared<AnalyticPTDHestonEngine>(
            model));

    Engine_model result = { engine_, model };
    return result;
}
Engine_model get_Analytic_engine(HestonParameters params,
    Handle<YieldTermStructure> riskFreeTS,
    Handle<YieldTermStructure> dividendYield,
    Handle<Quote> s0) {

    auto heston_process = ext::make_shared<HestonProcess>(
        riskFreeTS,
        dividendYield,
        s0,
        params.v0,
        params.kappa,
        params.theta,
        params.sigma,
        params.rho
        );
    auto heston_model = ext::make_shared<HestonModel>(heston_process);
    auto heston_prcing_engine = ext::make_shared<AnalyticHestonEngine>(heston_model
        );
    Engine_model result = { heston_prcing_engine ,heston_model };
    return result;
}



vector<Real> slice(vector<Real>& arr,
    int mat_pos)
{
    vector<Real> result;
    result.push_back(arr[mat_pos]);
    return result;
}
vector<vector<Real>> slice(vector<vector<Real>>& arr, int mat_pos)
{
    vector<vector<Real>> result;
    result.push_back(arr[mat_pos]);
    return result;
}
