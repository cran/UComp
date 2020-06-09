// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// UCompC
SEXP UCompC(SEXP commands, SEXP ys, SEXP us, SEXP models, SEXP periodss, SEXP rhoss, SEXP hs, SEXP tTests, SEXP criterions, SEXP ps, SEXP rubbish2s, SEXP rubbishs, SEXP verboses, SEXP stepwises, SEXP estimOks, SEXP p0s, SEXP vs, SEXP yFitVs, SEXP nonStationaryTermss, SEXP rubbish3s, SEXP harmonicss, SEXP criterias, SEXP cycleLimitss, SEXP betass, SEXP typeOutlierss);
RcppExport SEXP _UComp_UCompC(SEXP commandsSEXP, SEXP ysSEXP, SEXP usSEXP, SEXP modelsSEXP, SEXP periodssSEXP, SEXP rhossSEXP, SEXP hsSEXP, SEXP tTestsSEXP, SEXP criterionsSEXP, SEXP psSEXP, SEXP rubbish2sSEXP, SEXP rubbishsSEXP, SEXP verbosesSEXP, SEXP stepwisesSEXP, SEXP estimOksSEXP, SEXP p0sSEXP, SEXP vsSEXP, SEXP yFitVsSEXP, SEXP nonStationaryTermssSEXP, SEXP rubbish3sSEXP, SEXP harmonicssSEXP, SEXP criteriasSEXP, SEXP cycleLimitssSEXP, SEXP betassSEXP, SEXP typeOutlierssSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type commands(commandsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< SEXP >::type us(usSEXP);
    Rcpp::traits::input_parameter< SEXP >::type models(modelsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type periodss(periodssSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rhoss(rhossSEXP);
    Rcpp::traits::input_parameter< SEXP >::type hs(hsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type tTests(tTestsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type criterions(criterionsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ps(psSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rubbish2s(rubbish2sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rubbishs(rubbishsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type verboses(verbosesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type stepwises(stepwisesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type estimOks(estimOksSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p0s(p0sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type vs(vsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yFitVs(yFitVsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nonStationaryTermss(nonStationaryTermssSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rubbish3s(rubbish3sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type harmonicss(harmonicssSEXP);
    Rcpp::traits::input_parameter< SEXP >::type criterias(criteriasSEXP);
    Rcpp::traits::input_parameter< SEXP >::type cycleLimitss(cycleLimitssSEXP);
    Rcpp::traits::input_parameter< SEXP >::type betass(betassSEXP);
    Rcpp::traits::input_parameter< SEXP >::type typeOutlierss(typeOutlierssSEXP);
    rcpp_result_gen = Rcpp::wrap(UCompC(commands, ys, us, models, periodss, rhoss, hs, tTests, criterions, ps, rubbish2s, rubbishs, verboses, stepwises, estimOks, p0s, vs, yFitVs, nonStationaryTermss, rubbish3s, harmonicss, criterias, cycleLimitss, betass, typeOutlierss));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_UComp_UCompC", (DL_FUNC) &_UComp_UCompC, 25},
    {NULL, NULL, 0}
};

RcppExport void R_init_UComp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
