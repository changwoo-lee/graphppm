// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// uppertri_pcm
IntegerMatrix uppertri_pcm(IntegerMatrix zsamples, int n, int nsamples);
RcppExport SEXP _graphppm_uppertri_pcm(SEXP zsamplesSEXP, SEXP nSEXP, SEXP nsamplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type zsamples(zsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    rcpp_result_gen = Rcpp::wrap(uppertri_pcm(zsamples, n, nsamples));
    return rcpp_result_gen;
END_RCPP
}
// update_uppertri_pcm
void update_uppertri_pcm(IntegerMatrix uppertri_pcm, IntegerVector z, int n);
RcppExport SEXP _graphppm_update_uppertri_pcm(SEXP uppertri_pcmSEXP, SEXP zSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type uppertri_pcm(uppertri_pcmSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    update_uppertri_pcm(uppertri_pcm, z, n);
    return R_NilValue;
END_RCPP
}
// update_uppertri_logpsm
void update_uppertri_logpsm(NumericMatrix uppertri_logpsm, IntegerVector z, int n, double logweight);
RcppExport SEXP _graphppm_update_uppertri_logpsm(SEXP uppertri_logpsmSEXP, SEXP zSEXP, SEXP nSEXP, SEXP logweightSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type uppertri_logpsm(uppertri_logpsmSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type logweight(logweightSEXP);
    update_uppertri_logpsm(uppertri_logpsm, z, n, logweight);
    return R_NilValue;
END_RCPP
}
// runcppWilson
Rcpp::IntegerMatrix runcppWilson(int n, Rcpp::IntegerMatrix graph_edge_list, int rootvertex, Rcpp::IntegerVector ordering, int maxwalk);
RcppExport SEXP _graphppm_runcppWilson(SEXP nSEXP, SEXP graph_edge_listSEXP, SEXP rootvertexSEXP, SEXP orderingSEXP, SEXP maxwalkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type graph_edge_list(graph_edge_listSEXP);
    Rcpp::traits::input_parameter< int >::type rootvertex(rootvertexSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ordering(orderingSEXP);
    Rcpp::traits::input_parameter< int >::type maxwalk(maxwalkSEXP);
    rcpp_result_gen = Rcpp::wrap(runcppWilson(n, graph_edge_list, rootvertex, ordering, maxwalk));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_graphppm_uppertri_pcm", (DL_FUNC) &_graphppm_uppertri_pcm, 3},
    {"_graphppm_update_uppertri_pcm", (DL_FUNC) &_graphppm_update_uppertri_pcm, 3},
    {"_graphppm_update_uppertri_logpsm", (DL_FUNC) &_graphppm_update_uppertri_logpsm, 4},
    {"_graphppm_runcppWilson", (DL_FUNC) &_graphppm_runcppWilson, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_graphppm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
