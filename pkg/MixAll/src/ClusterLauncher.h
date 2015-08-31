/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2015  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 4 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file ClusterLauncher.h
 *  @brief In this file we define the ClusterLauncher class which
 *  construct properly a mixture model.
 **/


#ifndef CLUSTERLAUNCHER_H
#define CLUSTERLAUNCHER_H

#include "RDataHandler.h"

/**The ClusterLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
class ClusterLauncher : public STK::IRunnerBase
{
  public:
    /** constructor.
     * @param model a reference on the current model
     * @param strategy the strategy defined in R
     **/
    ClusterLauncher( SEXP model, SEXP nbCluster, SEXP modelNames, SEXP strategy, SEXP critName );
    /** constructor with a list of component.
     * @param model a reference on the current model
     * @param strategy the strategy defined in R
     **/
    ClusterLauncher( SEXP model, SEXP nbCluster, SEXP strategy, SEXP critName );
    /** destructor. */
    virtual ~ClusterLauncher();
    /** run the estimation */
    bool run();
    /** @return the model */
    inline Rcpp::S4 const& s4_model() const { return s4_model_;}

  protected:
    /** strategy from the R side */
    Rcpp::S4              s4_model_;
    /** strategy from the R side */
    Rcpp::S4              s4_strategy_;
    /** vector with the number of cluster to try */
    Rcpp::IntegerVector   v_nbCluster_;
    /** vector with the model names to try */
    Rcpp::CharacterVector v_modelNames_;
    /** character string with the model selection criterion name */
    std::string           criterion_;

  private:
    /** Select the best model among the modelNames and nbCluster given.
     *  @return the value of the best criteria.
     **/
    STK::Real selectSingleBestModel();
    /** Select the best model among the modelNames and nbCluster given.
     *  @return the value of the best criteria.
     **/
    STK::Real selectHeteroBestModel();
    /** create the mixtures in the given composer */
    void createMixtures(STK::MixtureComposer* p_composer);
    /** create the kernel mixtures in the given composer */
    void createKernelMixtures(STK::MixtureComposer* p_composer);
    /** get the parameters */
    void getParameters(Rcpp::S4& s4_component, std::string const& idData);
    /** get the diagonal Gaussian parameters */
    void getDiagGaussianParameters(Rcpp::S4& s4_component, std::string const& idData);
    /** get the Poisson parameters */
    void getPoissonParameters(Rcpp::S4& s4_component, std::string const& idData);
    /** get the gamma parameters */
    void getGammaParameters(Rcpp::S4& s4_component, std::string const& idData);
    /** get the gamma parameters */
    void getCategoricalParameters(Rcpp::S4& s4_component, std::string const& idData);
    /** get the kernel parameters */
    void getKernelParameters(Rcpp::S4& s4_component, std::string const& idData);

    /** data handler */
    RDataHandler handler_;

    /** diagonal Gaussian mixture models manager */
    STK::DiagGaussianMixtureManager<RDataHandler> diagGaussianManager_;
    /** Poisson mixture models manager */
    STK::PoissonMixtureManager<RDataHandler> poissonManager_;
    /** gamma mixture models manager */
    STK::GammaMixtureManager<RDataHandler> gammaManager_;
    /** categorical mixture models manager */
    STK::CategoricalMixtureManager<RDataHandler> categoricalManager_;
    /** kernel mixture models manager */
    STK::KernelMixtureManager<RDataHandler> kernelManager_;

    /** pointer on the main composer */
    STK::IMixtureComposer* p_composer_;

    /** Is the model with heterogeneous data ? */
    bool isHeterogeneous_;
};

#endif /* CLUSTERLAUNCHER_H */
