/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 4 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_ClusterLauncher.cpp
 *  @brief In this file we implement the ClusterLauncher which
 *  construct properly a mixture model.
 **/


#include "../inst/projects/MixAll/ClusterLauncher.h"
#include "../inst/projects/MixAll/ClusterFacade.h"

using namespace Rcpp;

namespace STK
{
/* facade design pattern.
 * The ClusterLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ClusterLauncher::ClusterLauncher( SEXP model, SEXP nbCluster, SEXP models, SEXP strategy, SEXP critName )
                                : IRunnerBase()
                                , s4_model_(model)
                                , s4_strategy_(strategy)
                                , v_nbCluster_(nbCluster)
                                , v_models_(models)
                                , criterion_(Rcpp::as<std::string>(critName))
                                , handler_()
                                , diagGaussianManager_(handler_)
                                , poissonManager_(handler_)
                                , gammaManager_(handler_)
                                , categoricalManager_(handler_)
                                , kernelManager_(handler_)
                                , p_composer_(0)
                                , isMixedData_(false)
{}
/* facade design pattern.
 * The ClusterLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ClusterLauncher::ClusterLauncher( SEXP model, SEXP nbCluster, SEXP strategy, SEXP critName )
                                : IRunnerBase()
                                , s4_model_(model)
                                , s4_strategy_(strategy)
                                , v_nbCluster_(nbCluster)
                                , v_models_()
                                , criterion_(Rcpp::as<std::string>(critName))
                                , handler_()
                                , diagGaussianManager_(handler_)
                                , poissonManager_(handler_)
                                , gammaManager_(handler_)
                                , categoricalManager_(handler_)
                                , kernelManager_(handler_)
                                , p_composer_(0)
                                , isMixedData_(true)
{}
/* destructor. */
ClusterLauncher::~ClusterLauncher() { if (p_composer_) delete p_composer_;}

/* run the estimation */
bool ClusterLauncher::run()
{
  // compute the best model
  STK::Real initCriter = s4_model_.slot("criterion");
  STK::Real criter = (isMixedData_) ? selectMixedBestModel()
                                        : selectSingleBestModel();

  // get result common part of the estimated model
  s4_model_.slot("criterion")      = criter;
  s4_model_.slot("nbCluster")      = p_composer_->nbCluster();
  s4_model_.slot("lnLikelihood")   = p_composer_->lnLikelihood();
  s4_model_.slot("nbFreeParameter")= p_composer_->nbFreeParameter();
  s4_model_.slot("pk")             = Rcpp::wrap(p_composer_->pk());
  s4_model_.slot("tik")            = Rcpp::wrap(p_composer_->tik());
  s4_model_.slot("zi")             = Rcpp::wrap(p_composer_->zi());
  NumericVector fi = s4_model_.slot("lnFi");
  NumericVector zi = s4_model_.slot("zi");
  for (int i=0; i< fi.length(); ++i)
  {
    fi[i] = p_composer_->computeLnLikelihood(i);
    zi[i] += (1 - STK::baseIdx);  // set base 1 for the class labels
  }
  if (criter == initCriter || !STK::Arithmetic<STK::Real>::isFinite(criter)) return false;
  return true;
}

/* get the parameters */
STK::Real ClusterLauncher::selectSingleBestModel()
{
  std::string idDataBestModel;
  // component
  Rcpp::S4 s4_component = s4_model_.slot("component");

  // wrap data matrix with Rcpp and wrap Rcpp matrix with STK++ matrix
  NumericMatrix m_data = s4_component.slot("data");
  STK::RMatrix<double> data(m_data);

  double critValue = s4_model_.slot("criterion");
  int nbSample     = s4_model_.slot("nbSample");
  STK::IMixtureComposer*  p_current =0;
  STK::IMixtureCriterion* p_criterion =0;

  try
  {
    // create criterion
    if (criterion_ == "BIC") { p_criterion = new STK::BICMixtureCriterion();}
    if (criterion_ == "AIC") { p_criterion = new STK::AICMixtureCriterion();}
    if (criterion_ == "ICL") { p_criterion = new STK::ICLMixtureCriterion();}

    // start the estimation process, should end with the best model according to
    // the criteria
    ClusterFacade facade(p_current);
    facade.createFullStrategy(s4_strategy_);
    // loop over all the models
    for (int l=0; l <v_models_.size(); ++l)
    {
      std::string idData = "model" + STK::typeToString<int>(l);
      std::string idModel(as<std::string>(v_models_[l]));
      // transform R model names to STK++ model names
      // check have been done on the R side so.... Let's go
      bool freeProp;
      STK::Clust::Mixture model = STK::Clust::stringToMixture(idModel, freeProp);
      // add Data set with the new model name, m_data is just a pointer on a SEXP
      // structure thus there is no difficulties in doing so
      handler_.addData(m_data, idData, STK::Clust::mixtureToString(model));
      // for the current model,
      for (int k=0; k <v_nbCluster_.length(); ++k)
      {
        int K = v_nbCluster_[k];
        // create composer
        if (freeProp) { p_current = new STK::MixtureComposer(nbSample, K);}
        else          { p_current = new STK::MixtureComposerFixedProp(nbSample, K);}

        // create current mixture and register it
        std::string idData = "model" + STK::typeToString<int>(l);
        createMixtures(static_cast<STK::MixtureComposer*>(p_current));

        // run estimation and get results if possible
        if (!facade.run()) { msg_error_ += facade.error();}
        // compute criterion and update model if necessary
        p_criterion->setModel(p_current);
        p_criterion->run();
        if (critValue > p_criterion->value())
        {
          if (p_composer_) { std::swap(p_current, p_composer_);}
          else             { p_composer_ = p_current; p_current = 0;}
          s4_component.slot("modelName") = idModel;
          idDataBestModel = idData;
          critValue = p_criterion->value();
        }
        // release current composer
        if (p_current) { delete p_current; p_current = 0;}
      }
    }
    // release
    delete p_criterion;
    // get specific parameters
    getParameters(s4_component, idDataBestModel);
    return critValue;
  }
  catch (STK::Exception const& e)
  {
    if (p_current) delete p_current;
    if (p_criterion) delete p_criterion;
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return STK::Arithmetic<STK::Real>::max();
}

/* select best mixed data model */
STK::Real ClusterLauncher::selectMixedBestModel()
{
  // list of the component
  Rcpp::List s4_list =s4_model_.slot("ldata");
  STK::Real criter =s4_model_.slot("criterion");
  int nbSample =s4_model_.slot("nbSample");
  // main pointer
  STK::IMixtureComposer* p_current =0;
  STK::IMixtureCriterion* p_criterion =0;
  try
  {
    // create facade and strategy
    ClusterFacade facade(p_current);
    facade.createFullStrategy(s4_strategy_);
    bool sameProp = true;
    // loop over the list of component and fil handler_
    for (int l=0; l <s4_list.size(); ++l)
    {
      // component
      Rcpp::S4 s4_component = s4_list[l];
      NumericMatrix m_data = s4_component.slot("data");
      // id of the data set and of the model
      std::string idData  = "model" + STK::typeToString<int>(l);
      std::string idModel = s4_component.slot("modelName");
       // register
      bool freeMixture;
      STK::Clust::Mixture model = STK::Clust::stringToMixture(idModel, freeMixture);
      handler_.addData(m_data, idData, STK::Clust::mixtureToString(model));
      // if one of the model is free proportion, then we use free proportion
      sameProp &= (!freeMixture);
    }
    // create criterion
    if (criterion_ == "BIC") { p_criterion = new STK::BICMixtureCriterion();}
    if (criterion_ == "AIC") { p_criterion = new STK::AICMixtureCriterion();}
    if (criterion_ == "ICL") { p_criterion = new STK::ICLMixtureCriterion();}

    for (int k=0; k <v_nbCluster_.length(); ++k)
    {
      int K = v_nbCluster_[k];
      // create composer
      if (!sameProp) { p_current = new STK::MixtureComposer(nbSample, K);}
      else           { p_current = new STK::MixtureComposerFixedProp(nbSample, K);}
      // create all mixtures
      createMixtures(static_cast<STK::MixtureComposer*>(p_current));
      // run estimation and get results if possible
      if (!facade.run()) { msg_error_ += facade.error();}
      // compute criterion and update model if necessary
      p_criterion->setModel(p_current);
      p_criterion->run();
      if (criter > p_criterion->value())
      {
        if (p_composer_) { std::swap(p_current, p_composer_);}
        else             { p_composer_ = p_current; p_current = 0;}
        criter = p_criterion->value();
      }
      // release current composer
      if (p_current) { delete p_current; p_current = 0;}
    }
    // release
    delete p_criterion;
    // get parameters
    for (int l=0; l <s4_list.size(); ++l)
    {
      // component
      Rcpp::S4 s4_component = s4_list[l];
      // id of the data set and of the model
      std::string idData  = "model" + STK::typeToString<int>(l);
      getParameters(s4_component, idData);
    }
    //
    return criter;
  }
  catch (STK::Exception const& e)
  {
    if (p_current) delete p_current;
    if (p_criterion) delete p_criterion;
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return STK::Arithmetic<STK::Real>::max();
}

/* create the mixtures in the given composer */
void ClusterLauncher::createMixtures(STK::MixtureComposer* p_composer)
{
  p_composer->createMixture(diagGaussianManager_);
  p_composer->createMixture(poissonManager_);
  p_composer->createMixture(gammaManager_);
  p_composer->createMixture(categoricalManager_);

  createKernelMixtures(p_composer);
}

/* create the kernel mixtures in the given composer */
void ClusterLauncher::createKernelMixtures(STK::MixtureComposer* p_composer)
{
  typedef STK::MixtureComposer::ConstMixtIterator ConstMixtIterator;
  p_composer->createMixture(kernelManager_);
  // loop over
  for (ConstMixtIterator it =  p_composer->v_mixtures().begin(); it != p_composer->v_mixtures().end(); it++)
  {
    std::string idData = (*it)->idData();
    std::string idModel;
    handler_.getIdModelName(idData, idModel);
    STK::Clust::Mixture typeModel = STK::Clust::stringToMixture(idModel);
    if (STK::Clust::mixtureToMixtureClass(typeModel) == STK::Clust::Kernel_)
    {
      Rcpp::S4 s4_component = s4_model_.slot("component");
      STK::RVector<double> dim((SEXP)s4_component.slot("dim"));
      double kdim = (dim.size()>0) ? dim[0] : 10;
      kernelManager_.setDim(p_composer->getMixture(idData), kdim);
    }
  }
}

/* fill the s4_component with the parameters */
void ClusterLauncher::getParameters(Rcpp::S4& s4_component, std::string const& idData)
{
  std::string rModelName = s4_component.slot("modelName");
  bool freeProp;
  switch (STK::Clust::mixtureToMixtureClass(STK::Clust::stringToMixture(rModelName, freeProp)))
  {
    case STK::Clust::Gaussian_:
      getDiagGaussianParameters(s4_component, idData);
      break;
    case STK::Clust::Poisson_:
      getPoissonParameters(s4_component, idData);
      break;
    case STK::Clust::Gamma_:
      getGammaParameters(s4_component, idData);
      break;
    case STK::Clust::Categorical_:
      getCategoricalParameters(s4_component, idData);
      break;
    case STK::Clust::Kernel_:
      getKernelParameters(s4_component, idData);
      break;
    case STK::Clust::unknown_mixture_class_:

    default:
      break;
  }
}

/* get the diagonal Gaussian parameters */
void ClusterLauncher::getDiagGaussianParameters(Rcpp::S4& s4_component, std::string const& idData)
{
  // get parameters
  STK::ArrayXX params;
  static_cast<STK::MixtureComposer*>(p_composer_)->getParameters(diagGaussianManager_,idData, params);
  // get dimensions
  int K = params.sizeRows()/2, nbVariable = params.sizeCols();
  // get results
  STK::ArrayXX mean(K, nbVariable), sigma(K, nbVariable);
  for (int k=0; k<K; ++k)
  {
    mean.row(k)  = params.row(2*k);
    sigma.row(k) = params.row(2*k+1);
  }
  // save results in s4_model
  s4_component.slot("mean")  = Rcpp::wrap(mean);
  s4_component.slot("sigma") = Rcpp::wrap(sigma);
  // get data
  s4_component.slot("data") = diagGaussianManager_.getData<double>(idData).matrix();
}

/* get the kernel parameters */
void ClusterLauncher::getKernelParameters(Rcpp::S4& s4_component, std::string const& idData)
{
  // get parameters
  STK::ArrayXX param;
  static_cast<STK::MixtureComposer*>(p_composer_)->getParameters(kernelManager_,idData, param);
  // save results in s4_model
  s4_component.slot("sigma") = Rcpp::wrap(param.col(0));
  s4_component.slot("dim")   = Rcpp::wrap(param.col(1));
  // get data -- not necessary for kernels--
  //s4_component.slot("data") = kernelManager_.getData<double>(idData).matrix();
}

/* get the diagonal Gaussian parameters */
void ClusterLauncher::getPoissonParameters(Rcpp::S4& s4_component, std::string const& idData)
{
  // get parameters
  STK::ArrayXX params;
  static_cast<STK::MixtureComposer*>(p_composer_)->getParameters(poissonManager_,idData, params);
  // save results in s4_model
  s4_component.slot("lambda")  = Rcpp::wrap(params);
  // get data
  s4_component.slot("data") = poissonManager_.getData<double>(idData).matrix();
}

/* get the gamma parameters */
void ClusterLauncher::getGammaParameters(Rcpp::S4& s4_component, std::string const& idData)
{
  // get parameters
  STK::ArrayXX params;
  static_cast<STK::MixtureComposer*>(p_composer_)->getParameters(gammaManager_,idData, params);
  // get dimensions
  int K = params.sizeRows()/2, nbVariable = params.sizeCols();
  // get results
  STK::ArrayXX shape(K, nbVariable), scale(K, nbVariable);
  for (int k=0; k<K; ++k)
  {
    shape.row(k) = params.row(2*k);
    scale.row(k) = params.row(2*k+1);
  }
  // save results in s4_model
  s4_component.slot("shape") = Rcpp::wrap(shape);
  s4_component.slot("scale") = Rcpp::wrap(scale);
  // get data
  s4_component.slot("data") = gammaManager_.getData<double>(idData).matrix();
}

/* get the diagonal Categorical parameters */
void ClusterLauncher::getCategoricalParameters(Rcpp::S4& s4_component, std::string const& idData)
{
  // get parameters
  STK::ArrayXX params;
  static_cast<STK::MixtureComposer*>(p_composer_)->getParameters(categoricalManager_,idData, params);
  params.shift(0,0);
  // save results in s4_model
  s4_component.slot("plkj") = Rcpp::wrap(params);
  // get data
  s4_component.slot("data") = categoricalManager_.getData<int>(idData).matrix();
}

} // namespace STK

