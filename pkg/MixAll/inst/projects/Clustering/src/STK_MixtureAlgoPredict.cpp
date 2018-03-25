/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_MixtureAlgo.cpp
 *  @brief In this file we implement the run method of the mixture algorithms.
 **/

#include <STKernel/include/STK_Exceptions.h>
#include "../include/MixtureAlgo/STK_MixtureAlgoPredict.h"
#include "../include/STK_IMixtureComposer.h"

namespace STK
{

bool EMPredict::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("-------------------------------\n");
  stk_cout << _T("Entering EMPredict::run() with:\n")
           << _T("nbIterMax_ = ") << nbIterMax_ << _T("\n")
           << _T("epsilon_ = ") << epsilon_ << _T("\n");
#endif
  if (p_model_->computeNbMissingValues() == 0) return predictWithNoMissingValues();
  bool result = true;
  try
  {
    if (!burnStep())
    {
      return false;
    }
    Real currentLnLikelihood = p_model_->lnLikelihood();
    int iter;
    for (iter = 0; iter < nbIterLong_; iter++)
    {
      p_model_->imputationStep();
      p_model_->eStep();
      Real lnLikelihood = p_model_->lnLikelihood();
      if ( (lnLikelihood - currentLnLikelihood) < epsilon_) // no abs as the likelihood should increase
      {
#ifdef STK_MIXTURE_VERY_VERBOSE
        stk_cout << _T("Terminating EMPredict::run() with:\n")
                 << _T("iter = ") << iter << _T("\n")
                 << _T("delta = ") << lnLikelihood - currentLnLikelihood << _T("\n");
#endif
        break;
      }
      currentLnLikelihood = lnLikelihood;
    }
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In EMPredict::run() iteration ") << iter << _T("terminated.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("An error occur in EMPredict::run():\n") << msg_error_ << _T("\n");
#endif
    result = false;
  }
  catch (Exception const& error)
  {
    msg_error_ = error.error();
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("An error occur in EMPredict::run():\n") << msg_error_ << _T("\n");
#endif
    result = false;
  }
  p_model_->finalizeStep();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Terminating EMPredict::run()\n");
  stk_cout << _T("----------------------------\n");
#endif
  return result;
}

bool SemiSEMPredict::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("------------------------------------\n");
  stk_cout << _T("Entering SemiSEMPredict::run() with:\n")
           << _T("nbIterMax_ = ") << nbIterMax_ << _T("\n")
           << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  bool result = true;
  try
  {
    if (!burnStep()) { return false;}
    int iter;
    for (iter = 0; iter < nbIterLong_; ++iter)
    {
      p_model_->samplingStep();
      p_model_->eStep();
      p_model_->storeIntermediateResults(iter+1); // store current parameters
    }
  }
  catch (Clust::exceptions const& error)
  {
    msg_error_ = Clust::exceptionToString(error);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("An exception occur in SemiSEMPredict::run(): ") << msg_error_ << _T("\n");
#endif
    result = false;
  }
  catch (Exception const& error)
  {
    msg_error_ = error.error();
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("An error occur in SemiSEMPredict::run():\n") << msg_error_ << _T("\n");
#endif
    result = false;
  }
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In SemiSEMPredict::run() iterations terminated.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  // set averaged parameters
  if (result)
  {
    p_model_->setParametersStep();
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("\nIn SemiSEMPredict::run(), setParameters done.\n")
             << _T("p_model_->lnLikelihood = ") << p_model_->lnLikelihood() << _T("\n");
#endif
  }
  else
  { p_model_->releaseIntermediateResults();}
  p_model_->finalizeStep();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Terminating SemiSEMPredict::run()\n");
  stk_cout << _T("---------------------------------\n");
#endif
  return result;
}


} // namespace STK
