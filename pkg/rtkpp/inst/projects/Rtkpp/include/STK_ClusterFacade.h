/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 4 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_ClusterFacade.h
 *  @brief In this file we define the ClusterFacade which construct properly a mixture model.
 **/


#ifndef STK_CLUSTERFACADE_H
#define STK_CLUSTERFACADE_H

#include <Rcpp.h>
#include "STK_RcppTraits.h"

namespace STK
{

class IMixtureComposer;
class IMixtureStrategy;

/** facade design pattern.
 * The ClusterFacade allow to create the strategy for estimating a mixture model
 * with less effort
 **/
class ClusterFacade
{
  public:
    /** constructor.
     * @param p_model a reference on the current model
     * @param R_strategy the strategy defined in R
     **/
    inline ClusterFacade( IMixtureComposer*& p_model)
                        : p_model_(p_model), p_strategy_(0)
    {}
    /** destructor. */
    ~ClusterFacade();
    /** set model in case it is needed after construction */
    inline void setModel(IMixtureComposer*& p_model) {p_model_ = p_model;};
    /** create a FullStrategy
     *  @param R_strategy the strategy defined in the R environment */
    void createFullStrategy(Rcpp::S4 R_strategy);
    /** run the strategy */
    bool run();

  protected:
    /** the mixture model to estimate */
    IMixtureComposer*& p_model_;
    /** the strategy to use in order to estimate the mixture model */
    IMixtureStrategy* p_strategy_;
};


}  // namespace STK

#endif /* STK_CLUSTERFACADE_H */
