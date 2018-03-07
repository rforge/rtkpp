/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016 Serge Iovleff

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
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file  STK_CategoricalParameters.cpp
 *  @brief In this file we implement the Parameters classes for categorical
 *  mixture models
 **/

#include "../include/CategoricalMixtureModels/STK_CategoricalParameters.h"

namespace STK
{

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Categorical_pjk_>::ModelParameters(int nbCluster): proba_(nbCluster) {}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Categorical_pjk_>::ModelParameters( ModelParameters const& param)
               : proba_(param.proba_)
{}
/* copy operator.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Categorical_pjk_>& ModelParameters<Clust::Categorical_pjk_>::operator=( ModelParameters const& param)
{
  proba_ = param.proba_;
  return *this;
}
/* destructor */
ModelParameters<Clust::Categorical_pjk_>::~ModelParameters() {}
/* resize the set of parameter */
void ModelParameters<Clust::Categorical_pjk_>::resize(Range const& rangeModalities, Range const& rangeCols)
{
  for (int k = proba_.begin(); k< proba_.end(); ++k)
  { proba_[k].resize(rangeModalities, rangeCols) = 1./rangeModalities.size();}
}

/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Categorical_pk_>::ModelParameters(int nbCluster): proba_(nbCluster) {}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Categorical_pk_>::ModelParameters( ModelParameters const& param)
               : proba_(param.proba_)
{}
/* destructor */
ModelParameters<Clust::Categorical_pk_>::~ModelParameters() {}
/* copy operator.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Categorical_pk_>& ModelParameters<Clust::Categorical_pk_>::operator=( ModelParameters const& param)
{
  proba_ = param.proba_;
  return *this;
}

/* resize the set of parameter */
void ModelParameters<Clust::Categorical_pk_>::resize(Range const& rangeModalities, Range const& range)
{
  for (int k = proba_.begin(); k< proba_.end(); ++k)
  { proba_[k].resize(rangeModalities) = 1./rangeModalities.size();}
}

} // namespace STK

