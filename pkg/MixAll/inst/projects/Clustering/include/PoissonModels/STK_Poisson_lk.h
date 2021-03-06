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

/** @file STK_Poisson_lk.h
 *  @brief In this file we implement the Poisson_lk class
 **/

#ifndef STK_POISSON_LK_H
#define STK_POISSON_LK_H

#include "../PoissonModels/STK_PoissonBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Poisson_lk;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the Poisson_lk traits policy. */
template<class Array_>
struct MixtureTraits< Poisson_lk<Array_> >
{
  typedef Array_                  Array;
  typedef typename Array::Type    Type;
  /** Type of the structure storing the parameters of a MixturePoisson_lk model*/
  typedef ModelParameters<Clust::Poisson_lk_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  The Poisson mixture model @c Poisson_lk has a probability function of the form
 * \f[
 *    P(\mathbf{x}=(n_1,\ldots,n_d)|\theta)
 *     = \sum_{k=1}^K p_k \prod_{j=1}^d e^{-\lambda_{k}} \frac{\lambda_{k}^{n_j}}{n_j!}.
 * \f]
 **/
template<class Array>
class Poisson_lk: public PoissonBase<Poisson_lk<Array> >
{
  public:
    typedef PoissonBase<Poisson_lk<Array> > Base;
    using Base::param_;
    using Base::p_data;
    

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Poisson_lk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Poisson_lk( Poisson_lk const& model): Base(model) {}
    /** destructor */
    ~Poisson_lk() {}
    /** @return the value of lambda of the kth cluster and jth variable */
    inline Real lambda(int k, int j) const { return param_.lambda_[k];}
    /** Initialize randomly the parameters of the Poisson mixture. */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** Compute the weighted probabilities. */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const { return this->nbCluster();}
};

/* Initialize randomly the parameters of the Poisson mixture. */
template<class Array>
void Poisson_lk<Array>::randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  Real m = p_data()->template cast<Real>().mean();
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  { param_.lambda_[k] = Law::Exponential::rand(m);}
}


/* Compute the modalities probabilities */
template<class Array>
bool Poisson_lk<Array>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.lambda_[k]= 0.;
    for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
    {
      for (int i=p_tik->beginRows(); i < p_tik->endRows(); ++i)
      { param_.lambda_[k] += p_data()->elt(i,j) * p_tik->elt(i,k);}
    }
    param_.lambda_[k] /= p_data()->sizeCols()*p_tk->elt(k);
  }
  return true;
}

} // namespace STK

#endif /* STK_Poisson_LK_H */
