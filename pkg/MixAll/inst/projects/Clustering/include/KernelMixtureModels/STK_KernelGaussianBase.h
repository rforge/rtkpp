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
 * created on: Oct 24, 2014
 * Author:   Serge Iovleff
 **/

/** @file STK_KernelGaussianBase.h
 *  @brief In this file we define the KernelGaussianBase base class
 **/

#ifndef STK_KERNELGAUSSIANBASE_H
#define STK_KERNELGAUSSIANBASE_H

#include "STK_KernelParameters.h"
#include "../STK_IMixtureDensity.h"

#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>
#include <STatistiK/include/STK_Kernel_IKernel.h>

namespace STK
{

/** @ingroup Clustering
 *  Base class for the Gaussian Kernel Mixture Models.
 **/
template<class Derived>
class KernelGaussianBase: public IMixtureDensity<Derived>
{
  public:
    typedef IMixtureDensity<Derived > Base;
    typedef typename hidden::MixtureTraits< Derived >::Array Array;

    using Base::param_;
    using Base::nbSample;
    using Base::nbCluster;

  protected:
    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    KernelGaussianBase( int nbCluster): Base(nbCluster), p_kernel_(0) {}
    /** copy constructor
     *  @param model Model to copy
     **/
    KernelGaussianBase( KernelGaussianBase const& model)
                            : Base(model)
                             , p_kernel_(model.p_kernel_)
                             , dik_(model.dik_)
    {}
    /** destructor */
    ~KernelGaussianBase() {}

  public:
    // access
    /** @return the value of sigma2 in the kth cluster */
    inline Real sigma2(int k) const { return param_.sigma2(k);}
    /** @return the value of dim in the kth */
    inline Real dim(int k) const { return param_.dim(k);}
    /** @return the pointer on the kernel */
    inline Kernel::IKernel const* p_kernel() const { return p_kernel_;}
    /** @return distance of the individual to the kth centroid */
    inline CArrayXX dik() const { return dik_;}

    // setter
    /** set the dimensions of the kernel mixture model using an unique value */
    inline void setDim(Real const& dim)  { param_.dim_ = dim;}
    /** set the dimension of the kernel mixture model */
    inline void setDim(CPointX const& dim)  { param_.dim_ = dim;}

    /** set the dimensions of the kernel mixture model using an unique value.
     *  call to this method will triger the call to @c initializeModelImpl
     *  (as setData will not be called)
     **/
    inline void setKernel(Kernel::IKernel const* p_kernel)
    {
      p_kernel_ = p_kernel;
      initializeModel();
    }
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param pk the probabilities of each class for the ith individual
     *  @note this method is not used as there is no missing values to impute
     **/
    template<class Weights>
    inline Real impute(int i, int j, Weights const& pk) const { return 0.;}
    /** @return a simulated value for the jth variable of the ith sample
     *  @note this method is not used as there is no missing values to simulate
     **/
    inline Real rand(int i, int j, int k) const { return 0.;}

    /** This function is used in order to get the current values of the
     *  parameters in an array.
     *  @param[out] params the array with the parameters of the mixture.
     */
    template<class ArrayParam>
    void getParameters(ArrayParam& params) const;
    /** This function can be used to write summary of parameters to the output stream.
     *  @param p_tik a constant pointer on the posterior probabilities
     *  @param os Stream where you want to write the summary of parameters.
     */
    void writeParameters(CArrayXX const* p_tik, ostream& os) const;

  protected:
    /** pointer on the kernel */
    Kernel::IKernel const* p_kernel_;
    /** Array of the intermediate results dik */
    CArrayXX dik_;

    /** @brief Initialize the model before its first use.
     * This function is triggered when kernel is set.
     * In this interface, the @c initializeModel() method
     *  - set the number of samples and variables of the mixture model
     *  - call the derived class implemented method
     * @code
     *   initializeModelImpl()
     * @endcode
     * for initialization of the specific model parameters if needed.
     **/
    void initializeModel();
    /** compute the distance of the ith individual to the kth centroid
     *  \f[ d_{ik} = \|\phi(x_i)-m_k\| \f]
     *  using the kernel trick.
     **/
    void compute_dik(CArrayXX const* p_tik, CPointX const* p_tk);

  private:
    /** @brief Set the data set.
     *  method is re-implemented for debug purpose
     *  @param data the data set to set
     **/
    void setData(Array const& data) { this->p_dataij_ = 0;}
};


/* compute the intermediate quantities \f$ d_{ik}^m = \|\phi(x_i)-m^m_k\| \f$
 *  using the kernel trick.
 **/
template<class Derived>
void KernelGaussianBase<Derived>::compute_dik(CArrayXX const* p_tik, CPointX const* p_tk)
{
#ifdef STK_KernelS_DEBUG
  stk_cout << _T("Entering KernelGaussianBase::compute_dik\n");
  stk_cout << _T("dik_.cols() =") << dik_.cols() << _T("\n");
  stk_cout << _T("dik_.rows() =") << dik_.rows() << _T("\n");
#endif
  // vector with values wik=\sum_{j=1}^n k(x_i,x_j) t_{jk}/t_{.k}, for k=1,..,K
  CVectorX wi(dik_.rows());
  for (int k=dik_.beginCols(); k<dik_.endCols(); ++k)
  {
    for (int i= wi.begin(); i < wi.end(); ++i)
    {
      wi[i] = 0.;
      for (int j= wi.begin(); j < wi.end(); ++j)
      { wi[i] += p_kernel_->comp(i, j) * p_tik->elt(j,k)/p_tk->elt(k);}
    }
    // compute dik_ = k(i,i) - 2 * wik + \sum_{i=1}^n t_{ik} w_{ik}/t_{.k}
    Real scal =p_tik->col(k).dot(wi)/p_tk->elt(k);
    for (int i= wi.begin(); i<wi.end(); ++i)
    { dik_(i,k) = p_kernel_->diag(i) - 2. * wi[i] + scal  ;}
  }
#ifdef STK_KernelS_DEBUG
  stk_cout << _T("KernelGaussianBase::compute_dik done\n");
#endif
}

/* @brief Initialize the model before its first use.
 * This function is triggered when data set is set.
 * In this interface, the @c initializeModel() method
 *  - set the number of samples and variables of the mixture model
 *  - call the derived class implemented method
 * @code
 *   initializeModelImpl()
 * @endcode
 * for initialization of the specific model parameters if needed.
 **/
template<class Derived>
void KernelGaussianBase<Derived>::initializeModel()
{
  // set dimensions
  this->setNbSample(p_kernel_->nbSample());
  this->setNbVariable(p_kernel_->nbVariable());
  dik_.resize(nbSample(), nbCluster());
  // call specific model initialization stuff (not really necessary)
  this->asDerived().initializeModelImpl();
  //compute_dik();
}

/* This function is used in order to get the current values of the
 *  parameters in an array.
 *  @param[out] params the array with the parameters of the mixture.
 */
template<class Derived>
template<class ArrayParam>
inline void KernelGaussianBase<Derived>::getParameters(ArrayParam& param) const
{
  param.resize(this->nbCluster(), 2);
  for (int k= param.beginRows(); k < param.endRows(); ++k)
  {
    param(k, baseIdx  ) = param_.sigma2(k);
    param(k, baseIdx+1) = param_.dim(k);
  }
}

/* This function can be used to write summary of parameters to the output stream.
 *  @param os Stream where you want to write the summary of parameters.
 */
template<class Derived>
inline void KernelGaussianBase<Derived>::writeParameters(CArrayXX const* p_tik, ostream& os) const
{
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    os << _T("---> Component ") << k << _T("\n");
    os << _T("sigma2 = ") << param_.sigma2(k) << _T("\n");
    os << _T("dim = ")    << param_.dim(k)    << _T("\n");
  }
}

} // namespace STK

#endif /* STK_KernelGAUSSIANBASE_H */
