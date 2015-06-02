/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015  Serge Iovleff

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
 * Project:  stkpp::STatistiK::Law
 * created on: 23 janv. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_UniformDiscrete.h
 *  @brief In this file we implement the uniform (discrete) law.
 **/

#ifndef STK_LAW_UNIFORMDISCRETE_H
#define STK_LAW_UNIFORMDISCRETE_H

#include "STK_Law_IUnivLaw.h"
#include "STKernel/include/STK_Integer.h"

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief class for the Uniform law distribution.
 *  In probability theory and statistics, the <emWdiscrete uniform distributionw/em>
 *  is a probability distribution whereby a finite number of values are equally
 *  likely to be observed; every one of @e n values has equal probability <em>1/n</em>.
 *  Another way of saying "discrete uniform distribution" would be "a known,
 *  finite number of outcomes equally likely to happen".
 *
 *  The probability density function of the discrete uniform distribution is:
 *  \f[
 *    f(x; a, b) = \frac{1}{b-a+1} 1_{ a \leq x \leq b}, \quad a,b,x\in\mathbb{N}.
 *  \f]
**/
class UniformDiscrete : public IUnivLaw<int>
{
  public:
    typedef IUnivLaw<int> Base;
    /** constructor.
     *  @param a,b the lower and upper bounds
     **/
    UniformDiscrete( int a, int b);
    /** copy constructor.
     *  @param law the law to copy
     **/
    UniformDiscrete( UniformDiscrete const& law);
    /** destructor. */
	  virtual ~UniformDiscrete();
    /** @return the lower bound */
    inline int const& a() const { return a_;}
    /** @return the upper bound */
    inline int const& b() const { return b_;}
    /** @return the value b-a+1 */
    inline Real const& n() const { return n_;}
    /** @param a set the lower bound */
    inline void setA(int a) { a_ =a; n_ = b_-a_+1;}
    /** @param b set the upper bound */
    inline void setB(int b){ b_ =b; n_ = b_ - a_ + 1;}

    /** Generate a pseudo Uniform random variate. */
    virtual int rand() const;
    /** Give the value of the pdf at x.
     *  @param x a real value
     **/
    virtual Real pdf( int const& x) const;
    /** Give the value of the log-pdf at x.
     *  @param x a real value
     **/
    virtual Real lpdf( int const& x) const;
    /** The cumulative distribution function is
     * \f[
     *  F(t; a,b)= \frac{t - a}{b-a}
     * \f]
     *  @param t a real value
     **/
    virtual Real cdf( Real const& t) const;
    /** The inverse cumulative distribution function is
     * \f[
     * F^{-1}(p; \lambda) = p (b-a) + a.
     * \f]
     *  @param p a probability
     **/
    virtual int icdf( Real const& p) const;

    /** Generate a pseudo Uniform random variate.
     *  @param a,b the lower and upper bounds
     **/
    static int rand( int a, int b);
    /** Give the value of the pdf at x.
     *  @param x a real value
     *  @param a,b the lower and upper bounds
     **/
    static Real pdf( Real const& x, int a, int b);
    /** Give the value of the log-pdf at x.
     *  @param p a probablility
     *  @param a,b the lower and upper bounds
     **/
    static Real lpdf( Real const& p, int a, int b);
    /** Give the value of the cdf at t.
     *  @param t a real value
     *  @param a,b the lower and upper bounds
     **/
    static Real cdf( Real const& t, int a, int b);
    /** Give the value of the quantile at @e p.
     *  @param p a probability
     *  @param a,b the lower and upper bounds
     **/
    static int icdf( Real const& p, int a, int b);

  protected:
    /** The lower bound. */
    int a_;
    /** The upper bound. */
    int b_;

  private:
    Real n_;
};

} // namespace Law

} // namespace STK

#endif /*STK_LAW_UNIFORMDISCRETE_H*/
