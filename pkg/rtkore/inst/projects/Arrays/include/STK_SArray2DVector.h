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
 * Project: stkpp::Arrays
 * Purpose:  Define the SArray2DVector class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_SArray2DVector.h
  * @brief A SArray2DVector is a one dimensional horizontal container
  *
  * An SArray2DVector is an implementation of the interface IArray2D.
  * It's a one dimensional horizontal container.
 **/

#ifndef STK_SARRAY2DVECTOR_H
#define STK_SARRAY2DVECTOR_H

#include "STK_SArray1D.h"

#include "STK_IArray2D.h"
#include "STK_IArray2DSlicers.h"
#include "STK_IArray2DModifiers.h"

namespace STK
{

template<typename> class SArray2DPoint;
template<typename> class SArray2DVector;


/** @ingroup Arrays
  * @brief final class for a Real vertical container.
  *
  * A Vector is a column oriented 1D container of Real.
  **/
typedef SArray2DVector<Real>   SVectorX;
typedef SArray2DVector<double> SVectorXd;
typedef SArray2DVector<int>    SVectorXi;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for the SArray2DVector class.
 **/
template<class Type_>
struct Traits< SArray2DVector<Type_> >
{
  typedef SArray2DPoint<Type_>  Row;
  typedef SArray2DVector<Type_> Col;
  typedef SArray2DPoint<Type_>  SubRow;
  typedef SArray2DVector<Type_> SubCol;
  typedef SArray2DVector<Type_> SubArray;
  typedef SArray2DVector<Type_> SubVector;

  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type const& TypeConst;

  enum
  {
    structure_ = Arrays::vector_,
    orient_    = Arrays::by_col_,
    sizeCols_  = 1,
    sizeRows_  = UnknownSize,
    size_      = UnknownSize,
    storage_   = Arrays::sparse_
  };
  typedef SArray1D<Type, UnknownSize, UnknownSize> ColVector;
};

} // namespace hidden

/** @ingroup Arrays
 *  @brief template one dimensional horizontal Array.
 *
 *  An SArray2DVector is a Vertical container of a single column.
 *
 *  By default the index of the first element is 0 but this can be
 *  modified using the appropriate constructor or using the method @c shift.
 *
 *  @sa SArray2DPoint
 **/
template<class Type_>
class SArray2DVector: public IArray2D< SArray2DVector<Type_> >
{
  public:
    typedef IArray2D< SArray2DVector > Base;
    typedef ArrayBase< SArray2DVector > LowBase;

    typedef typename hidden::Traits<SArray2DVector >::Row Row;
    typedef typename hidden::Traits<SArray2DVector >::Col Col;
    typedef typename hidden::Traits<SArray2DVector >::SubRow SubRow;
    typedef typename hidden::Traits<SArray2DVector >::SubCol SubCol;
    typedef typename hidden::Traits<SArray2DVector >::SubVector SubVector;
    typedef typename hidden::Traits<SArray2DVector >::SubArray SubArray;

    typedef typename hidden::Traits<SArray2DVector >::Type Type;
    typedef typename hidden::Traits<SArray2DVector >::TypeConst TypeConst;

    enum
    {
      structure_ = hidden::Traits< SArray2DVector >::structure_,
      orient_    = hidden::Traits< SArray2DVector >::orient_,
      sizeCols_  = hidden::Traits< SArray2DVector >::sizeCols_,
      sizeRows_  = hidden::Traits< SArray2DVector >::sizeRows_,
      size_      = hidden::Traits< SArray2DVector >::size_,
      storage_   = hidden::Traits< SArray2DVector >::storage_
    };
    using Base::elt;
    /** Default constructor */
    SArray2DVector(): Base( Range(), Range(1)) {}
    /** constructor with specified range.
     *  @param I range of the container
     **/
    SArray2DVector( Range const& I) :Base(I, Range(1)) {}
    /** constructor with specified range, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    SArray2DVector( Range const& I, Type const& v): Base(I, Range(1))
    { LowBase::setValue(v);}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if this is a wrapper of T
     **/
    SArray2DVector( const SArray2DVector &T, bool ref =false)
                 : Base(T, ref) {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the columns range to wrap
     **/
    SArray2DVector( const SArray2DVector& T, Range const& I)
                 : Base(T, I, T.cols())
    {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the range of the data to wrap
     *  @param col the index of the column to wrap
     **/
    template<class OtherArray>
    SArray2DVector( IArray2D<OtherArray> const& T, Range const& I, int col)
                 : Base(T, I, Range(col, 1))
    {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    SArray2DVector( ExprBase<OtherDerived> const& T): Base( Range(), Range(1))
    { LowBase::operator=(T);}
    /** constructor by reference, ref_=1.
     *  @param p_data a pointer on the data to wrap
     *  @param I the range of the data to wrap
     *  @param col the index of the column to wrap
     **/
     SArray2DVector( Type** p_data, Range const& I, int col)
                 : Base(p_data, I, Range(col, 1))
    {}
    /** destructor. */
    ~SArray2DVector() {}
    /** @return a constant reference on the ith element
     *  @param i index of the element (const)
     **/
    Type const& elt1Impl( int i) const { return this->elt(i, this->beginCols());}
    /** New first index for the object.
     *  @param rbeg the index of the first row to set
     **/
    void shift1D( int rbeg) { Base::shift(rbeg, this->beginCols());}
    /**  Resize the container.
     *  @param I the range to set to the container
     **/
    SArray2DVector<Type>& resize1D( Range const& I)
    { Base::resize(I, this->cols()); return *this;}
    /** Set value at position i
     *  @param i,v position and value to set
     **/
    void setValue1D( int i, TypeConst v)
    { this->setValue(i, this->beginCols(), v);}

    /** Add n elements to the container.
     *  @param n number of elements to add
     **/
    void pushBack( int n=1) { Base::pushBackRows(n);}
    /** Delete n elements at the pos index to the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
    void erase( int pos, int const& n=1) { Base::eraseRows(pos, n);}
    /** Insert n elements at the position pos of the container. The bound
     *  end_ should be modified at the very end of the insertion as pos
     *  can be a reference to it.
     *  @param pos index where to insert elements
     *  @param n number of elements to insert (default 1)
     **/
    void insertElt(int pos, int const& n =1)
    { Base::insertRows(pos, n);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    SArray2DVector& operator=(ExprBase<Rhs> const& T) { return LowBase::operator=(T);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    SArray2DVector& operator=(SArray2DVector const& T) { return LowBase::assign(T);}
    /** set the container to a constant value.
     *  @param v the value to set
     **/
    SArray2DVector& operator=(Type const& v) { return LowBase::setValue(v);}
};

} // namespace STK

#endif // STK_SARRAY2DVECTOR_H
