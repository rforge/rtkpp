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
 * Project:  stkpp::Arrays
 * created on: 13 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ArrayBase.h
 *  @brief In this file we define the base class for Arrays. Derived ArrayBase can be
 *  a lhs
 **/

#ifndef STK_ARRAYBASE_H
#define STK_ARRAYBASE_H

#include "STK_ExprBase.h"
#include "STatistiK/include/STK_Law_IUnivLaw.h"


namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  Utility class allowing to know if in an assignment the dimensions are
 *  correct
 **/
template<class Derived, int Structure_, int RhsStucture_>
struct CheckAssign;

/** @ingroup hidden
 *  Utility class allowing to know if in an assignment the destination must
 *   be resized or shifted
 **/
template<class Derived, int Structure_>
struct CheckShift;

/** @ingroup hidden
 *  Specialization for general array2D_
 **/
template<class Derived,int RhsStructure_>
struct CheckAssign<Derived, Arrays::array2D_, RhsStructure_>
{
  // all range are authorized for array2D_
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
/** @ingroup hidden
 *  Specialization for square_
 **/
template<class Derived, int RhsStructure_>
struct CheckAssign<Derived, Arrays::square_, RhsStructure_>
{
  // same ranges for square_ arrays
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return I==J;}
};

/** @ingroup hidden
 *  Specialization for upper_triangular_
 **/
template<class Derived>
struct CheckAssign<Derived, Arrays::upper_triangular_, Arrays::upper_triangular_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
/** @ingroup hidden
 *  Specialization for lower_triangular_
 **/
template<class Derived>
struct CheckAssign<Derived, Arrays::lower_triangular_, Arrays::lower_triangular_>
{
  // all range are authorized for lower_triangular_
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};

/** @ingroup hidden
 *  Specialization for diagonal_
 **/
template<class Derived>
struct CheckAssign<Derived, Arrays::diagonal_, Arrays::diagonal_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
template<class Derived>
struct CheckAssign<Derived, Arrays::diagonal_, Arrays::vector_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
template<class Derived>
struct CheckAssign<Derived, Arrays::diagonal_, Arrays::point_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};

/** @ingroup hidden
 *  Specialization for vector_
 **/
template<class Derived>
struct CheckAssign<Derived, Arrays::vector_, Arrays::diagonal_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
template<class Derived>
struct CheckAssign<Derived, Arrays::vector_, Arrays::vector_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
template<class Derived>
struct CheckAssign<Derived, Arrays::vector_, Arrays::point_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
/** @ingroup hidden
 *  Specialization for point_
 **/
template<class Derived>
struct CheckAssign<Derived, Arrays::point_, Arrays::diagonal_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
template<class Derived>
struct CheckAssign<Derived, Arrays::point_, Arrays::vector_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
template<class Derived>
struct CheckAssign<Derived, Arrays::point_, Arrays::point_>
{
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
};
// for number_
template<class Derived>
struct CheckAssign<Derived, Arrays::number_, Arrays::number_>
{
  // same range only for diagonal_ arrays
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return (I.size() == 1 && J.size() == 1);}
};

/** @ingroup hidden
 *  Specialization for general array2D_
 **/
template<class Derived>
struct CheckShift<Derived, Arrays::array2D_>
{
  // all range are authorized for array2D_
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
  // check if resize is necessary
  static bool resize(Derived const& array, Range const& I, Range const& J)
  { return (array.rows() != I || array.cols() != J);}
  // check if shift is necessary
  static bool shift(Derived const& array, int beginRow, int beginCol)
  { return (array.beginRows() != beginRow || array.beginCols() != beginCol);}
  // check if resize is necessary
  static bool resize(Derived const& array, Range const& I)
  { return (array.rows() != I || array.cols() != I);}
  // check if shift is necessary
  static bool shift(Derived const& array, int begin)
  { return (array.beginRows() != begin || array.beginCols() != begin);}
};
/** @ingroup hidden
 *  Specialization for upper_triangular_
 **/
template<class Derived>
struct CheckShift<Derived, Arrays::upper_triangular_>
{
  // all range are authorized for upper_triangular_
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
  static bool resize(Derived const& array, Range const& I, Range const& J)
  { return (array.rows() != I || array.cols() != J);}
  static bool shift(Derived const& array, int beginRow, int beginCol)
  { return (array.beginRows() != beginRow || array.beginCols() != beginCol);}
  static bool resize(Derived const& array, Range const& I)
  { return (array.rows() != I || array.cols() != I);}
  static bool shift(Derived const& array, int begin)
  { return (array.beginRows() != begin || array.beginCols() != begin);}
};
/** @ingroup hidden
 *  Specialization for lower_triangular_
 **/
template<class Derived>
struct CheckShift<Derived, Arrays::lower_triangular_>
{
  // all range are authorized for lower_triangular_
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return true;}
  static bool resize(Derived const& array, Range const& I, Range const& J)
  { return (array.rows() != I || array.cols() != J);}
  static bool shift(Derived const& array, int beginRow, int beginCol)
  { return (array.beginRows() != beginRow || array.beginCols() != beginCol);}
  static bool resize(Derived const& array, Range const& I)
  { return (array.rows() != I || array.cols() != I);}
  static bool shift(Derived const& array, int begin)
  { return (array.beginRows() != begin || array.beginCols() != begin);}
};

/** @ingroup hidden
 *  Specialization for square_
 **/
template<class Derived>
struct CheckShift<Derived, Arrays::square_>
{
  // same range only for square_ arrays
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return I==J;}
  static bool resize(Derived const& array, Range const& I, Range const& J)
  { return (array.rows() != I || array.cols() != J);}
  static bool shift(Derived const& array, int beginRow, int beginCol)
  { return (array.beginRows() != beginRow || array.beginCols() != beginCol);}
  static bool resize(Derived const& array, Range const& I)
  { return (array.range() != I);}
  static bool shift(Derived const& array, int begin)
  { return (array.beginRows() != begin || array.beginCols() != begin);}
};

/** @ingroup hidden
 *  Specialization for diagonal_
 **/
template<class Derived>
struct CheckShift<Derived, Arrays::diagonal_>
{
  // same range only for diagonal_ arrays
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return I==J;}
  static bool resize(Derived const& array, Range const& I, Range const& J)
  { return (array.rows() != I || array.cols() != J);}
  static bool shift(Derived const& array, int beginRow, int beginCol)
  { return (array.beginRows() != beginRow || array.beginCols() != beginCol);}
  static bool resize(Derived const& array, Range const& I)
  { return (array.range() != I);}
  static bool shift(Derived const& array, int begin)
  { return (array.beginRows() != begin || array.beginCols() != begin);}
};

// for vectors
template<class Derived>
struct CheckShift<Derived, Arrays::vector_>
{
  // same range only for vector_ arrays
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return J.size() == 1;}
  static bool resize(Derived const& array, Range const& I)
  { return (array.range() != I);}
  static bool shift(Derived const& array, int begin)
  { return (array.begin() != begin);}
};

// for point
template<class Derived>
struct CheckShift<Derived, Arrays::point_>
{
  // same range only for diagonal_ arrays
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return I.size() == 1;}
  static bool resize(Derived const& array, Range const& I)
  { return (array.range() == I) ? false : true;}
  static bool shift(Derived const& array, int begin)
  { return (array.begin() == begin) ? false : true;}
};
// for point
template<class Derived>
struct CheckShift<Derived, Arrays::number_>
{
  // same range only for diagonal_ arrays
  static bool isAllowed(Derived const& array, Range const& I, Range const& J)
  { return (I.size() == 1 && J.size() == 1);}
  static bool resize(Derived const& array, Range const& I)
  { return (array.range() == I) ? false : true;}
  static bool shift(Derived const& array, int begin)
  { return (array.begin() == begin) ? false : true;}
};

} // namespace hidden
/** @ingroup Arrays
 *  @brief base class for template arrays.
 *
 * This class is the base that is inherited by all containers storing
 * values (matrix, vector, point). Expressions are not arrays. Any derived
 * class can be a lhs in an expression.
 *
 * The common API for these objects is contained in this class.
 *
 *  @tparam Derived is the derived type, e.g., a matrix, vector, point type or
 *  an expression.
 **/
template< class Derived>
class ArrayBase :  public ExprBase<Derived>
{
  public:
    typedef ExprBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::ReturnType ReturnType;

    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;
    typedef typename hidden::Traits<Derived>::SubRow SubRow;
    typedef typename hidden::Traits<Derived>::SubCol SubCol;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;
    typedef typename hidden::Traits<Derived>::SubArray SubArray;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };

  protected:
    /** Default constructor. Default values are cols=(1:0) and rows=(1:0). */
    ArrayBase() : Base() {}
    /** destructor */
    ~ArrayBase() {}

  public:
    /** Apply the Visitor @c visitor to the whole coefficients of the array.
      * The template parameter @c Visitor is the type of the visitor and provides
      * the following interface:
      * @code
      * struct MyVisitor {
      *   // called for all  coefficients
      *   void operator() (Type& value);
      * };
      * @endcode
      *
      * @note visitors offer automatic unrolling for small fixed size matrix.
      *
      * @sa setValue, setOnes(), setZeros()
      */
    template<typename Visitor>
    void apply(Visitor& visitor);
    /** set random values to this using a uniform law. @sa apply(), randGauss() */
    Derived& randUnif();
    /** set random values to this using a standard gaussian law.
     *  @sa randUnif(), rand(Law::IUnivLaw<Type> const& law), apply()
     **/
    Derived& randGauss();
    /** set random values to this using a law given by the user. @sa randGauss(), randUnif(), apply() */
    Derived& rand( Law::IUnivLaw<Type> const& law);
    /** set one to this using a Visitor. @sa apply(), setValue(), setZeros() */
    Derived& setOnes();
    /** set zero to this using a Visitor. @sa apply(), setOnes(), setValue()*/
    Derived& setZeros();
    /** set one to this using a Visitor. @sa apply(), setValue(), zeros() */
    Derived& ones();
    /** set zero to this using a Visitor. @sa apply(), ones(), setValue()*/
    Derived& zeros();

    /** set a value to this container. @sa apply(), setOnes(), setZeros()
     *  @param value the value to set
     **/
    Derived& setValue(Type const& value);
    /** @return a copy of @c rhs inside @c this object.
     *  If the ranges of @c this and @c rhs are not exactly the same, the assign
     *  method will call the resize method on this.
     *
     *  @note If @c this is a reference, it cannot be resized and thus an
     *  exception will be thrown.
     **/
    template<class Rhs> Derived& assign(ExprBase<Rhs> const& rhs);

    /** @return the matrix or vector obtained by setting this constant*/
    Derived& operator=( Type const& rhs) { return setValue(rhs);}
    /** @return the matrix or vector obtained by evaluating this expression */
    Derived& operator=( Derived const& rhs) { return assign(rhs);}
    /** @return the matrix or vector obtained by evaluating this expression */
    template<typename Rhs>
    inline Derived& operator=( ExprBase<Rhs> const& rhs)
    { return assign(rhs.asDerived());}
    /** Adding a Rhs to this. */
    template<typename Rhs>
    inline Derived& operator+=( ExprBase<Rhs> const& other);
    /** subtract a Rhs to this. */
    template<typename Rhs>
    inline Derived& operator-=( ExprBase<Rhs> const& other);
    /** divide this by Rhs. */
    template<typename Rhs>
    inline Derived& operator/=( ExprBase<Rhs> const& other);
    /** multiply this by Rhs. */
    template<typename Rhs>
    inline Derived& operator*=( ExprBase<Rhs> const& other);
    /** Adding a constant to this. */
    inline Derived& operator+=( Type const& other);
    /** Substract a constant to this. */
    inline Derived& operator-=( Type const& other);
    /** product of this by a constant. */
    inline Derived& operator*=( Type const& other);
    /** dividing this by a constant. */
    inline Derived& operator/=( Type const& other);

    /** @return the jth column of this */
    inline Col col(int j) const { return this->asDerived().colImpl(j);}
    /** @return the jth column of this in the range I*/
    inline SubCol col(Range const& I, int j) const
    {
      STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Derived)
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > I.begin())
      { STKOUT_OF_RANGE_2ARG(ArrayBase::col, I, j, beginRows() > I.begin());}
      if (this->endRows() < I.end())
      { STKOUT_OF_RANGE_2ARG(ArrayBase::col, I, j, endRows() < I.end());}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(ArrayBase::col, I, j, beginCols() > j);}
      if (this->endCols() <= j)
      { STKOUT_OF_RANGE_2ARG(ArrayBase::col, I, j, endCols() <= j);}
#endif
      return this->asDerived().colImpl(I, j);
    }
    /** @return the sub array with the column in the range J */
    inline SubArray col(Range const& J) const
    {
      STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Derived)
#ifdef STK_BOUNDS_CHECK
      if (this->beginCols() > J.begin())
      { STKOUT_OF_RANGE_1ARG(ArrayBase::col, J, beginCols() > J.begin());}
      if (this->endCols() < J.end())
      { STKOUT_OF_RANGE_1ARG(ArrayBase::col, J, endCols() < J.end());}
#endif
      return this->asDerived().colImpl(J);
    }
    /** @return the ith row of this */
    inline Row row(int i) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_1ARG(ArrayBase::row, i, beginRows() > i);}
      if (this->endRows() <= i)
      { STKOUT_OF_RANGE_1ARG(ArrayBase::row, i, endRows() <= i);}
#endif
      return this->asDerived().rowImpl(i);
    }
    /** @return the ith row of this in the range J */
    inline SubRow row(int i, Range const& J) const
    {
      STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Derived)
      return this->asDerived().rowImpl(i,J);
    }
    /** @return the sub array with the rows in the range J */
    inline SubArray row(Range const& I) const
    {
      STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Derived)
      return this->asDerived().rowImpl(I);
    }
    /** @return the sub-array (I,J)*/
    inline SubArray sub(Range const& I, Range const& J) const
    {
      STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Derived)
      return this->asDerived().subImpl(I, J);
    }
    /** @return the sub-vector in the range I*/
    inline SubVector sub(Range const& I) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived)
      return this->asDerived().subImpl(I);
    }
    /** @return the element (i,j) of the 2D container.
     *  @param i,j indexes of the row and column
     **/
    inline Type& elt(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_2ARG(Type& ArrayBase::elt, i, j, beginRows() > i);}
      if (this->endRows() <= i)
      { STKOUT_OF_RANGE_2ARG(Type& ArrayBase::elt, i, j, endRows() <= i);}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(Type& ArrayBase::elt, i, j, beginCols() > j);}
      if (this->endCols() <= j)
      { STKOUT_OF_RANGE_2ARG(Type& ArrayBase::elt, i, j, endCols() <= j);}
#endif
      return this->asDerived().elt2Impl(i,j);
    }
    /** @return a constant reference on element (i,j) of the 2D container
     *  @param i, j indexes of the row and of the column
     **/
    inline Type const& elt(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_2ARG(Type const& ArrayBase::elt, i, j, beginRows() > i);}
      if (this->endRows() <= i)
      { STKOUT_OF_RANGE_2ARG(Type const& ArrayBase::elt, i, j, endRows() <= i);}
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_2ARG(Type const& ArrayBase::elt, i, j, beginCols() > j);}
      if (this->endCols() <= j)
      { STKOUT_OF_RANGE_2ARG(Type const& ArrayBase::elt, i, j, endCols() <= j);}
#endif
      return this->asDerived().elt2Impl(i,j);
    }
    /** @return a reference on the ith element
     *  @param i index of the ith element
     **/
    inline Type& elt(int i)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived)
#ifdef STK_BOUNDS_CHECK
      if (this->begin() > i)
      { STKOUT_OF_RANGE_1ARG(Type& ArrayBase::elt, i, begin() > i);}
      if (this->end() <= i)
      { STKOUT_OF_RANGE_1ARG(Type& ArrayBase::elt, i, end() <= i);}
#endif
      return this->asDerived().elt1Impl(i);
    }
    /** @return the constant ith element
     *  @param i index of the ith element
     **/
    inline Type const& elt(int i) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived)
#ifdef STK_BOUNDS_CHECK
      if (this->begin() > i)
      { STKOUT_OF_RANGE_1ARG(Type const& ArrayBase::elt, i, begin() > i);}
      if (this->end() <= i)
      { STKOUT_OF_RANGE_1ARG(Type const& ArrayBase::elt, i, end() <= i);}
#endif
      return this->asDerived().elt1Impl(i);
    }
    /** @return a reference on the number */
    inline Type& elt()
    {
      STK_STATIC_ASSERT_ZERO_DIMENSION_ONLY(Derived)
      return this->asDerived().elt0Impl();
    }
    /** @return a constant reference on the number */
    inline Type const& elt() const
    {
      STK_STATIC_ASSERT_ZERO_DIMENSION_ONLY(Derived)
      return this->asDerived().elt0Impl();
    }
    // overloaded operators
    /** @return a reference on the ith element
     *  @param i index of the ith element
     **/
    inline Type& operator[](int i) { return elt(i);}
    /** @return the ith element
     *  @param i index of the ith element
     **/
    inline Type const& operator[](int i) const { return elt(i);}
    /** @return the ith element
     *  @param I range to get
     **/
    inline SubVector operator[](Range const& I) const { return sub(I);}
    /** @return a reference on the element (i,j) of the 2D container.
     *  @param i, j indexes of the row and of the column
     **/
    inline Type& operator()(int i, int j) { return elt(i,j);}
    /** @return a constant reference on the element (i,j) of the 2D container.
     *  @param i,j indexes of the row and column
     **/
    inline Type const& operator()(int i, int j) const { return elt(i,j);}
    /** @return a constant reference on the number */
    inline Type const& operator()() const { return elt();}
    /** @return the number */
    inline Type& operator()() { return elt();}
    /** @param I range of the index of the rows
     *  @param j index of the column
     *  @return a Vertical container containing the column @c j of this
     *  in the range @c I
     **/
    inline SubCol operator()(Range const& I, int j) const { return col(I, j);}
    /** @param i index of the row
     *  @param J range of the columns
     *  @return an Horizontal container containing the row @c i of this
     *  in the range @c J
     **/
    inline SubRow operator()(int i, Range const& J) const { return row(i, J);}
    /** @param I,J range of the rows and of the columns
     *  @return a 2D container containing this in the range @c I, @c J
     **/
    inline SubArray operator()(Range const& I, Range const& J) const { return sub(I, J);}
    //
    /** @return the first element */
    inline Type& front()
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived)
      return elt(this->begin());
    }
    /** @return the last element */
    inline Type& back()
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived)
      return elt(this->lastIdx());
    }
    /** @return the first element */
    inline Type const& front() const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived)
      return elt(this->begin());
    }
    /** @return the last element */
    inline Type const& back() const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived)
      return elt(this->lastIdx());
    }
    /** overwrite @c this with @c src.
     *  @note this method does not take care of the possibility of overlapping
     *  @param src the array to copy
     **/
    template<class OtherDerived>
    Derived& copy( ExprBase<OtherDerived> const& rhs);
    /** Convenient operator to set the coefficients of a matrix.
     *
     * The coefficients must be provided in the row/column order and exactly
     * match the size of the matrix. Otherwise an exception is throwed.
     **/
    ArrayInitializer<Derived> operator<<(Type const& s);

    /** \sa operator<<(Type const&) */
    template<typename Rhs>
    ArrayInitializer<Derived> operator<<(ArrayBase<Rhs> const& other);
};

} // namespace STK

#endif /* STK_ARRAYBASE_H_ */
