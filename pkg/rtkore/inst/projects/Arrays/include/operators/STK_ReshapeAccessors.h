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
 * created on: 17 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ReshapeAccessors.h
 *  @brief In this file we implement the DiagonalAccessor, DiagonalGetterAccessor,
 *   UpperTriangularizeAccessor,LowerTriangularizeAccessor, SymmetrizeAccessor,
 *   UpperSymmetrizeAccessor, LowerSymmetrizeAccessor classes.
 **/

#ifndef STK_RESHAPEACCESSORS_H
#define STK_RESHAPEACCESSORS_H


#define EGAL(arg1, arg2) ((arg1::structure_ == int(Arrays::arg2)))

namespace STK
{
// forward declaration
template< typename Lhs> class DiagonalizeAccessor;
template< typename Lhs> class DiagonalGetterAccessor;
template< typename Lhs> class UpperTriangularizeAccessor;
template< typename Lhs> class LowerTriangularizeAccessor;
template< typename Lhs> class SymmetrizeAccessor;
template< typename Lhs> class UpperSymmetrizeAccessor;
template< typename Lhs> class LowerSymmetrizeAccessor;

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for DiagonalizeAccessor operator
 */
template<typename Lhs>
struct Traits< DiagonalizeAccessor <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = (Lhs::structure_== int(Arrays::diagonal_) )
                 ? Lhs::sizeRows_ != UnknownSize ? Lhs::sizeRows_ : Lhs::sizeCols_
                 : Lhs::structure_== int(Arrays::point_) ? int(Lhs::sizeCols_) : int(Lhs::sizeRows_),
    sizeCols_  = (Lhs::structure_== int(Arrays::diagonal_) )
                 ? Lhs::sizeRows_ != UnknownSize ? Lhs::sizeRows_ : Lhs::sizeCols_
                 : Lhs::structure_== int(Arrays::point_) ?  Lhs::sizeCols_ : Lhs::sizeRows_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<DiagonalizeAccessor < Lhs> > Row;
  typedef ColOperator<DiagonalizeAccessor < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::TypeConst TypeConst;
};

} // namespace hidden


/** @ingroup Arrays
 *  @class DiagonalizeAccessor
  *
  * @brief Generic expression when a one dimensional vector/point/idagonal
  * expression is "diagonalized".
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * DiagonalizeAccessor operator.
  *
  * This class represents an expression where a DiagonalizeAccessor operator is
  * applied to a vector/point/diagonal array. It is the return type of the
  * diagonalizeize() operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalizeAccessor type explicitly.
  */
template< typename Lhs>
class DiagonalizeAccessor: public ArrayBase< DiagonalizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< DiagonalizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalizeAccessor<Lhs> >::TypeConst TypeConst;
    enum
    {
        structure_ = hidden::Traits< DiagonalizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalizeAccessor<Lhs> >::storage_,
        // this is safe as we can use DiagonalizeAccesor only on 1D container
        size_      = (sizeRows_ != UnknownSize) ? sizeRows_ : sizeCols_
    };
    /** Type of the Range for the rows */
    typedef TRange<size_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<size_> ColRange;
    /** Constructor */
    inline DiagonalizeAccessor(Lhs& lhs): Base(), lhs_(lhs)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Lhs);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().range();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().range();}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline TypeConst elt2Impl(int i, int j) const { return (lhs_.elt(i, j));}
    /** @return a constant reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline TypeConst elt1Impl(int i) const { return (lhs_.elt(i));}
    /** @return a constant reference on the element of the expression */
    inline TypeConst elt0Impl() const { return (lhs_.elt());}

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
      if (i != j) { STKOUT_OF_RANGE_2ARG(DiagonalizeAccessor::elt, i, j, i != j);}
#endif
      return (lhs_.elt(i));
    }
    /** @return a reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline Type& elt1Impl(int i) { return (lhs_.elt(i));}
    /** @return a reference on the element of the expression */
    inline Type& elt0Impl() { return (lhs_.elt());}

  protected:
    Lhs& lhs_;
};


namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for DiagonalGetterAccessor operator
 */
template<typename Lhs>
struct Traits< DiagonalGetterAccessor <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = ((Lhs::sizeRows_ < Lhs::sizeCols_)) ?  Lhs::sizeRows_ : Lhs::sizeCols_,
    sizeCols_  = sizeRows_,
    storage_   = Lhs::storage_,
    isValid_   = ( hidden::Traits<Lhs>::structure_==(int)Arrays::array2D_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::square_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::diagonal_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::lower_triangular_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::upper_triangular_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::symmetric_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::lower_symmetric_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::upper_symmetric_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::number_),
    use_       = (sizeRows_ != UnknownSize) ? Arrays::useLhsSize_ : Arrays::useLhsOtherSize_,
    size_      = (sizeRows_ != UnknownSize) ? sizeRows_ : sizeCols_
  };
  typedef RowOperator<DiagonalGetterAccessor < Lhs> > Row;
  typedef ColOperator<DiagonalGetterAccessor < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::TypeConst TypeConst;
};


} // namespace hidden


/** @ingroup Arrays
 *  @class DiagonalGetterAccessor
  *
  * @brief Generic expression when we want to get the diagonal of a
  * two-dimensional square expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * DiagonalGetterAccessor operator.
  *
  * This class represents an expression where a diagonal operator is applied to
  * an expression. It is the return type of the diagonal operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalGetterAccessor type explicitly.
  */
template< typename Lhs>
class DiagonalGetterAccessor: public ArrayBase< DiagonalGetterAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< DiagonalGetterAccessor< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalGetterAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalGetterAccessor<Lhs> >::TypeConst TypeConst;

    enum
    {
        structure_ = hidden::Traits< DiagonalGetterAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalGetterAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalGetterAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalGetterAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalGetterAccessor<Lhs> >::storage_,
        isValid_   = hidden::Traits< DiagonalGetterAccessor<Lhs> >::isValid_,
        use_       = hidden::Traits< DiagonalGetterAccessor<Lhs> >::use_,
        size_      = hidden::Traits< DiagonalGetterAccessor<Lhs> >::size_
    };
    typedef hidden::DiagonalRangeImpl<Lhs, size_, use_> RangeImpl;
    /** Type of the Range for the rows */
    typedef TRange<size_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<size_> ColRange;

    /** Constructor */
    inline DiagonalGetterAccessor( Lhs& lhs): Base(), lhs_(lhs)
    {
      STK_STATIC_ASSERT(isValid_,YOU_TRIED_CALLING_A_MATRIX_METHOD_ON_A_VECTOR);
#ifdef STK_BOUNDS_CHECK
      if (lhs.rows()!=lhs.cols())
        STKRUNTIME_ERROR_NO_ARG(DiagonalGetterAccessor,lhs.rows()!=lhs.cols());
#endif
    }
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return RangeImpl::rangeImpl(lhs_);}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return RangeImpl::rangeImpl(lhs_);}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline TypeConst elt2Impl(int i, int j) const { return (lhs_.elt(i, j));}
    /** @return a constant reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline TypeConst elt1Impl(int i) const { return (lhs_.elt(i,i));}
    /** @return a constant reference on the element of the expression */
    inline TypeConst elt0Impl() const { return (lhs_.elt());}

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline Type& elt2Impl(int i, int j) { return (lhs_.elt(i, j));}
    /** @return a reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline Type& elt1Impl(int i) { return (lhs_.elt(i,i));}
    /** @return a reference on the element of the expression */
    inline Type& elt0Impl() { return (lhs_.elt());}

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for UpperTriangularizeAccessor operator
 */
template<typename Lhs>
struct Traits< UpperTriangularizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::upper_triangular_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< UpperTriangularizeAccessor< Lhs> > Row;
  typedef ColOperator< UpperTriangularizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::TypeConst TypeConst;
};

} // namespace hidden

/** @ingroup Arrays
 *  @class UpperTriangularizeAccessor
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * UpperTriangularizeAccessor operator.
  *
  * This class represents an expression where an UpperTriangularizeOperator
  * operator is applied to an expression. It is the return type of the
  * upperTriangularize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name UpperTriangularizeAccessor type explicitly.
  */
template< typename Lhs>
class UpperTriangularizeAccessor: public ArrayBase< UpperTriangularizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< UpperTriangularizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< UpperTriangularizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< UpperTriangularizeAccessor<Lhs> >::TypeConst TypeConst;

    enum
    {
        structure_ = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline UpperTriangularizeAccessor( Lhs& lhs): Base(), lhs_(lhs) {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline TypeConst elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (j<i)
        STKRUNTIME_ERROR_2ARG(UpperTriangularizeAccessor::elt2Impl,i,j,use of the lower part);
#endif
      return (lhs_.elt(i, j));
    }

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j) {
#ifdef STK_BOUNDS_CHECK
      if (i>j)
        STKRUNTIME_ERROR_2ARG(UpperTriangularizeAccessor::elt2Impl,i,j,use of the lower part);
#endif
      return (lhs_.elt(i, j));
    }

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for LowerTriangularizeAccessor operator
 */
template<typename Lhs>
struct Traits< LowerTriangularizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::lower_triangular_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< LowerTriangularizeAccessor< Lhs> > Row;
  typedef ColOperator< LowerTriangularizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::TypeConst TypeConst;
};

} // namespace hidden

/** @ingroup Arrays
 *  @class LowerTriangularizeAccessor
  *
  * @brief Generic expression when we want to get the lower-part of a
  * two-dimensional expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * LowerTriangularizeAccessor operator.
  *
  * This class represents an expression where an LowerTriangularizeAccessor
  * operator is applied to an expression. It is the return type of the
  * lowerTriangularize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name LowerTriangularizeAccessor type explicitly.
  */
template< typename Lhs>
class LowerTriangularizeAccessor: public ArrayBase< LowerTriangularizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< LowerTriangularizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< LowerTriangularizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< LowerTriangularizeAccessor<Lhs> >::TypeConst TypeConst;

    enum
    {
        structure_ = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline LowerTriangularizeAccessor( Lhs& lhs): Base(), lhs_(lhs) {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline TypeConst elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (j>i)
        STKRUNTIME_ERROR_2ARG(LowerTriangularizeAccessor::elt2Impl,i,j,use of the upper part);
#endif
      return (lhs_.elt(i, j));
    }

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
      if (j>i)
        STKRUNTIME_ERROR_2ARG(LowerTriangularizeAccessor::elt2Impl,i,j,use of the upper part);
#endif
      return (lhs_.elt(i, j));
    }

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for SymmetrizeAccessor operator
 */
template<typename Lhs>
struct Traits< SymmetrizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::symmetric_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< SymmetrizeAccessor< Lhs> > Row;
  typedef ColOperator< SymmetrizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::TypeConst TypeConst;
};

} // namespace hidden

/** @ingroup Arrays
 *  @class SymmetrizeAccessor
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional symmetric expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * SymmetrizeAccessor operator.
  *
  * This class represents an expression where an SymmetrizeAccessor
  * operator is applied to an expression. It is the return type of the
  * symmetrize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name SymmetrizeAccessor type explicitly.
  */
template< typename Lhs>
class SymmetrizeAccessor: public ArrayBase< SymmetrizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< SymmetrizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< SymmetrizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< SymmetrizeAccessor<Lhs> >::TypeConst TypeConst;

    enum
    {
        structure_ = hidden::Traits< SymmetrizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< SymmetrizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< SymmetrizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< SymmetrizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< SymmetrizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline SymmetrizeAccessor( Lhs& lhs): Base(), lhs_(lhs)
    { STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Lhs);}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline TypeConst elt2Impl(int i, int j) const
    { return (lhs_.elt(i, j));}

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j) { return (lhs_.elt(i, j));}

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for UpperSymmetrizeAccessor operator
 */
template<typename Lhs>
struct Traits< UpperSymmetrizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::upper_symmetric_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< UpperSymmetrizeAccessor< Lhs> > Row;
  typedef ColOperator< UpperSymmetrizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::TypeConst TypeConst;
};

} // namespace hidden

/** @ingroup Arrays
 *  @class UpperSymmetrizeAccessor
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional symmetric expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * UpperSymmetrizeAccessor operator.
  *
  * This class represents an expression where an UpperSymmetrizeAccessor
  * operator is applied to an expression. It is the return type of the
  * upperSymmetrize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name UpperSymmetrizeAccessor type explicitly.
  */
template< typename Lhs>
class UpperSymmetrizeAccessor: public ArrayBase< UpperSymmetrizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< UpperSymmetrizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::TypeConst TypeConst;

    enum
    {
        structure_ = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline UpperSymmetrizeAccessor( Lhs& lhs): Base(), lhs_(lhs)
    { STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Lhs);}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline TypeConst elt2Impl(int i, int j) const
    { return ((j<i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j)
    { return ((j<i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for LowerSymmetrizeAccessor operator
 */
template<typename Lhs>
struct Traits< LowerSymmetrizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::lower_symmetric_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< LowerSymmetrizeAccessor< Lhs> > Row;
  typedef ColOperator< LowerSymmetrizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::TypeConst TypeConst;
};

} // end namespace hidden


/** @ingroup Arrays
 *  @class LowerSymmetrizeAccessor
  *
  * @brief Generic expression when we want to get the lower-part of a
  * two-dimensional symmetric expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * LowerSymmetrizeAccessor operator.
  *
  * This class represents an expression where an LowerTriangularizeAccessor
  * operator is applied to an expression. It is the return type of the
  * lowerSymmetrize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name LowerSymmetrizeAccessor type explicitly.
  */
template< typename Lhs>
class LowerSymmetrizeAccessor: public ArrayBase< LowerSymmetrizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< LowerSymmetrizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::TypeConst TypeConst;

    enum
    {
        structure_ = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline LowerSymmetrizeAccessor( Lhs& lhs): Base(), lhs_(lhs)
    { STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Lhs);}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline TypeConst elt2Impl(int i, int j) const
    { return ((j>i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j)
    { return ((j>i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

  protected:
    Lhs& lhs_;
};

} // namespace STK


#undef EGAL

#endif /* STK_RESHAPEACCESSORS_H */
