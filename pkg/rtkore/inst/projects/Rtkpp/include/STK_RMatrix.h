/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2015  Serge Iovleff, University Lille 1, Inria

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  rtkpp
 * created on: 25 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_RMatrix.h
 *  @brief In this file we define the wrapper of the the Rcpp matrices.
 **/


#ifndef STK_RMATRIX_H
#define STK_RMATRIX_H

#include <Arrays/include/STK_ExprBaseVisitor.h>
#include <Arrays/include/STK_ExprBaseDot.h>
#include <Arrays/include/STK_ExprBaseProduct.h>

#include <Arrays/include/STK_ArrayBaseApplier.h>
#include <Arrays/include/STK_ArrayBaseAssign.h>
#include <Arrays/include/STK_ArrayBaseInitializer.h>


namespace STK
{

// forward declaration
template <typename Type_> class RMatrix;
template <typename Type_> class RowRMatrix;
template <typename Type_> class ColRMatrix;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of Traits class for STK::RMatrix class
 **/
template<typename Type_>
struct Traits< RMatrix<Type_> >
{
  private:
    class Void {};
  public:
    typedef RowRMatrix< Type_ > Row;
    typedef ColRMatrix< Type_ > Col;

    typedef RowOperator< RMatrix<Type_> > SubRow;
    typedef ColOperator< RMatrix<Type_> > SubCol;

    typedef Void SubVector;
    typedef SubOperator< RMatrix<Type_>, UnknownSize, UnknownSize > SubArray;
    typedef Void Number;

    typedef Type_ Type;
    typedef Type const& ReturnType;
    typedef Type const& ConstReturnType;

   enum
    {
      structure_ = Arrays::array2D_,
      orient_    = Arrays::by_col_,
      sizeRows_  = UnknownSize,
      sizeCols_  = UnknownSize,
      storage_   = Arrays::dense_
    };
};

/** @ingroup hidden
 *  @brief Specialization of Traits class for STK::RowRMatrix class
 **/
template<typename Type_>
struct Traits< RowRMatrix<Type_> >
{
  private:
    class Void {};
  public:
    typedef RowOperator< RowRMatrix<Type_> > Row;
    typedef RowOperator< RowRMatrix<Type_> > SubRow;
    typedef ColOperator< RowRMatrix<Type_> > Col;
    typedef ColOperator< RowRMatrix<Type_> > SubCol;
    typedef Void SubVector;
    typedef SubOperator< RowRMatrix<Type_>, UnknownSize, UnknownSize > SubArray;
    typedef Void Number;

    typedef Type_                Type;
    typedef typename RemoveConst<Type>::Type const& ConstReturnType;

   enum
    {
      structure_ = Arrays::array2D_,
      orient_    = Arrays::by_col_,
      sizeRows_  = 1,
      sizeCols_  = UnknownSize,
      storage_   = Arrays::dense_
    };
};

/** @ingroup hidden
 *  @brief Specialization of Traits class for STK::ColRMatrix class
 **/
template<typename Type_>
struct Traits< ColRMatrix<Type_> >
{
  private:
    class Void {};
  public:
    typedef RowOperator< ColRMatrix<Type_> > Row;
    typedef RowOperator< ColRMatrix<Type_> > SubRow;
    typedef ColOperator< ColRMatrix<Type_> > Col;
    typedef ColOperator< ColRMatrix<Type_> > SubCol;
    typedef Void SubVector;
    typedef SubOperator< ColRMatrix<Type_>, UnknownSize, UnknownSize > SubArray;
    typedef Void Number;

    typedef Type_                Type;
    typedef typename RemoveConst<Type>::Type const& ConstReturnType;

   enum
    {
      structure_ = Arrays::array2D_,
      orient_    = Arrays::by_col_,
      sizeRows_  = UnknownSize,
      sizeCols_  = 1,
      storage_   = Arrays::dense_
    };
};

} // namespace hidden

template <typename Type_>
class RMatrix: public ArrayBase< RMatrix<Type_> >, public TRef<1>
{
  public:
    typedef typename hidden::Traits<RMatrix<Type_> >::Type Type;
    typedef typename hidden::Traits<RMatrix<Type_> >::ConstReturnType ConstReturnType;
    enum
    {
      structure_ = hidden::Traits<RMatrix<Type_> >::structure_,
      orient_    = hidden::Traits<RMatrix<Type_> >::orient_,
      sizeRows_  = hidden::Traits<RMatrix<Type_> >::sizeRows_,
      sizeCols_  = hidden::Traits<RMatrix<Type_> >::sizeCols_,
      storage_   = hidden::Traits<RMatrix<Type_> >::storage_,

      Rtype_ = hidden::RcppTraits<Type_>::Rtype_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    typedef typename hidden::Traits<RMatrix<Type_> >::Row Row;
    typedef typename hidden::Traits<RMatrix<Type_> >::Col Col;

    /** Default Constructor */
    inline RMatrix(): matrix_(),rows_(), cols_() {}
    /** Constructor with given dimension */
    inline RMatrix(int nRow, int nCol): matrix_(nRow, nCol),rows_(nRow), cols_(nCol) {}
    /** Constructor with SEXP */
    inline RMatrix( SEXP robj)
                  : matrix_(robj),rows_(0, matrix_.rows()), cols_(0, matrix_.cols())
    {}
    /** Constructor */
    inline RMatrix( Rcpp::Matrix<Rtype_> matrix)
                  : matrix_(matrix), rows_(0, matrix_.rows()), cols_(0, matrix_.cols())
    {}
    /** Copy Constructor. @c ref is only there for compatibility */
    inline RMatrix( RMatrix const& matrix, bool ref = true)
                  : matrix_(matrix), rows_(0, matrix_.rows()), cols_(0, matrix_.cols())
    {}

    /** overwrite the RMatrix with mat, using Rcpp::operator=.
     *  @param mat the matrix to copy
     **/
    template<class OtherType>
    inline RMatrix& operator=( RMatrix<OtherType> const& mat)
    {
      matrix_ = mat.matrix_;  rows_ = mat.rows_; cols_ = mat.cols_;
      return *this;
    }
    /** overwrite the RMatrix with mat, using Rcpp::operator=.
     *  @param mat the matrix to copy
     **/
    template<int OtherRtype>
    inline RMatrix& operator=( Rcpp::Matrix<OtherRtype> const& mat)
    {
      matrix_ = mat; rows_ = mat.rows_; cols_ = mat.cols_;
      return *this;
    }

    /** @return the underlying Rcpp matrix */
    inline Rcpp::Matrix<Rtype_> const& matrix() const { return matrix_;}
    /** cast operator */
    inline operator Rcpp::Matrix<Rtype_>() const { return matrix_;}

    /** @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /**@return the range of the columns */
    inline ColRange const& colsImpl() const { return cols_;}

    /** set Matrix .
     *  @param matrix the Rcpp matrix to wrap
     *  @note cannot be passed as const& due to a bug from the Rcpop side
     * */
    inline void setMatrix( Rcpp::Matrix<Rtype_> matrix)
    { matrix_ = matrix;
      rows_ = RowRange(0, matrix_.rows());
      cols_ = RowRange(0, matrix_.cols());
    }

    /** @return a reference on column j
      *  @param j index of the column to wrap
      **/
    inline Col col( int j) const { return Col(*this, j);}
    /** @return a reference on column j
      *  @param j index of the column to wrap
      **/
    inline Row row( int i) const { return Row(*this, i);}
   /** @return a constant reference on element (i,j)
     *  @param i, j indexes of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    { return static_cast<ConstReturnType>(matrix_(i,j));}
    /** @return a reference on the element (i,j)
     *  @param i, j indexes of the row and of the column
     **/
    inline Type& elt2Impl(int i, int j) { return (matrix_(i,j));}
    /** overwrite the RMatrix with mat using Rcpp::operator=.
     *  @param mat the RMatrix to copy
     **/
    inline RMatrix& operator=( RMatrix const& mat)
    {
      matrix_ = mat.matrix_;  rows_ = mat.rows_; cols_ = mat.cols_;
      return *this;
    }
  private:
    Rcpp::Matrix<Rtype_> matrix_;
    RowRange rows_;
    ColRange cols_;
};

template <typename Type_>
class RowRMatrix: public ArrayBase< RowRMatrix<Type_> >, public TRef<1>
{
  public:
    typedef typename hidden::Traits<RowRMatrix<Type_> >::Type Type;
    typedef typename hidden::Traits<RowRMatrix<Type_> >::ReturnType ReturnType;
    enum
    {
      structure_ = hidden::Traits<RowRMatrix<Type_> >::structure_,
      orient_    = hidden::Traits<RowRMatrix<Type_> >::orient_,
      sizeRows_  = hidden::Traits<RowRMatrix<Type_> >::sizeRows_,
      sizeCols_  = hidden::Traits<RowRMatrix<Type_> >::sizeCols_,
      storage_   = hidden::Traits<RowRMatrix<Type_> >::storage_,

      Rtype_ = hidden::RcppTraits<Type_>::Rtype_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Reference Constructor. Create a slice of a RMatrix */
    inline RowRMatrix( RMatrix<Type> const& mat, int i)
                     : matrix_(mat.matrix()), rows_(i, 1), cols_(mat.cols())
    {}
    /** Reference Constructor. @c ref is only there for compatibility */
    inline RowRMatrix( RowRMatrix<Type> const& rowmat, bool ref =true)
                     : matrix_(rowmat.matrix_), rows_(rowmat.rows()), cols_(rowmat.cols())
    {}
  private:
    Rcpp::Matrix<Rtype_> matrix_;
    RowRange rows_;
    ColRange cols_;
};

template <typename Type_>
class ColRMatrix: public ArrayBase< ColRMatrix<Type_> >, public TRef<1>
{
  public:
    typedef typename hidden::Traits<ColRMatrix<Type_> >::Type Type;
    typedef typename hidden::Traits<ColRMatrix<Type_> >::ReturnType ReturnType;
    enum
    {
      structure_ = hidden::Traits<ColRMatrix<Type_> >::structure_,
      orient_    = hidden::Traits<ColRMatrix<Type_> >::orient_,
      sizeRows_  = hidden::Traits<ColRMatrix<Type_> >::sizeRows_,
      sizeCols_  = hidden::Traits<ColRMatrix<Type_> >::sizeCols_,
      storage_   = hidden::Traits<ColRMatrix<Type_> >::storage_,

      Rtype_ = hidden::RcppTraits<Type_>::Rtype_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Reference Constructor. Create a slice of a RMatrix */
    inline ColRMatrix( RMatrix<Type> const& mat, int j)
                     : matrix_(mat.matrix()), rows_(mat.rows()), cols_(j,1)
    {}
    /** Reference Constructor. @c ref is only there for compatibility */
    inline ColRMatrix( ColRMatrix<Type> const& colmat, bool ref =true)
                     : matrix_(colmat.matrix_), rows_(colmat.rows()), cols_(colmat.cols())
    {}
  private:
    Rcpp::Matrix<Rtype_> matrix_;
    RowRange rows_;
    ColRange cols_;
};

} // namespace STK


#endif /* STK_RMATRIX_H */
