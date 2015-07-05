/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

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
 * created on: 25 déc. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ProductRaw.h
 *  @brief In this file we implement the raw static methods used by the products methods.
 **/

#ifndef STK_PRODUCTRAW_H
#define STK_PRODUCTRAW_H

namespace STK
{
/* size of the block and panels used in the product algorithm */
const int blockSize = 4;
const int panelSize = 64;
const int vectorSize = 256;

namespace hidden
{
/** @ingroup hidden
 *  This structure encapsulate the data allocated for a panel.
 **/
template<class Type>
struct Panel
{
  Type panel[blockSize*panelSize];
  inline Type const& operator[](int i) const { return panel[i];}
  inline Type& operator[](int i) { return panel[i];}
};

/** @ingroup hidden
 *  This structure encapsulate the data allocated for a block.
 **/
template<class Type>
struct Block
{
  Type block[blockSize*blockSize];
  inline Type const& operator[](int i) const { return block[i];}
  inline Type& operator[](int i) { return block[i];}
};

/** @ingroup hidden
 *  This structure encapsulate the data allocated for a block.
 **/
template<class Type>
struct RawVec
{
  Type vec[panelSize];
  inline Type const& operator[](int i) const { return vec[i];}
  inline Type& operator[](int i) { return vec[i];}
};

/** @ingroup hidden
 *  This structure regroup the methods to used after block multiplication in
 *  order to perform the product of the remaining rows and columns.
 **/
template<typename Lhs, typename Rhs, typename Result>
struct MultCoefImpl
{
  typedef typename Result::Type Type;
  enum
  {
    sizeRows_  = Result::sizeRows_,
    sizeCols_  = Result::sizeCols_,
    orient_    = Result::orient_,
    storage_   = Result::storage_
  };
  /** dot product. general by general*/
  static void dot( Lhs const& lhs, Rhs const& rhs, Result& res, int iRow, int jCol)
  {
    res.elt(iRow, jCol) = Type(0);
    Range const dotRange = Range::inf(lhs.rangeColsInRow(iRow), rhs.rangeRowsInCol(jCol));
    for (int k=dotRange.begin(); k< dotRange.end(); ++k)
      res.elt(iRow, jCol) += lhs.elt(iRow, k) * rhs.elt(k, jCol);
  }
  /** dot product. general by vector */
  static void dot( Lhs const& lhs, ITContainer<Rhs, Arrays::vector_> const& rhs
                 , ICAllocator<Result>& res, int iRow)
  {
    res.elt(iRow) = Type(0);
    Range const dotRange = Range::inf(lhs.rangeColsInRow(iRow), rhs.range());
    for (int k=dotRange.begin(); k< dotRange.end(); ++k)
      res.elt(iRow) += lhs.elt(iRow, k) * rhs.elt(k);
  }
  /** dot product. general by vector */
  static void dot( ITContainer<Lhs, Arrays::point_> const& lhs
                 , Rhs const& rhs, ICAllocator<Result>& res, int jCol)
  {
    res.elt(jCol) = Type(0);
    Range const dotRange = Range::inf(rhs.rangeRowsInCol(jCol), lhs.range());
    for (int k=dotRange.begin(); k< dotRange.end(); ++k)
      res.elt(jCol) += lhs.elt(k) * rhs.elt(k, jCol);
  }
  /** multiplication with one sized vectors */
  static void mult1( Lhs const& lhs, Rhs const& rhs, Result& res
                   , int lhsCol, int rhsRow)
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
        res.elt(i,j) += lhs.elt(i, lhsCol) * rhs.elt(rhsRow, j);
  }
  /** multiplication with two sized vectors */
  static void mult2( Lhs const& lhs, Rhs const& rhs, Result& res
                   , int lhsCol, int rhsRow)
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
        res.elt(i,j) += lhs.elt(i, lhsCol  ) * rhs.elt(rhsRow,   j)
                      + lhs.elt(i, lhsCol+1) * rhs.elt(rhsRow+1, j);
  }
  /** multiplication with three sized vectors */
  static void mult3( Lhs const& lhs, Rhs const& rhs, Result& res
                   , int lhsCol, int rhsRow)
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
        res.elt(i,j) += lhs.elt(i, lhsCol  ) * rhs.elt(rhsRow, j)
                      + lhs.elt(i, lhsCol+1) * rhs.elt(rhsRow+1, j)
                      + lhs.elt(i, lhsCol+2) * rhs.elt(rhsRow+2, j);
  }
  /** multiplication with one sized vectors */
  static void multVec1( Lhs const& lhs, Rhs const& rhs, Result& res
                      , int lhsCol, int rhsRow)
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      res.elt(i) += lhs.elt(i, lhsCol) * rhs.elt(rhsRow);
  }
  /** multiplication with two sized vectors */
  static void multVec2( Lhs const& lhs, Rhs const& rhs, Result& res
                      , int lhsCol, int rhsRow)
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      res.elt(i) += lhs.elt(i, lhsCol) * rhs.elt(rhsRow)
                  + lhs.elt(i, lhsCol+1) * rhs.elt(rhsRow+1);
  }
  /** multiplication with three sized vectors */
  static void multVec3( Lhs const& lhs, Rhs const& rhs, Result& res
                      , int lhsCol, int rhsRow)
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); ++i)
      res.elt(i) += lhs.elt(i, lhsCol) * rhs.elt(rhsRow)
                  + lhs.elt(i, lhsCol+1) * rhs.elt(rhsRow+1)
                  + lhs.elt(i, lhsCol+2) * rhs.elt(rhsRow+2);
  }
  /** multiplication with one sized vectors */
  static void multPoint1( Lhs const& lhs, Rhs const& rhs, Result& res
                        , int lhsCol, int rhsRow)
  {
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      res.elt(j) += lhs.elt(lhsCol) * rhs.elt(rhsRow, j);
  }
  /** multiplication with two sized vectors */
  static void multPoint2( Lhs const& lhs, Rhs const& rhs, Result& res
                        , int lhsCol, int rhsRow)
  {
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      res.elt(j) += lhs.elt(lhsCol) * rhs.elt(rhsRow, j)
                  + lhs.elt(lhsCol+1) * rhs.elt(rhsRow+1, j);
  }
  /** multiplication with three sized vectors */
  static void multPoint3( Lhs const& lhs, Rhs const& rhs, Result& res
                        , int lhsCol, int rhsRow)
  {
    for (int j=rhs.beginCols(); j< rhs.endCols(); ++j)
      res.elt(j) += lhs.elt(lhsCol) * rhs.elt(rhsRow, j)
                  + lhs.elt(lhsCol+1) * rhs.elt(rhsRow+1, j)
                  + lhs.elt(lhsCol+2) * rhs.elt(rhsRow+2, j);
  }
};


} // namespace hidden

} // namespace STK

namespace STK
{

namespace hidden
{

} // namespace hidden

} // namespace STK

#endif /* STK_PRODUCTIMPL_H */
