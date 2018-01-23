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
 * Project:  stkpp::Sdk
 * created on: 9 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_StaticAssert.h
 *  @brief In this file we define a StaticAssert convenience structure
 *  for error checking at compile time.
 **/


#ifndef STK_STATICASSERT_H
#define STK_STATICASSERT_H

#include "STK_MetaTemplate.h"

#define STK_SINGLE_ARG2(A,B) A,B
#define STK_SINGLE_ARG3(A,B,C) A,B,C
#define STK_SINGLE_ARG4(A,B,C,D) A,B,C,D
#define STK_SINGLE_ARG5(A,B,C,D,E) A,B,C,D,E

#if(__cplusplus < 201103L) // use old way

namespace STK
{

template<int condition> struct StaticAssert;
// if the condition is false the compiler will complain it does not find the message
template<> struct StaticAssert<0> {};
template<> struct StaticAssert<1>
{
  enum
  {
    YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX,
    YOU_TRIED_CALLING_A_MATRIX_METHOD_ON_A_VECTOR,
    YOU_TRIED_CALLING_A_POINT_METHOD_ON_SOMETHING_ELSE,
    YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_SOMETHING_ELSE,
    YOU_TRIED_TO_CONSTRUCT_A_REFERENCE_WITH_A_DIFFERENT_ORIENTATION,
    YOU_TRIED_TO_CONSTRUCT_A_SQUAREMATRIX_WITH_DIFFERENT_DIMENSIONS,
    YOU_TRIED_TO_CONSTRUCT_A_VECTOR_WITH_MORE_THAN_ONE_COLUM,
    YOU_TRIED_TO_CONSTRUCT_A_POINT_WITH_MORE_THAN_ONE_ROW,
    YOU_TRIED_TO_CONSTRUCT_A_SCALAR_WITH_MORE_THAN_ONE_ROW_OR_ONE_COLUMN,
    YOU_TRIED_TO_APPLY_A_BINARY_OPERATOR_WITH_WRONG_ROWS_SIZE,
    YOU_TRIED_TO_APPLY_A_BINARY_OPERATOR_WITH_WRONG_COLS_SIZE,
    YOU_TRIED_TO_APPLY_A_BINARY_OPERATOR_BETWEEN_NOT_COMPATIBLE_ARRAYS,
    YOU_TRIED_TO_ASSIGN_A_NOT_COMPATIBLE_ARRAY,
    YOU_TRIED_TO_APPLY_A_METHOD_BETWEEN_NOT_COMPATIBLE_ARRAYS,
    YOU_TRIED_TO_APPLY_A_PRODUCT_OPERATOR_BETWEEN_NOT_COMPATIBLE_ARRAYS,
    YOU_TRIED_TO_APPLY_A_BINARY_OPERATOR_WITH_WRONG_DIMENSIONS_SIZE,
    YOU_TRIED_TO_WRAP_A_CONTAINER_WITH_THE_WRONG_TYPE_AS_A_VARIABLE,
    YOU_TRIED_TO_WRAP_A_CONTAINER_WHICH_IS_NOT_1D_AS_A_VARIABLE,
    YOU_TRIED_TO_USE_A_UNIDIMENSIONAL_METHOD_ON_A_MATRIX,
    YOU_CANNOT_USED_AN_UNIDIMENSIONAL_METHOD_ON_THIS_OBJECT,
    YOU_CANNOT_USED_A_ZERODIMENSIONAL_METHOD_ON_THIS_OBJECT,
    YOU_CANNOT_USED_THIS_TYPE_OF_DATA_WITH_THIS_OBJECT,
    YOU_CANNOT_USED_THIS_METHOD_WITH_THIS_KIND_OF_ARRAY,
    YOU_HAVE_TO_USE_TWO_ARRAYS_WITH_THE_SAME_ORIENTATION,
    YOU_HAVE_TO_USE_A_VECTOR_OR_POINT_IN_THIS_METHOD,
    YOU_HAVE_TO_USE_A_SQUARE_MATRIX_IN_THIS_METHOD,
    YOU_HAVE_TO_USE_A_MATRIX_OR_SQUARE_MATRIX_IN_THIS_METHOD,
    INVALID_VECTOR_VECTOR_PRODUCT,
    INVALID_POINT_POINT_PRODUCT
  };
};

} // namespace STK

#ifndef IS_RTKPP_LIB
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#endif
#define STK_STATIC_ASSERT(COND,MSG) \
if (STK::StaticAssert<bool(COND)>::MSG) {}

#else /* C++11  */

// if native static_assert is enabled, let's use it
#include <type_traits>
#define STK_STATIC_ASSERT(COND,MSG) static_assert(bool(COND),#MSG);

#endif // C++11

#define STK_STATIC_ASSERT_SAME_ORIENTATION(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_CONSTRUCT_A_REFERENCE_WITH_A_DIFFERENT_ORIENTATION)

#define STK_STATIC_ASSERT_DIMENSIONS_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_APPLY_A_METHOD_BETWEEN_NOT_COMPATIBLE_ARRAYS)

#define STK_STATIC_ASSERT_SQUARE_SIZES_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_CONSTRUCT_A_SQUAREMATRIX_WITH_DIFFERENT_DIMENSIONS)

#define STK_STATIC_ASSERT_VECTOR_SIZECOLS_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_CONSTRUCT_A_VECTOR_WITH_MORE_THAN_ONE_COLUM)

#define STK_STATIC_ASSERT_POINT_SIZEROWS_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_CONSTRUCT_A_POINT_WITH_MORE_THAN_ONE_ROW)

#define STK_STATIC_ASSERT_SCALAR_SIZE_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_CONSTRUCT_A_SCALAR_WITH_MORE_THAN_ONE_ROW_OR_ONE_COLUMN)

#define STK_STATIC_ASSERT_ROWS_DIMENSIONS_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_APPLY_A_BINARY_OPERATOR_WITH_WRONG_ROWS_SIZE)

#define STK_STATIC_ASSERT_COLS_DIMENSIONS_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_APPLY_A_BINARY_OPERATOR_WITH_WRONG_COLS_SIZE)

#define STK_STATIC_ASSERT_BINARY_OPERATOR_DIMENSIONS_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_APPLY_A_BINARY_OPERATOR_WITH_WRONG_DIMENSIONS_SIZE)

#define STK_STATIC_ASSERT_BINARY_OPERATOR_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_APPLY_A_BINARY_OPERATOR_BETWEEN_NOT_COMPATIBLE_ARRAYS)

#define STK_STATIC_ASSERT_PRODUCT_OPERATOR_MISMATCH(COND) \
    STK_STATIC_ASSERT(COND,YOU_TRIED_TO_APPLY_A_PRODUCT_OPERATOR_BETWEEN_NOT_COMPATIBLE_ARRAYS)

#define STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(EXPR) \
STK_STATIC_ASSERT( (int(hidden::Traits<EXPR>::structure_)==(int)Arrays::array2D_) \
                 ||(int(hidden::Traits<EXPR>::structure_)==(int)Arrays::square_) \
                 ||(int(hidden::Traits<EXPR>::structure_)==(int)Arrays::diagonal_) \
                 ||(int(hidden::Traits<EXPR>::structure_)==(int)Arrays::number_) \
                 ||(int(hidden::Traits<EXPR>::structure_)==(int)Arrays::lower_triangular_) \
                 ||(int(hidden::Traits<EXPR>::structure_)==(int)Arrays::upper_triangular_) \
                 ||(int(hidden::Traits<EXPR>::structure_)==(int)Arrays::symmetric_) \
                 ||(int(hidden::Traits<EXPR>::structure_)==(int)Arrays::lower_symmetric_) \
                 ||(int(hidden::Traits<EXPR>::structure_)==(int)Arrays::upper_symmetric_) \
                 ,YOU_TRIED_CALLING_A_MATRIX_METHOD_ON_A_VECTOR)

#define STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(EXPR) \
STK_STATIC_ASSERT(  STK_SINGLE_ARG2(hidden::IsEqual< int(hidden::Traits<EXPR>::structure_), int(Arrays::vector_)>::value_) \
                  ||STK_SINGLE_ARG2(hidden::IsEqual< int(hidden::Traits<EXPR>::structure_), int(Arrays::point_)>::value_) \
                  ||STK_SINGLE_ARG2(hidden::IsEqual< int(hidden::Traits<EXPR>::structure_), int(Arrays::diagonal_)>::value_) \
                  ||STK_SINGLE_ARG2(hidden::IsEqual< int(hidden::Traits<EXPR>::structure_), int(Arrays::number_)>::value_) \
                  ,YOU_CANNOT_USED_AN_UNIDIMENSIONAL_METHOD_ON_THIS_OBJECT)

#define STK_STATIC_ASSERT_SAME_ORIENT(LHS, RHS) \
STK_STATIC_ASSERT( hidden::Traits<LHS>::structure_== hidden::Traits<RHS>::structure_ \
                 ,YOU_HAVE_TO_USE_TWO_ARRAYS_WITH_THE_SAME_ORIENTATION)

#define STK_STATIC_ASSERT_POINT_ONLY(EXPR) \
STK_STATIC_ASSERT(int(hidden::Traits<EXPR>::structure_) == int(Arrays::point_) \
                 ,YOU_TRIED_CALLING_A_POINT_METHOD_ON_SOMETHING_ELSE)

#define STK_STATIC_ASSERT_VECTOR_ONLY(EXPR) \
STK_STATIC_ASSERT(STK_SINGLE_ARG2(hidden::IsEqual< int(hidden::Traits<EXPR>::structure_), int(Arrays::vector_)>::value_) \
                 ,YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_SOMETHING_ELSE)

#define STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(EXPR) \
STK_STATIC_ASSERT(  STK_SINGLE_ARG2(hidden::IsEqual< int(hidden::Traits<EXPR>::structure_), int(Arrays::vector_)>::value_) \
                  ||STK_SINGLE_ARG2(hidden::IsEqual< int(hidden::Traits<EXPR>::structure_), int(Arrays::point_)>::value_) \
                  ,YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_SOMETHING_ELSE)

#define STK_STATIC_ASSERT_ZERO_DIMENSION_ONLY(EXPR) \
STK_STATIC_ASSERT(STK_SINGLE_ARG2(hidden::IsEqual< int(hidden::Traits<EXPR>::structure_), int(Arrays::number_)>::value_) \
                 ,YOU_CANNOT_USED_A_ZERODIMENSIONAL_METHOD_ON_THIS_OBJECT)

#endif /* STK_STATICASSERT_H */
