/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

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
 * Project:  stkpp:stkernel::base
 * Purpose:  Define the fundamental type Integer.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Integer.h
 *  @brief In this file we define the fundamental type Integer.
 **/

#ifndef STK_INTEGER_H
#define STK_INTEGER_H

// for building
#ifdef RTKPP_LIB
#include <Rcpp.h>
#endif

#include "STK_String.h"
#include <map>

namespace STK
{
/** @ingroup Base
  *  @brief STK fundamental type of discrete values.
  *
  *  The type Interger is defined for the numerical computation and the
  *  internal representation of the discrete variables.
  **/
typedef int Integer ;

template<> struct Arithmetic<Integer>;
template<> struct IdTypeImpl<Integer>;

/** @ingroup Arithmetic
 *  @brief Specialization for Integer (long).
 * 
 *  We are using the largest element of the underlying
 *  Type for representing NA (not available) discrete numbers.
 */
template<>
struct Arithmetic<Integer>  : public std::numeric_limits<Integer>
{
#ifdef RTKPP_LIB
  enum
  {
    Rtype_ = hidden::RcppTraits<Integer>::Rtype_
  };
#endif
  /** True if the type has a representation for a "Not Available". */
  static const bool hasNA = true;
  /** Adding a Non Avalaible (NA) special number. */
#ifdef RTKP_LIB
  static inline Integer NA() throw()  { return Rcpp::traits::get_na<Rtype_>();}
#else
  static inline Integer NA() throw()
  { return std::numeric_limits<Integer>::max();}
  /** We are using the maximal value of the Integer Type for NA values. */
  static inline Integer max() throw()
  { return std::numeric_limits<Integer>::max() -1; }
#endif
  /** Test if x is a Non Avalaible (NA) special number.
   *  @param x the Integer number to test.
   **/
#ifdef RTKP_LIB
  static bool isNA(Integer const& x) throw() { Rcpp::traits::is_na<Rtype_>(x);}
#else
  static inline bool isNA(Integer const& x) throw()
  { return (x==std::numeric_limits<Integer>::max());}
#endif
  /** @return @c true if x is  infinite : always false for Integer.
   *  @param x the Integer number to test.
   **/
  static inline bool isInfinite(Integer const& x) throw() { return false; }
  /** @return @¢ true if x is  finite : i.e. if x is not a NA value.
   *  @param x the value to test.
   **/
  static inline bool isFinite(Integer const& x) throw() { return (!isNA(x));}
};

/** @ingroup RTTI 
 *  @brief Specialization of the IdTypeImpl for the type Integer.
 **/
template<>
struct IdTypeImpl<Integer>
{
  /** Give the IdType of the type Integer. */
  static inline IdType returnType() { return(integer_);}
};

/** @ingroup Base
 *  @brief Convert a String to an Integer.
 *  @param str the String we want to convert
 *  @return the Integer represented by the String @c str. If the string
 *  does not match any known name, the NA value is returned.
 **/
Integer stringToInt( String const& str);

/** @ingroup Base
 *  Convert a String to an Integer using a map.
 *  @param str the String we want to convert
 *  @param mapping the mapping between the string and the Int
 *  @return the Int represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_ type is returned.
 **/
Integer stringToInt( String const& str, std::map<String, Integer> const& mapping);

/** @ingroup Base
 *  Convert an Integer to a String.
 *  @param value the Integer we want to convert
 *  @param f format, by default write every number in decimal
 *  @return the string associated to this value.
 **/
String intToString( Integer const& value, std::ios_base& (*f)(std::ios_base&) = std::dec);

/** @ingroup Base
 *  Convert an Integer to a String.
 *  @param value the Integer we want to convert
 *  @param mapping the mapping between Integer and String
 *  @return the String associated to this value.
 **/
String intToString( Integer const& value, std::map<Integer, String> const& mapping);

/** @ingroup Base
 *  @brief specialization for Integer
 *  @param s the String to convert
 *  @return The value to get from the String
 **/
template<>
inline Integer stringToType<Integer>( String const& s)
{ return stringToInt(s);}

/** @ingroup Base
 *  @brief specialization for Integer
 *  @param t The Int to convert to String
 *  @param f format, by default write every number in decimal
 **/
template<>
inline String typeToString<Integer>( Integer const& t, std::ios_base& (*f)(std::ios_base&))
{ return intToString(t, f);}

} // namespace STK

#endif /*STK_INTEGER_H*/
