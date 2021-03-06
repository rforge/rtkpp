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
 * Project:  stkpp::DManager
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Variable.h
 *  @brief Implement the interface class IVariable as Variable.
 **/

#ifndef STK_VARIABLE_H
#define STK_VARIABLE_H

#include <Arrays/include/STK_IArray1D.h>
#include "STK_IVariable.h"

namespace STK
{

template< class Type_> class Variable;

namespace hidden
{
/** @ingroup hidden
 *  Specialization of the Traits struct for Variable class.
 **/
template<class Type_>
struct Traits< Variable<Type_> >
{
  typedef Variable<Type_> Row;
  typedef Variable<Type_> Col;
  typedef Variable<Type_> SubRow;
  typedef Variable<Type_> SubCol;
  typedef Variable<Type_> SubArray;
  typedef Variable<Type_> SubVector;

  typedef Type_          Type;
  typedef typename RemoveConst<Type_>::Type const& TypeConst;

  enum
  {
    structure_ = Arrays::vector_,
    orient_    = Arrays::by_col_,
    sizeCols_  = 1,
    sizeRows_  = UnknownSize,
    size_      = UnknownSize,
    storage_   = Arrays::dense_ // always dense
  };

  typedef MemAllocator<Type_, UnknownSize> Allocator;

  typedef TRange<size_> RowRange;
  typedef TRange<1>     ColRange;

  typedef int Index;

  typedef DenseRandomIterator<Variable<Type_> > Iterator;
  typedef ConstDenseRandomIterator<Variable<Type_> > ConstIterator;

  typedef std::reverse_iterator<Iterator> ReverseIterator;
  typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;
};

} // namespace hidden
/** @ingroup DManager
  * @brief Variable is an implementation of the Base class IVariable
  * using The Array1D class for storing the data.
  * It implements all purely virtual methods defined in the IVariable base
  * class.
 **/
template<class Type_>
class Variable: public IVariable
              , public IArray1D< Variable<Type_> >
{
  public:
    typedef Variable<Type_> Row;
    typedef Variable<Type_> Col;
    typedef Variable<Type_> SubRow;
    typedef Variable<Type_> SubCol;
    typedef Variable<Type_> SubArray;
    typedef Variable<Type_> SubVector;

    typedef typename hidden::Traits<Variable<Type_> >::Type Type;
    typedef typename hidden::Traits<Variable<Type_> >::TypeConst TypeConst;

    enum
    {
      structure_ = hidden::Traits< Variable<Type_> >::structure_,
      orient_    = hidden::Traits< Variable<Type_> >::orient_,
      size_      = hidden::Traits< Variable<Type_> >::size_,
      sizeCols_  = hidden::Traits< Variable<Type_> >::sizeCols_,
      sizeRows_  = hidden::Traits< Variable<Type_> >::sizeRows_,
      storage_   = hidden::Traits< Variable<Type_> >::storage_
    };

    typedef typename hidden::Traits< Variable<Type_> >::Allocator Allocator;

    //typedef MemAllocator<Type_, size_> Allocator;
    typedef IArray1D< Variable<Type_> > Base;

    using Base::elt;
    /** default constructor
     *  @param name name of the variable
     **/
    explicit Variable( String const& name = stringNa)
                     : IVariable(IdTypeImpl<Type_>::returnType(), name)
                     , Base()
    {}
    /** constructor
     *  @param size,name size and name of the variable
     **/
    explicit Variable( int size, String const& name = stringNa)
                     : IVariable(IdTypeImpl<Type_>::returnType(), name)
                     , Base(size)
    {}
    /** Default constructor
     *  @param I,name range and name of the variable
     **/
    explicit Variable( Range const& I, String const& name = stringNa)
                     : IVariable(IdTypeImpl<Type_>::returnType(), name)
                     , Base(I)
    {}
    /** constructor with specified initial value
     *  @param I,v range of the data and initial value
     *  @param name name of the variable
     **/
    Variable( Range const& I, Type_ const& v, String const& name)
            : IVariable(IdTypeImpl<Type_>::returnType(), name)
            , Base(I)
    { this->setValue(v);}
    /** copy/reference constructor.
     *  @param V the Variable to copy
     *  @param ref @c true if we want to wrap V
     **/
    explicit Variable( Variable const& V, bool ref = false)
                     : IVariable(V)
                     , Base(V, ref)
    {}
    /** reference constructor
     *  @param V,I Variable and range to wrap
     *  @param ref @c true if we want to wrap V
     **/
    explicit Variable( Variable const& V, Range const& I, bool ref = true)
                     : IVariable(V), Base(V, I, true)
    {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     **/
    template<class OtherArray>
    Variable( IArray1D<OtherArray> const& T)
            : IVariable(IdTypeImpl<Type_>::returnType(), stringNa)
            , Base(T, T.range(), true)
    {}
    /** destructor. */
    ~Variable() {}
    /** clone return a ptr on a copy of the Object.
     *  @param ref true if we want just a reference
     **/
    virtual Variable* clone( bool ref = false) const
    { return new Variable(*this, ref);}
    /** New first index for the object.
     *  @param rbeg the index of the first row to set
     **/
    void shift1D(int rbeg) { Base::shift(rbeg, this->beginCols());}
    /**  Resize the container.
     *  @param I the range to set to the container
     **/
    inline Variable<Type_>& resize1D(Range const& I)
    { Base::resize(I, this->cols()); return *this;}

    /** remove n elements to the end of the container
     *  @param n number of element to remove */
    virtual void popBack(int n =1) { Base::popBack(n);}
    /** Add n elements to the container.
     *  @param n number of elements to add
     **/
    void pushBack( int n=1) { Base::pushBack(n);}
    /** Add an element to the container.
     *  @param v the element to add
     **/
    void push_back( Type_ const& v)
    { Base::pushBack();
      this->back() = v;
    }
    /** Insert n elements at the position pos of the container. The bound
     *  end_ should be modified at the very end of the insertion as pos
     *  can be a reference to it.
     *  @param pos,n index and number of elements to insert
     **/
    void insertElt(int pos, int n =1)
    { Base::insertElt(pos, n);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param V the container to copy
     **/
    inline Variable& operator=(Variable const& V)
    { // copy IVariable part
      this->name_ = V.name_;
      return Base::assign(V);
    }
    /** set the container to a constant value.
     *  @param v the value to set
     **/
    inline Variable& operator=(Type_ const& v) { this->setValue(v); return *this;}

    /** move the variable in this
     *  @param V variable to move in this
     **/
    inline void move(Variable const& V)
    { Base::move(V);
      name_ = V.name_;
    }
    /** encode values as ints. Not used yet. */
    void encode()
    { int code = baseIdx;
      std::pair< typename std::map<Type_, int>::iterator, bool> ret;
      for (int i=this->begin(); i< this->end(); i++)
      { ret=coding_.insert(std::pair<Type_, int>(this->elt(i), i));
        if (ret.second==true) { code++;}
      }
    }
    /** @return the maximal size of all the fields as String in the variable */
    int maxLength(bool with_name) const
    {
      typename String::size_type maxlength = with_name ? this->name().size() : 0;
      // loop over the values
      for (int i=this->begin(); i<this->end(); i++)
      { maxlength = std::max(maxlength, this->template eltAsString<Type_>(i).size() );}
      return int(maxlength);
    }
    /** @return the number of missing values in the variable */
    int nbMiss() const;

    // IVariable part
    /** @return the number of sample in the variable */
    virtual int size() const { return Base::size();}
    /** clear Container from all elements and memory allocated. */
    virtual void clear() { Base::clear();}
    /** Delete n elements at the @c pos index from the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
     **/
    virtual void erase(int pos, int n=1) { Base::erase(pos,n);}
    /** New first index for the object.
     *  @param beg the index of the first column to set
     **/
    virtual void shift(int beg) { Base::shift(beg);}
    /** push back n NA values.
     *  @param n number of NA values to add
     **/
    virtual inline void pushBackNAValues(int n=1);
    /** Overwrite the variable V by converting the data into strings.
     *  @param V Variable of String
     **/
    virtual inline void exportAsString( Variable< String >& V) const;
    /** Operator << : overwrite the Variable by converting the Strings
     *  contained in V into the Type_.
     *  @param V the Variable of string to import
     **/
    virtual inline Variable& operator<<( Variable< String > const& V);
    /** Operator >> : convert the Variable V into strings.
     *  @param V Variable of String
     **/
    virtual inline Variable<Type_> const& operator>>(Variable< String >& V) const;

    /** overwrite the Variable by converting the strings
     *  contained in V into the Type_.
     *  @param V Variable of String
     *  @param f io flags
     *  @return number of successful conversion
     **/
    inline int importFromString( Variable< String > const& V
                               , std::ios_base& (*f)(std::ios_base&) = std::dec
                               );
  protected:
    /** store the map String <-> int */
    std::map<Type_, int> coding_;
};

/** @return the number of missing values in the variable */
template<class Type_>
int Variable<Type_>::nbMiss() const
{
  int nbMiss = 0;
  // loop over the values
  for (int i=this->begin(); i<this->end(); i++)
  { if (Arithmetic<Type_>::isNA(this->elt(i))) nbMiss++;}
  return nbMiss;
}
/** push back n NA values. Specialization for Type_ = String.
 *  @param n number of NA values to add
 **/
template<>
inline void Variable<String>::pushBackNAValues(int n)
{ int first = this->end(), end = first+n;
  this->insertElt(this->end(), n);
  for (int i=first; i<end; i++)
  this->elt(i) = stringNa;
}
/** push back n NA values.
 *  @param n number of NA values to add
 **/
template<class Type_>
inline void Variable<Type_>::pushBackNAValues(int n)
{
    int first = this->end(), end = first+n;
    this->insertElt(this->end(), n);
    for (int i=first; i<end; i++)
    this->elt(i) =  Arithmetic<Type_>::NA();
}
/** overwrite the Variable by converting the strings
 *  contained in V into the Type_. Give the number of success.
 *  @param V Variable of String
 *  @param f io flags
 *  @return number of successful conversion
 **/
template<>
inline int Variable<String>::importFromString( Variable< String > const& V
                                             , std::ios_base& (*f)(std::ios_base&)
                                             )
{ *this = V; return V.size();}
/** overwrite the Variable by converting the strings
 *  contained in V into the Type_.
 *  @param V Variable of String
 *  @param f io flags
 *  @return number of successful conversion
 **/
template<class Type_>
inline int Variable<Type_>::importFromString( Variable< String > const& V
                                , std::ios_base& (*f)(std::ios_base&)
                                )
{
  this->resize(V.range());
  this->setName(V.name());
  int nSuccess = V.size();
  for (int i=V.begin(); i<V.end(); i++)
    if ( (Arithmetic<String>::isNA(V[i])) || (V[i]==stringNa) ) // not Available
      this->elt(i) = Arithmetic<Type_>::NA();
    else
    if (!stringToType<Type_>(this->elt(i), V[i], f)) nSuccess--;
  return nSuccess;
}
/** Overwrite the variable V by converting the data into strings.
 *  @param V Variable of String
 **/
template<>
inline void Variable<String>::exportAsString( Variable< String >& V) const
{ V = *this;}
/** Overwrite the variable V by converting the data into strings.
 *  @param V Variable of String
 **/
template<class Type_>
inline void Variable<Type_>::exportAsString( Variable< String >& V) const
{
  V.resize(this->range());
  V.setName(this->name());
  for (int i=this->begin(); i<=this->lastIdx(); i++)
  { V[i] = this->template eltAsString<Type_>(i);}
}
/** Operator << : overwrite the Variable by converting the strings
 *  contained in V into the String.
 *  @param V the Variable of string to import
 **/
template<>
inline Variable<String>& Variable<String>::operator<<( Variable< String > const& V)
{
  this->resize(V.range());
  this->setName(V.name());
  for (int i=V.begin(); i<V.end(); i++) this->elt(i) = V[i];
  return *this;
}
/** Operator << : overwrite the Variable by converting the Strings
 *  contained in V into the Type_.
 *  @param V the Variable of string to import
 **/
template<class Type_>
inline Variable<Type_>& Variable<Type_>::operator<<( Variable< String > const& V)
{
  this->resize(V.range());
  this->setName(V.name());
  for (int i=V.begin(); i<V.end(); i++) this->elt(i) = stringToType<Type_>(V[i]);
  return *this;
}
/** Operator >> : convert the Variable V into strings.
 *  @param V Variable of String
 **/
template<>
inline Variable<String> const& Variable<String>::operator>>(Variable< String >& V) const
{
  V.resize(this->range());
  V.setName(this->name());
  for (int i=this->begin(); i<this->end(); i++) V[i] = this->elt(i);
  return *this;
}
/** Operator >> : convert the Variable V into strings.
 *  @param V Variable of String
 **/
template<class Type_>
inline Variable<Type_> const& Variable<Type_>::operator>>(Variable< String >& V) const
{
  V.resize(this->range());
  V.setName(this->name());
  for (int i=this->begin(); i<this->end(); i++) V[i] = typeToString<Type_>(this->elt(i));
  return *this;
}
/** ostream for Variable. */
template<typename Type_>
inline ostream& operator<<(ostream& s, Variable<Type_> const& V)
{
  s << V.name() << _T("\n");
  return out1D(s, V);
}

} // namespace STK

#endif //STK_VARIABLE_H
