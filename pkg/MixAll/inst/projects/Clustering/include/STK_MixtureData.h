/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013 Serge Iovleff

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
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MixtureData.h
 *  @brief In this file we define the data manager class associated to the MixtureBridge class.
 **/

#ifndef STK_MIXTUREDATA_H
#define STK_MIXTUREDATA_H

#include <DManager/include/STK_IDataBridge.h>

namespace STK
{

/** @ingroup Clustering
 *  @brief bridge a data set for a MixtureBridge.
 *
 * @tparam Data The data bridged by the MixtureData class
 */
template<class Data>
class MixtureData: public IDataBridge
{
  public:
    typedef IDataBridge Base;
    typedef typename Data::Type Type;
    /** default constructor. */
    inline MixtureData(std::string const& idData): Base(idData), dataij_(), nbVariable_(0) {}
    /** copy constructor (Warning: will copy the data set)
     *  @param manager the MixtureData to copy
     **/
    MixtureData( MixtureData const& manager)
               : Base(manager)
               , dataij_(manager.dataij_)
               , nbVariable_(manager.nbVariable_) {}
    /** destructor */
    virtual ~MixtureData() {}
    /** getter. @return the number of variables (the number of columns of the data)  */
    inline int nbVariable() const { return nbVariable_;}
    /** getter. @return a constant reference on the data set */
    inline Data const& dataij() const { return dataij_;}

    /** setter. @return a reference on the number of variables */
    inline int& nbVariable() { return nbVariable_;}
    /** setter. @return a reference on the data set */
    inline Data& dataij() { return dataij_;}

   /** get the (imputed) missing values of a data set.
    *  @note In C++11, it will be possible to use a tuple rather that this
    *  pair of pair...
    *  @param data the array to return with the missing values
    **/
   template<typename Type_>
   void getMissingValues( std::vector< std::pair<std::pair<int,int>, Type_ > >& data) const;

 protected:
   /** data set */
   Data dataij_;
   /** number of variables in the data set */
   int nbVariable_;
   /** utility function for lookup the data set and find missing values
    *  coordinates.
    *  @return the number of missing values
    **/
  virtual std::vector< std::pair<int,int> >::size_type findMissing();
};

template<class Data>
std::vector< std::pair<int,int> >::size_type MixtureData<Data>::findMissing()
{
#ifdef STK_DMANAGER_VERBOSE
  stk_cout << _T("IDataBridge::Entering findMissing()\n");
#endif
  for (int j=dataij_.beginCols(); j< dataij_.endCols(); ++j)
  {
    for (int i=dataij_.beginRows(); i< dataij_.endRows(); ++i)
    {
      if (Arithmetic<Type>::isNA(dataij_(i,j)))
      { v_missing_.push_back(std::pair<int,int>(i,j));}
    }
  }
 return v_missing_.size();
#ifdef STK_DMANAGER_VERBOSE
  stk_cout << _T("IDataBridge::findMissing() terminated, nbMiss= ") << v_missing_.size() << _T("\n");
#endif
}

template<class Data>
template<typename Type_>
void MixtureData<Data>::getMissingValues( std::vector< std::pair<std::pair<int,int>, Type_ > >& data) const
{
  data.resize(v_missing_.size());
  for(size_t i = 0; i< v_missing_.size(); ++i)
  {
    data[i].first  = v_missing_[i];
    data[i].second = dataij_(v_missing_[i].first, v_missing_[i].second);
  }
}
} // namespace STK

#endif /* STK_MIXTUREDATA_H */
