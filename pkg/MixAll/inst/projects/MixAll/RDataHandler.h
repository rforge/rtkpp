/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

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
 * Project:  MixAll
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file  RDataHandler.h
 *  @brief In this file we define the RDataHandler class.
 **/


#ifndef STK_RDATAHANDLER_H
#define STK_RDATAHANDLER_H

#include <RTKpp.h>
#include <Clustering.h>

// forward declaration
namespace STK
{ class RDataHandler;
}

namespace STK
{
namespace hidden
{
/** @ingroup hidden
 *  Specialization of the  DataHandlerTraits for DataHandler
 **/
template<typename Type>
struct DataHandlerTraits< STK::RDataHandler, Type>
{
  typedef CArray<Type> Data;
};

} // namespace hidden

} // namespace STK

namespace STK
{

/** @ingroup DManager
 *  The RDataHandler class allow to store various Rcpp::Matrix in a
 *  Rcpp::List with the corresponding InfoMap and InfoType.
 */
class RDataHandler: public DataHandlerBase<RDataHandler>
{
  public:
    typedef DataHandlerBase<RDataHandler> Base;
    typedef Base::InfoMap InfoMap;
    typedef std::map<std::string, int> InfoType;
    /** default constructor */
    inline RDataHandler(): Base(), nbSample_(0), nbVariable_(0) {}
    /** destructor */
    inline ~RDataHandler() {}
    /** @return the number of sample (the number of rows of the data) */
    inline int nbSample() const { return nbSample_;}
    /** @return the number of variables (the number of columns of the data) */
    inline int nbVariable() const { return nbVariable_;}

    /** Add a Matrix to the existing data sets
     *  @param idData Id of the data set
     *  @param idModel Id of the model to use
     *  @param data a R list with the numeric data sets
     *  @note cannot be passed as const& due to a bug from the Rcpp side (fixed
     *  in the next release of Rcpp).
     **/
    template<int Rtype>
    inline void addData( Rcpp::Matrix<Rtype> const& data, std::string idData, std::string const& idModel )
    {
      if (addInfo(idData, idModel))
      {
        data_.push_back(data, idData);
        addType(idData, Rtype);
        nbVariable_ += data.cols();
      }
    }
    /** return in an CArray the copied data with the given idData */
    template<typename Type>
    void getData(std::string const& idData, CArray<Type>& data, int& nbVariable) const
    {
      enum
      {
        Rtype_ = hidden::RcppTraits<Type>::Rtype_
      };
      Rcpp::Matrix<Rtype_> Rdata = data_[idData];
      RMatrix<Type> aux(Rdata);
      data = aux;
      nbVariable = data.sizeCols();
    }

  private:
    /** @brief Add the Rtype associated to a data set identified by an Id.
     *  @param idData can be any string given by the user
     *  @param Rtype represent the id of a Rtype (as defined in Rinternal.h)
     *
     *  @return @c false if there exists already idData associated to an other
     *  Rtype, @c true otherwise.
     **/
    bool addType(std::string const& idData, int Rtype);
    /** List of R data */
    Rcpp::List data_;
    /** Store the Rtype of a given data set by pair:
     * - idData: an arbitrary idData for a model,
     * - Rtype: an integer (defined in Rinternal.h) giving the type of the data
     **/
    InfoType infoType_;
    /** Number of sample */
    int nbSample_;
    /** Number of variable */
    int nbVariable_;
};

} // namespace STK

#endif /* STK_RDATAHANDLER_H */
