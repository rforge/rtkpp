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
 * Project:  stkpp::STatistiK
 * Purpose:  Compute factors of a set of variables.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_Factor.h
 *  @brief In this file we define and implement the Factor class.
 **/

#ifndef STK_STAT_FACTOR_H
#define STK_STAT_FACTOR_H

#include <Sdk/include/STK_IRunner.h>

namespace STK
{
namespace Stat
{
/** @ingroup StatDesc
 *  @brief Computation of the Factors of a 2D Container.
 *
 *  The class @c Factor is a factory class for computing the factors of an array.
 *  The values can be of any type. the Coding is performed from the previous type
 *  in integer. The mapping is stored and preserved in an array of map array.
 **/
template <class Array>
class Factor: public IRunnerWithData<Array>
{
  public:
    typedef IRunnerWithData<Array> Base;
    typedef typename Array::Row RowVector;
    typedef typename Array::Col ColVector;
    typedef typename Array::Type Type;

    typedef std::map<Type, int> EncodingMap;
    typedef std::map<int, Type> DecodingMap;

    typedef CArrayPoint<EncodingMap> Encoder;
    typedef CArrayPoint<DecodingMap> Decoder;

    using Base::p_data_;

    /** Default Constructor. */
    Factor();
    /** Constructor.
     *  @param data a reference on the data set
     **/
    Factor( Array const& data);
    /** Constructor.
     *  @param p_data a pointer on the data set
     **/
    Factor( Array const* p_data);
    /** copy constructor.
     *  @param f the Factor to copy
     **/
    Factor( Factor const& f);
    /** virtual destructor.*/
    inline virtual ~Factor() {}

    /** clone pattern */
    inline virtual Factor* clone() const { return new Factor(*this);}

    /** @return array with the factors encoded as integers */
    inline CArrayXXi const& asInteger() const { return asInteger_;}
    /** @return array with for each variables the levels */
    inline CArrayPoint< Array2DVector<Type> > const& levels() const {return levels_;}
    /** @return array with for each variables the counts of the levels */
    inline CArrayPoint< VectorXi > const& counts() const {return counts_;}
    /** @return array with for each variables the number of levels */
    inline CPointXi const& nbLevels() const { return nbLevels_;}
    /** @return array with the encoding maps factor to int */
    inline Encoder const& encoder() const { return encoder_;}
    /** @return array with the encoding maps factor to int */
    inline Decoder const& decoder() const { return decoder_;}

    /** run the estimation of the Factor statistics. **/
    virtual bool run();

  protected:
    /** Array of the data size with the levels of each variables in an integer format*/
    CArrayXXi asInteger_;
    /** Number of levels of each variables */
    CPointXi nbLevels_;
    /** Array with the levels of each variables */
    CArrayPoint< Array2DVector<Type> > levels_;
    /** Array with the counts of each factor */
    CArrayPoint< VectorXi > counts_;
    /** encoder of the levels */
    Encoder encoder_;
    /** decoder of the levels */
    Decoder decoder_;

    /** udpating method in case we set a new data set */
    virtual void update();
};

template <class Array>
Factor<Array>::Factor(): Base(), asInteger_(), nbLevels_(), levels_(), counts_(), encoder_() {}

template <class Array>
Factor<Array>::Factor( Array const& data): Base(data)
                                         , asInteger_(p_data_->rows(),p_data_->cols())
                                         , nbLevels_(p_data_->cols(), baseIdx - 1)
                                         , levels_(p_data_->cols())
                                         , counts_(p_data_->cols())
                                         , encoder_(p_data_->cols())
                                         , decoder_(p_data_->cols())
{}

/* Constructor.
 *  @param p_data a pointer on the data set
 **/
template <class Array>
Factor<Array>::Factor( Array const* p_data): Base(p_data)
                                           , asInteger_()
                                           , nbLevels_()
                                           , levels_()
                                           , encoder_()
                                           , decoder_()
{
  if (p_data_)
  {
    asInteger_.resize(p_data_->rows(),p_data_->cols());
    nbLevels_.resize(p_data_->cols()).setValue(baseIdx - 1);
    levels_.resize(p_data_->cols());
    counts_.resize(p_data_->cols());
    encoder_.resize(p_data_->cols());
    decoder_.resize(p_data_->cols());
  }
}

template <class Array>
Factor<Array>::Factor( Factor const& f): Base(f), asInteger_(f.asInteger_), nbLevels_(f.nbLevels_)
                                       , levels_(f.levels_), counts_(f.counts_)
                                       , encoder_(f.encoder_)
                                       , decoder_(f.decoder_)
{}

template <class Array>
void Factor<Array>::update()
{
   // if there is no data there is nothing to update
   if (p_data_)
   {
     asInteger_.resize(p_data_->rows(),p_data_->cols());
     nbLevels_.resize(p_data_->cols()).setValue(baseIdx - 1);
     for(int j=encoder_.begin(); j<std::min(encoder_.end(), p_data_->cols().end()); ++j)
     {
       levels_[j].clear();
       counts_[j].clear();
       encoder_[j].clear();
       decoder_[j].clear();
     }
     levels_.resize(p_data_->cols());
     counts_.resize(p_data_->cols());
     encoder_.resize(p_data_->cols());
     decoder_.resize(p_data_->cols());
   }
}

template <class Array>
bool Factor<Array>::run()
{
  if (!p_data_)
  { this->msg_error_ = STKERROR_NO_ARG(FactorArray::run,data is not set);
    return false;
  }
  try
  {
    for (int j=p_data_->beginCols(); j< p_data_->endCols(); ++j)
    {
      for (int i=p_data_->beginRows(); i< p_data_->endRows(); ++i)
      {
        // find coding
        Type idData = p_data_->elt(i,j);
        typename EncodingMap::const_iterator it = encoder_[j].find(idData);
        if (it != encoder_[j].end()) // levels already exist, just update the levels array
        { asInteger_(i,j) = it->second;
          counts_[j][it->second]++;  // add one to this level
        }
        else // find a new level to add
        {
          // create a new level and set it
          asInteger_(i,j) = (++nbLevels_[j]);
          encoder_[j].insert(std::pair<Type, int>(idData, nbLevels_[j]));
          decoder_[j].insert(std::pair<int, Type>(nbLevels_[j], idData));
          levels_[j].push_back(idData);
          counts_[j].push_back(1); // start counting for this new level
        }
      }
    }
  }
  catch (Exception const& error)
  {
    this->msg_error_ += _T("Error in Factor::run():\nWhat: ");
    this->msg_error_ += error.error();
    return false;
  }
  // no error
  return true;
}

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_FACTOR_H */
