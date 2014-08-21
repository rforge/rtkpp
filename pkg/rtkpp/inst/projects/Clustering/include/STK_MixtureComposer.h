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
 * Project:  stkpp::Clustering
 * created on: 14 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_MixtureComposer.h
 *  @brief In this file we define the class MixtureComposer.
 **/

#ifndef STK_MIXTURECOMPOSER_H
#define STK_MIXTURECOMPOSER_H

#include <vector>
#include "STK_IMixtureComposer.h"

namespace STK
{

class IMixture;
/** @ingroup Clustering
 *  Main class for handling composed mixture models.
 *  A composed mixture model on some composed space
 *  \f$  \mathbb{X} = \subset \mathcal{X}^1\times \ldots  \times \mathcal{X}^L \f$
 *  is a density of the form
 * \f[
 *     f(\mathbf{x}|\boldsymbol{\theta})
 *     = \sum_{k=1}^K p_k \prod_{l=1}^L f^l(\mathbf{x}^l;\boldsymbol{\lambda}^l_k,\boldsymbol{\alpha}^l)
 *     \quad \mathbf{x} \in \mathbb{X}.
 * \f]
 * The \f$ p_k > 0\f$ with  \f$ \sum_{k=1}^K p_k =1\f$ are the mixing proportions.
 * The density \e f is called the component of the model. The parameters
 * \f$\boldsymbol{\lambda}^l_k, \, k=1,\ldots K \f$ are the cluster specific parameters
 * and the parameters \f$ \boldsymbol{\alpha}^l \f$ are the shared parameters.
 * */
class MixtureComposer : public IMixtureComposer
{
  public:
    typedef std::vector<IMixture*>::const_iterator ConstMixtIterator;
    typedef std::vector<IMixture*>::iterator MixtIterator;
    /** Constructor.
     * @param nbCluster,nbSample,nbVariable number of clusters, samples and Variables
     */
    MixtureComposer( int nbSample, int nbVariable, int nbCluster);
    /** copy constructor.
     *  @param composer the composer to copy
     */
    MixtureComposer(MixtureComposer const& composer);
    /** The registered mixtures will be deleted there.*/
    virtual ~MixtureComposer();

    /** @return a constant reference on the vector of mixture */
    inline std::vector<IMixture*> const& v_mixtures() const { return v_mixtures_;}
    /** Utility lookup function allowing to find a Mixture from its idData
     *  @param idData the id name of the data
     *  @return a pointer on the mixture, NULL if the mixture is not found
     **/
     IMixture* getMixture( String const& idData) const;

    /** Create a composer, but reinitialize the mixtures parameters. */
    virtual MixtureComposer* create() const;
    /** Create a clone of the current model, with mixtures parameters preserved. */
    virtual MixtureComposer* clone() const;

    /** initialize randomly the parameters of the components of the model */
    virtual void randomInit();
    /** Compute the proportions and the model parameters given the current tik
     *  mixture parameters.
     **/
    virtual void mStep();
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i index of the sample
     *  @param k index of the component
     **/
    virtual Real lnComponentProbability(int i, int k);
    /** write the parameters of the model in the stream os. */
    virtual void writeParameters(ostream& os) const;
    /** @brief compute the number of free parameters of the model.
     *  lookup on the mixtures and sum the nbFreeParameter.
     **/
    virtual int computeNbFreeParameters() const;
    /**@brief This step can be used by developer to initialize any thing which
     * is not the model. It will be called before running the estimation
     * algorithm. In this class, the @c initializeSterp method
     *  - set the number of free parameters using the pure virtual function @¢ computeNbFreeParameters()
     *  - Compute the
     *  - call the initializeStep() of all the mixtures,
     *  First initialization of the parameters of the model.
     *  This method is called in order to initialize the parameters.
     *
     **/
    virtual void initializeStep();
    /** @brief Impute the missing values.
     **/
    virtual void imputationStep();
    /** @brief Simulation of all the latent variables and/or missing data
     *  excluding class labels.
     */
    virtual void samplingStep();
    /**@brief This step can be used to signal to the mixtures that they must
     * store results. This is usually called after a burn-in phase.
     **/
    virtual void storeIntermediateResults(int iteration);
    /**@brief This step can be used by developer to finalize any thing. It will
     *  be called only once after we finish running the estimation algorithm.
     **/
    virtual void finalizeStep();
    /** @brief register a mixture to the composer.
     *  When a mixture is registered, the composer:
     *  - assign composer pointer (itself) to the mixture
     *  - add it to v_mixtures_
     * @param p_mixture a pointer on the mixture
     **/
    void registerMixture(IMixture* p_mixture);
    /** Utility method allowing to create all the mixtures using the DataHandler
     *  info of the manager.
     **/
    template<class MixtureManager>
    void createMixtures(MixtureManager& manager)
    { manager.createMixtures(*this, nbCluster());}
    /** Utility method allowing to create a specific mixture and register it.
     *  @param manager the manager with the responsibility of the creation
     *  @param idModel the id of the mixture we want to create
     *  @param idData the id name of the data
     **/
    template<class MixtureManager>
    void createMixture(MixtureManager& manager, Clust::Mixture idModel, String const& idData)
    {
      IMixture* p_mixture = manager.createMixture(idModel, idData, nbCluster());
      if (p_mixture) registerMixture(p_mixture);
    }
    /** Utility method allowing to create a mixture with a given data set
     *  and register it. The Mixture Manager will find the associated model
     *  to use with this data set.
     *  @param manager the manager with the responsibility of the creation.
     *  @param idData the id name of the data to modelize.
     **/
    template<class MixtureManager>
    void createMixture(MixtureManager& manager, String const& idData)
    {
      IMixture* p_mixture = manager.createMixture( idData, nbCluster());
      if (p_mixture) registerMixture(p_mixture);
    }
    /** Utility method allowing to get the parameters of a specific mixture.
     *  @param manager the manager with the responsibility of the parameters
     *  @param idData the Id of the data we want the parameters
     *  @param data the structure which will receive the data
     **/
    template<class ParametersManager, class Param>
    void getParameters(ParametersManager const& manager, String const& idData, Param& data)
    { manager.getParameters(getMixture(idData), idData, data);}

  protected:
    /** @brief Create the composer using existing data handler and mixtures.
     * This method is essentially used by the create() method and can be
     * reused in derived classes.
     * @sa MixtureComposerFixedProp
     **/
    void createComposer( std::vector<IMixture*> const& v_mixtures_);
  private:
    /** vector of pointers to the mixtures components */
    std::vector<IMixture*> v_mixtures_;
};

/** @brief specialization of the composer for the fixed proportion case.
 **/
class MixtureComposerFixedProp : public MixtureComposer
{
  public:
    /** Constructor.
     * @param nbCluster,nbSample,nbVariable number of clusters, samples and Variables
     */
    inline MixtureComposerFixedProp( int nbSample, int nbVariable, int nbCluster)
                                   : MixtureComposer( nbSample, nbVariable, nbCluster) {}
    /** copy constructor.
     *  @param model the model to copy
     */
    inline MixtureComposerFixedProp( MixtureComposer const& model)
                                   : MixtureComposer(model) {}
    /** destructor */
    inline virtual ~MixtureComposerFixedProp() {}
    /** Create a composer, but reinitialize the mixtures parameters. */
    virtual MixtureComposerFixedProp* create() const;
    /** Create a clone of the current model, with mixtures parameters preserved. */
    virtual MixtureComposerFixedProp* clone() const;
    /** @return the number of free parameters of the mixture */
    int computeNbFreeParameters() const;
    /** overloading of the computePropotions() method.
     * Let them initialized to 1/K. */
    virtual void pStep();
};

} /* namespace STK */

#endif /* STK_MIXTURECOMPOSER_H */
