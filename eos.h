/**
 * \file   eos.h
 *
 * \brief  This header file contains classes for the calculation of the Equation of State (EOS)
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *
 * \version 1.0 
 *
 * \date 2014/05/09
 *
 * Created on: 2014/05/01 
 *
 */


#ifndef __EOS_H__
#define __EOS_H__

#include <vector>


/**
 * \class EOS
 * 
 * \brief An abstract class for the calculation of energy and sound speed based on different EOS models
 *
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *
 * \version 1.0 
 *
 * \date 2014/05/09 
 *
 * Created on: 2014/05/01 
 *
 */
class EOS {
protected:
	int m_iEOSChoice; ///< The eos choice: 1=Polytropic gas; 2=Stiffened Polytropic gas
public:
	/// Destructor
	virtual ~EOS() {};

	/// Getter function of the protected data member m_iEOSChoice
	int getEOSChoice() {return m_iEOSChoice;}
	
	/**
	 * \brief       Getter function of all the parameters specified in the construtor argument list 
	 * \param [out] params All the parameters of this EOS 
	 * \return      None
	 */
	virtual void getParameters(std::vector<double>& params) = 0;

	/**
	 * \brief       Calculates energy based on this EOS and the input pressure and density values  
	 * \param [in]  pressure the input pressure value
	 * \param [in]  density the input density value
	 * \return      the calculated energy value
	 */
	virtual double getEnergy(double pressure, double density) = 0;
	
	/**
	 * \brief       Calculates sound speed based on this EOS and the input pressure and density values  
	 * \param [in]  pressure the input pressure value
	 * \param [in]  density the input density value
	 * \return      the calculated sound speed value
	 */
	virtual double getSoundSpeed(double pressure, double density) = 0;
};










/**
 * \class PolytropicGasEOS
 * 
 * \brief Calculates energy and sound speed based on the Polytropic gas eos model
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *
 * \version 1.0 
 *
 * \date 2014/05/09 
 *
 * Created on: 2014/05/01 
 *
 */
class PolytropicGasEOS : public EOS {
protected:
	double m_fGamma; ///< The parameter \e gamma

public:
    /**
	 * \brief       Constructor
	 * \param [in]  gamma The parameter \e gamma
	 */
	PolytropicGasEOS(double gamma) : m_fGamma(gamma) {m_iEOSChoice=1;}
	
	// Destructor
	virtual ~PolytropicGasEOS() {}	

	/**
	 * \brief       Getter function of all the parameters specified in the construtor argument list 
	 * \param [out] params the params vector contains m_fGamma 
	 * \return      None
	 */
	virtual void getParameters(std::vector<double>& params) {params.push_back(m_fGamma);}

	/**
	 * \brief       Calculates energy based on the Polytropic gas eos and the input pressure and density values  
	 * \param [in]  pressure the input pressure value
	 * \param [in]  density the input density value
	 * \return      the calculated energy value
	 */
	virtual double getEnergy(double pressure, double density);

	/**
	 * \brief       Calculates sound speed based on the Polytropic gas eos and the input pressure and density values  
	 * \param [in]  pressure the input pressure value
	 * \param [in]  density the input density value
	 * \return      the calculated sound speed value
	 */
	virtual double getSoundSpeed(double pressure, double density);
};












/**
 * \class StiffPolytropicGasEOS
 * 
 * \brief Calculates energy and sound speed based on the Stiffened Polytropic gas eos model
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *
 * \version 1.0 
 *
 * \date 2014/05/09 
 *
 * Created on: 2014/05/01 
 *
 */
class StiffPolytropicGasEOS : public EOS {
protected:
	double m_fGamma; ///< The parameter \e gamma
	double m_fPinf; ///< The parameter pressure infinity 
	double m_fEinf; ///< The parameter energy infinity

public:
	/**
	 * \brief       Constructor
	 * \param [in]  gamma The parameter \e gamma
	 * \param [in]  pinf  The parameter \e pressure infinity
	 * \param [in]  einf  The parameter \e energy infinity
	 */
	StiffPolytropicGasEOS(double gamma, double pinf, double einf):
		m_fGamma(gamma), m_fPinf(pinf), m_fEinf(einf) {m_iEOSChoice=2;}
	
	/// Destructor
	virtual ~StiffPolytropicGasEOS() {}
	
	/**
	 * \brief       Getter function of all the parameters specified in the construtor argument list 
	 * \param [out] params the params vector contains m_fGamma, m_fPinf, and m_fEinf 
	 * \return      None
	 */
	virtual void getParameters(std::vector<double>& params) { 
		params.push_back(m_fGamma);
		params.push_back(m_fPinf);
		params.push_back(m_fEinf);
	}

	/**
	 * \brief       Calculates energy based on the Stiffened Polytropic gas eos 
	 *              and the input pressure and density values  
	 * \param [in]  pressure the input pressure value
	 * \param [in]  density the input density value
	 * \return      the calculated energy value
	 */
	virtual double getEnergy(double pressure, double density);
	
	/**
	 * \brief       Calculates sound speed based on the Stiffened Polytropic gas eos 
	 *              and the input pressure and density values  
	 * \param [in]  pressure the input pressure value
	 * \param [in]  density the input density value
	 * \return      the calculated sound speed value
	 */
	virtual double getSoundSpeed(double pressure, double density);
};



#endif // __EOS_H__
