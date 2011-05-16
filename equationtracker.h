#ifndef _EQUATIONTRACKER_H
#define _EQUATIONTRACKER_H

#include "threevector.h"
#include "basetracking.h"
#include "basegeometry.h"

/**
 * @class EquationTracker
 * Calculate the track of a particle by solving the equation of motion.
 * This class is an abstract base class, implementations have to overwrite
 * the derivs method.
 */
class EquationTracker : public Basetracking {
	public:
		EquationTracker(Random *ran, Basegeometry *geo);
		void makeTrack(double t_start, double h);
		void initialize();
		Threevector getPosition(double time);

	private:
		void rkStep(const double &h, const double &t, Threevector &x, Threevector &v);
		double fStepSize; ///< stepsize of the Runge-Kutta integration

		// State for makeTrack()
		double fTime;
		Threevector fPos;
		Threevector fVel;
	
	protected:
		/**
		 * Provides the derivatives @p xdot and @p vdot of time and velocity.
		 * Override this pure virtual method to provide motional equations.
		 *
		 * @attention The parameters @p xdot and @p vdot are not initialized
		 * to zero, so you have to overwrite all components!
		 *
		 * @param[in]  t    current time
		 * @param[in]  x    current location
		 * @param[in]  v    current velocity
		 * @param[out] xdot \f$\dot{x}\f$, usually just \f$\dot{x} = v\f$, <b>not initialized to zero</b>
		 * @param[out] vdot \f$\dot{v}\f$, usually \f$\dot{v} = \frac{F}{m}\f$ (Newton), <b>not initialized to zero</b>
		 */
		virtual void derivs(const double t, const Threevector x, const Threevector v, Threevector &xdot, Threevector &vdot) = 0;
	};

#endif // _EQUATIONTRACKER_H
