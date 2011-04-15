#include "equationtracker.h"
#include "basegeometry.h"

class GravitationTracker : public EquationTracker {
	public:
		/**
		 * Construct a new GravitationTracker and set the acceleration to @p g.
		 * @param g acceleration due to gravitation, e.g. g = 9.81
		 */
		explicit GravitationTracker(Random *ran, Basegeometry *geo, const double g) :
			EquationTracker(ran, geo), minus_g(-g)
		{}

		/**
		 * This class represents the equation of motion for a particle in the
		 * gravitational field \f$ F = -m \cdot g \cdot \hat{e}_z \f$.
		 *
		 * The motion of equation
		 * \f[ \ddot{\mathbf{x}} = -g\,\hat{\mathbf{e}}_{z} \f]
		 *
		 * is transformed to
		 *
		 * \f{eqnarray*}{
		 * \dot{\mathbf{x}} & = & \mathbf{v}\\
         * \dot{\mathbf{v}} & = & -g\,\hat{\mathbf{e}}_{z}
		 * \f}
		 */
		inline void derivs(const double t, const Threevector x, const Threevector v, Threevector &xdot, Threevector &vdot) {
			/// The derivative @p xdot is the velocity and hence set to @p v.	
			xdot = v;

			/// Newton says that @p vdot is equal to the acceleration thich is g in z direction here.
			/// The x and y components must explicitly be set to zero because @p xdot and @p vdot are not
			/// zeroed before they are given to derivs.
			vdot[0] = 0;
			vdot[1] = 0;
			vdot[2] = minus_g;
		}

	private:
		double minus_g; ///< Acceleration \f$\dot{v}=-g\f$

		GravitationTracker(); // No default constructor
};
