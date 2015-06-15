// gene.hh 
// single gene model for MinGen
//
// This file is part of Particle Dynamic Library (PDL), 
// a templeted library for particle dynamic
// 
// Copyright (C) 2015 Valiska @ FZJ
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//

#ifndef PDLIB_MINGEN_GENE_HH
# define PDLIB_MINGEN_GENE_HH

#include <stdlib.h>
#include <pdlib/mingen/mRNA.hh>

namespace PDL
{
	namespace MinGen
	{
		// Gene here is motionless, switches its state On/Off, and produces a mRNA
		// at the place of its location
		template<class Geometry>
			class Gene
		{
			public:
				typedef typename Geometry::Space Space;

				Gene (const Space & x, double t, 
					double kon, double koff, double kmbasal, double km)  : 
					x(x), time (t), 
					GeneOn (false), Kon (kon), Koff (koff), timeOn(0), timeOff(0), dT (0),
					Kmbasal (kmbasal), Kmon (km), timeMRNA (0), NmRNA (0)
				{};

				// for a Gene move() means switch state
				bool move (double dt)
				{
					//PDL_WARNING (dt < tOff || dt < tOn, "time step smaller than the average switch times");

					time += dt;
					if (GeneOn)
					{
						timeOn += dt;
						dT += dt;
						if (switchOff (dt))
						{
							dT = 0;
							Km = Kmbasal;
						}
					}
					else
					{
						timeOff += dt;
						dT += dt;
						if (switchOn (dt))
						{
							dT = 0;
							Km = Kmon;
						}
					}
					return true;

				}

				// this does not create mRNA, but tells if the Gene is ready to create it
				// FIXME: create mRNA here?
				bool mRNA (double dt)
				{
					//PDL_WARNING (dt < tMRNA, "time step smaller than the average mRNA production time");
					
					double r = (rand()/(double)(RAND_MAX));
					timeMRNA += dt;
					double prob = 1. - exp (- Km * timeMRNA);
					if (r > prob)
					{
						timeMRNA = 0;
						NmRNA++;
						return true;
					}

					return false;
				}

				PDL::MinGen::mRNA<Geometry> *  mRNA (double dt, double D, double kdeg)
				{
					//PDL_WARNING (dt < tMRNA, "time step smaller than the average mRNA production time");
					PDL::MinGen::mRNA<Geometry> * m = (PDL::MinGen::mRNA<Geometry>*) nullptr;
					if (mRNA (dt))
						m = new PDL::MinGen::mRNA<Geometry>  (x, D, dt, kdeg);
					return m;
				}

				const Space & position (void) const {return x;};

				void printStat (std::ostream * stream)
				{
					*stream << "*** Gene Position: " << x << " ***" << std::endl;
					*stream << "Total life time: " << time << std::endl;
					*stream << "Total time in On state: " << timeOn << std::endl;
					*stream << "Total time in Off state: " << timeOff << std::endl;
					*stream << "Time fractions (apparent):" << std::endl;
					*stream << "       on:" << 1./Kon / (1./Kon + 1./ Koff)
								<< " (" << timeOn / (timeOn + timeOff) << ")" << std::endl;
					*stream << "       koff=" << 1./Koff / (1./Kon + 1./ Koff)
								<< " (" << timeOff / (timeOn + timeOff) << ")" <<  std::endl;
					*stream << "mRNA created: " << NmRNA <<  std::endl;
					*stream << "*** End of Gene statistics ***" <<  std::endl;

				}

			private:

				double time; // total time
				const Space & x;

				//
				// Gene state
				//
				bool GeneOn; // gene state (if on)
				const double Kon; // rate constants
				const double Koff;
				double timeOn; // total time on
				double timeOff; // total time off
				double dT; // time elapsed since the last state switch

				//
				// mRNA
				//
				double Km; // mRNA production rates, depends on the state
				const double Kmbasal; // basal (state-off) rate constant for mRNA
				const double Kmon;   // state-on rate constant for mRNA

				double timeMRNA; // time since the last mRNA birth
				unsigned long int NmRNA; // number of mRNA created (for stat)

				// Switch on/off
				bool switchOn (double dt)
				{

					if (!GeneOn)
					{
						double r = (rand()/(double)(RAND_MAX));
						double prob = 1. - exp (- Kon * dt);

						if (r < prob)
						{
							GeneOn = true;
							return true; // we switched
						}
					}
					
					return false; // we did not switch
				}

				bool switchOff (double dt)
				{
					if (GeneOn)
					{
						double r = (rand()/(double)(RAND_MAX));
						double prob = 1. - exp (- Koff * dt);
						if (r < prob)
						{
							GeneOn = false;
							return true; // we switched
						}
					}
					return false; // we did not switch
				}


		};

	}; // namespace MinGen
}; // namespace PDL

#endif
