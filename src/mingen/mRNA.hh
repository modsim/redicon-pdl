// mRNA.hh 
// mRNA model for MinGen
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

#ifndef PDLIB_MINGEN_MRNA_HH
# define PDLIB_MINGEN_MRNA_HH

#include <stdlib.h>
#include <pdlib/particles/pbp.hh>

namespace PDL
{
	namespace MinGen
	{
		// mRNA is both a Brownian point-particle 
		// and a 'reaction' which is a first order degradation reaction
		template<class Geometry>
			class mRNA : public BrownianParticle<Geometry>
		{
			public:

				mRNA (const typename BrownianParticle<Geometry>::Space & x, double D, double dt, double kdeg)  : 
					BrownianParticle<Geometry>::BrownianParticle (x, D, dt, 1),
					kdeg (kdeg) {};

				// This does NOT delete mRNA, but tells if it is ready to degrade by returning true
				bool degrade (double dt)
				{
					//PDL_WARNING (dt < tdeg, "time step smaller than the rate degradation constant");

					time += dt;
					double r = (rand()/(double)(RAND_MAX + 1));
					double prob = 1. - exp (- kdeg * time);
					if (r > prob)
						return true; // I will be deleted by the system as I leave the list empty
				}

			private:

				double time;
				const double kdeg;

		};

	}; // namespace MinGen
}; // namespace PDL

#endif
