// pbp.hh 
// point-like Brownian Particle
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

#ifndef PDLIB_PBP_HH
# define PDLIB_PBP_HH

#include <fstream>

#include <chrono>
#include <random>

#define MAX_ATTEMPT 100

namespace PDL
{
	template<class Geometry>
		class BrownianParticle
		{

			public:
				typedef typename Geometry::Space Space;

				BrownianParticle (const Space & x, double D, double dt)  : 
					x(x), D(D), dt(dt), s(sqrt (2. * D * dt)), t (-1), number(-1),
					engine (std::chrono::system_clock::now().time_since_epoch().count()),
					distribution (0, s)
				{
					static_assert(std::is_same<typename Geometry::type, double>::value, "only double supported");
				};

				BrownianParticle (const Space & x, double D, double dt, int type)  : 
					x(x), D(D), dt(dt), s(sqrt (2. * D * dt)), t (type), number(-1),
					engine (std::chrono::system_clock::now().time_since_epoch().count()),
					distribution (0, s)
				{
					static_assert(std::is_same<typename Geometry::type, double>::value, "only double supported");
				};

				int type () const {return t;};

				bool move (const double dt, Geometry & g, std::vector<class BrownianParticle*> plist)
				{
					return move (dt, g);
				}

				bool move (const double dt, Geometry & g)
				{
					if (dt != this->dt)
					{
						std::cerr << "Cannot move: Time step must be the same" << std::endl;
						return false;
					}
#ifdef DEBUG		
					std::cerr << "mean=" << distribution.mean() << ", std=" << distribution.stddev() << std::endl;
#endif
					Space xnew = x;
					for (int attempt = 0; attempt < MAX_ATTEMPT; attempt++)
					{
						for (int i = 0; i < Geometry::dimension; i++)
						{
							double a =  distribution (engine);
#ifdef DEBUG					
							std::cerr << "dX" << i << "=" << a << std::endl;
#endif					
							xnew[i] += a;
						}

						if (g.inside(xnew))
						{
							x = xnew;
#ifdef DEBUG
							std::cerr << "x=" << x 
								<< " xnew=" << xnew << std::endl;
#endif
							return true;
						} 
					}
					return false;
				};

				Space position () const {return x;};
				double getDt () {return dt;};

				void print (const std::string & name)
				{
					std::ofstream * stream = new std::ofstream (name);
					print (stream);
					delete stream;
				}

				void printApp (const std::string & name)
				{
					std::ofstream stream;
					stream.open(name, std::ios_base::app);
					print (&stream);
				}

				void print (std::ofstream * stream)
				{
					*stream << position() << std::endl;
				}

				void setNumber (int n) {number = n;} ;
				int getNumber (void) const { return number;} ;

			private:

				Space  x;
				double D;
				double dt;
				double s;

				const int t; // type
				int number; // particle number in a system, set by a system (FIXME: make friends?)

				std::mt19937 engine;
				std::normal_distribution<typename Geometry::type>  distribution;
		};

}; // namespace PDL

#endif
