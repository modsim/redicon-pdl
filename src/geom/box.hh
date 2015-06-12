// box.hh 
// rectangular box for system.hh
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

#ifndef PDLIB_GEO_BOX_HH
# define PDLIB_GEO_BOX_HH

#include <dune/common/fvector.hh>

namespace PDL 
{
	template<typename ctype, int dim>
		class GeometryBox
		{
			public:
				enum {dimension = dim};
				typedef ctype type;
				typedef typename Dune::FieldVector<ctype,dim>  Space;

				GeometryBox (const Space & x0, const Space & H)
					: x0(x0), H(H)
				{
					for (int i = 0; i < dimension; i++)
					{
						xl[i] = x0[i] - 0.5 * H[i];
						xr[i] = x0[i] + 0.5 * H[i];
					}
#ifdef DEBUG
					std::cerr << "left corner: " << xl << std::endl;
					std::cerr << "right corner: " << xr << std::endl;
#endif
					srand(time(NULL)); // Seed the time

				};

				Space randomPoint () 
				{
					Space x;
					for (int i = 0; i < dim; i++)
					{
						x[i] = xl[i] + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(H[i])));
					}
					return x;
				}

				bool inside (const Space & x) const
				{
#ifdef DEBUG
					std::cerr << "x= " << x << std::endl;
					std::cerr << "left corner: " << xl << std::endl;
					std::cerr << "right corner: " << xr << std::endl;
#endif
					for (int i = 0; i < dimension; i++)
						if ( (x[i] > xr[i]) || (x[i] < xl[i]) )
							return false;
					return true;
				};
				double toBoundary (const Space & x, const Space & n) const
				{
					double t = 0;
					double norm = n.two_norm();
#ifdef DEBUG
					std::cerr << "norm = " << norm << std::endl;
#endif
					for (int i = 0; i < dimension; i++)
					{
						if (n[i] != 0.)
						{
							double tl = norm * (xl[i] - x[i])  / n[i];
							std::cerr << "i=" << i << " t=" << t
								<< " tl=" << tl 
								<< " xl=" << xl[i] << " x=" << x[i] << std::endl;
							if (tl > 0.) 
							{
								if ((t == 0.0) || (tl < t))
									t = tl;
							}

							double tr = norm * (xr[i] - x[i])  / n[i];
							std::cerr << "i=" << i << " t=" << t
								<< " tr=" << tr << " xr=" 
								<< xr[i] << " x=" << x[i] << std::endl;
							if (tr > 0.) 
							{
								if ((t == 0.0) || (tr < t))
									t = tr;
							}
						}
					}

					return t;
				};

			private:

				const Space & x0;
				const Space & H;

				Space xl;
				Space xr;
		};

}; // namespace PDL

#endif

