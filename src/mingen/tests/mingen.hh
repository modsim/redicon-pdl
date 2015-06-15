/*  mingen.hh  2015-06-15 MinGen Particle and its Factory
 *
 * Copyright (C) 2015 Svyatoslav Kondrat (Valiska)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <pdlib/mingen/mRNA.hh>
#include <pdlib/mingen/gene.hh>
#include <vector>

typedef enum {MINGEN_MRNA, MINGEN_GENE} MinGenParticleType;

template <class Geometry>
class MinGenParticle
{
	public:

		typedef typename Geometry::Space Space;

		MinGenParticle (const Space & x, double t, double kon, double koff, double kmbasal, double km)
			: _type(MINGEN_GENE), g (x, 0, kon, koff, kmbasal, km), number (0) {};

		MinGenParticle (const Space & x, double D, double dt, double kdeg)
			: _type(MINGEN_MRNA), m (x, D, dt, kdeg), number (0) {};

		bool move (double dt, Geometry geo)
		{
			switch (_type)
			{
				case MINGEN_MRNA: m.move(dt, geo);
					return true;

				case MINGEN_GENE: g.move (dt);
					return true;

				default:
					throw "MinGenParticle: Unknown particle type ";
			}
		};

		Space position (void) const
		{
			switch (_type)
			{
				case MINGEN_MRNA:
					return m.position();

				case MINGEN_GENE: 
					return g.position ();

				default:
					throw "MinGenParticle: Unknown particle type ";
			}
		};

		void setNumber (int n)
		{
			number = n;
		}

		// this is a public variable to tell if this particle needs to be removed 
		// by a system after 1st order reaction (mRNA needs, Gene does not)
		bool remove (void)
		{
			if (_type == MINGEN_MRNA)
				return true;

			return false;
		}; 

		MinGenParticleType type (void) const {return _type;};

		PDL::MinGen::mRNA<Geometry> * getMRNA (void)
		{
			if (_type == MINGEN_MRNA)
				return &m;
			else
				return nullptr;
		}

		PDL::MinGen::Gene<Geometry> * getGene (void)
		{
			if (_type == MINGEN_GENE)
				return &g;
			else
				return nullptr;
		}

	private:

		const MinGenParticleType _type;
		const union 
		{
			PDL::MinGen::mRNA<Geometry>  m;
			PDL::MinGen::Gene<Geometry>  g;
		};
		long int number;

};

template <class Geometry>
class MinGenFactory
{
	public:

		typedef MinGenParticle<Geometry> Particle;

		MinGenFactory (double kon, double koff, double kmbasal, double km,
			double D, double dt, double kdeg)
		: kon (kon), koff (koff), kmbasal (kmbasal), km(km),
			D(D), dt (dt) , kdeg (kdeg) {};

		MinGenParticle<Geometry> * createParticle (const typename Geometry::Space & x, MinGenParticleType type) const
		{
			switch (type)
			{
				case MINGEN_GENE: 
					return new MinGenParticle<Geometry> (x, 0, kon, koff, kmbasal, km);

				case MINGEN_MRNA:
					return new MinGenParticle<Geometry> (x, D, dt, kdeg);

				default:
					throw "MinGenParticle: Unknown particle type";
			}
		};

	private:
		const double kon, koff, kmbasal, km;
		const double D, dt, kdeg;
};

template<class Geometry>
class MinGenReaction
{
	public:
		const int order = 1;

		MinGenReaction (const MinGenFactory<Geometry> & F) : F(F) {};

		bool apply (MinGenParticle<Geometry> * p, double dt, std::vector<MinGenParticle<Geometry>*> * l)
		{
			switch (p->type())
			{
				case MINGEN_MRNA:
				{
					PDL::MinGen::mRNA<Geometry> * m = p->getMRNA ();
					if (m->degrade(dt))
					{
						std::cerr << "I am dying...." << std::endl;
						return true;
					}
					break;
				}
				case MINGEN_GENE: 
				{
					PDL::MinGen::Gene<Geometry> * g = p->getGene ();

					if (g->mRNA (dt))
					{
						MinGenParticle<Geometry> * mnew 
							= F.createParticle (g->position(), MINGEN_MRNA);
						l->push_back (mnew);
						return true;
					}
					break;
				}

				default:
					throw "MinGenReaction: Unknown particle type ";
			}
			return false;
		};
	private:
		const MinGenFactory<Geometry> & F;
};

