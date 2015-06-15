// system.hh 
// system of particles and reactions
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

#ifndef PDLIB_SYSTEM_HH
# define PDLIB_SYSTEM_HH

#include <vector>
#include <string>
#include <fstream>

namespace PDL 
{

template<class P>
class null_reaction
{
	public:
		const int order = -1;
		bool apply (const P * , double dt, std::vector<P*> *) {return false;}
		//{std::vector<P*> l; return l;};
};

template<class Geometry, class Factory, class Reaction = null_reaction<typename Factory::Particle>>
class System
{
	public:

		typedef typename Factory::Particle Particle;

		System (Geometry & G, Factory & F) : geo(G), F(F), number (0), t (0)
		{
//			static_assert (typename Particle::Space == typename Geometry::Space, "Particle and geometry inconsistent");
		};

		~System () 
		{
			for (int i = getNParticles() - 1 ; i >= 0; i--)
			{
#ifdef DEBUG
				std::cerr << "removing particle " << i << std::endl;
#endif				
				delParticle (i);
			}
		}

		// type is `abstract' particle type to make it possible to create
		// particles of different types (types depend on the implementation)
		template<typename... Args>
		bool addParticle (const typename Geometry::Space & x, Args... args)
		{
			Particle * p = nullptr;
			if (geo.inside (x))
			{
				p = F.createParticle(x, args...);
				return addParticle (p);
			}

			return false;
		}
		template<typename... Args>
		bool addParticle (Args... args)
		{
			typename Geometry::Space x = geo.randomPoint();
			Particle * p = F.createParticle(x, args...);
			return addParticle (p);
		}

		bool addParticle (Particle * p)
		{
			p->setNumber(number++); // to keep track of particles
			plist.push_back(p);
			return true;
		}


		void delParticle (int number)
		{
			delete plist.at (number);
			//std::cerr << "particle # " << number << " of " << getNParticles() << std::endl;
			plist.erase (plist.begin() + number);
			//std::cerr << "particle # " << number << " of " << getNParticles() << std::endl;
		}

		Particle * getParticle (int i) const {return plist.at (i);};
		int getNParticles () const {return plist.size ();};

		typename Geometry::Space particlePosition (int i) 
		{
			Particle * p = plist.at(i); 
			return p->position(); 
		};

		void print (const std::string & fname) 
		{
			std::ofstream * stream = new std::ofstream (fname);
			for (typename std::vector<Particle*>::iterator ps = plist.begin(); ps != plist.end(); ++ps)
			{
				(*ps)->print (stream);
			}
			delete stream;
		}

		void printPositions (const std::string & fname) 
		{	
			std::ostream * stream = new std::ofstream (fname);
			for (typename std::vector<Particle*>::iterator ps = plist.begin(); ps != plist.end(); ++ps)
			{
					*stream << (*ps)->position () << std::endl;
			}
			delete stream;

		};

		// add a reaction: no need to remove
		// standard calls: 
		void addReaction (const Reaction & r)
		{
			return rxnlist.push_back(r);
		}


		bool evolve (const double dt)
		{
			// first check reactions
			for (typename std::vector<Reaction>::iterator rxn = rxnlist.begin(); rxn != rxnlist.end(); ++rxn)
			{
				if ((*rxn).order == 1)
				{
					std::vector<Particle*> newplist;

					for (int i = 0; i < getNParticles(); i++)
					{
						Particle * p = getParticle (i);
						if ( (*rxn).apply (p, dt, &newplist) )
						{
#ifdef DEBUG						
							std::cerr << "Removing particle " << i << std::endl;
#endif							
							if (p->remove())
								delParticle (i); // this is particle number
						}
					}
					for (typename std::vector<Particle*>::iterator newp = newplist.begin(); newp != newplist.end(); ++newp)
					{
						addParticle (*newp);
//						plist.push_back(*newp);
					}
				}
			}

			// now make a dynamic move
			for (typename std::vector<Particle*>::iterator ps = plist.begin(); ps != plist.end(); ++ps)
			{
			//	(*ps)->move (dt, geo, plist);
			//		return false;
				if (!(*ps)->move (dt, geo))
					return false;
			}
			t += dt;
			return true;
		};

		std::vector<Particle*> particleList () {return plist;};

		double time () const {return t;};
	private:

		Geometry & geo;
		Factory & F;
		std::vector<Particle*> plist;
		std::vector<Reaction> rxnlist;
		int number; // total number of particles
		double t; // current time

};

}; // End of namespace PDL

#endif
