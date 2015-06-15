/*  test MinGen class  2015-06-15
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

#include <pdlib/system.hh>
#include <pdlib/geom/box.hh>

#include "mingen.hh"

int main (int argc, char ** argv) 
{
	Dune::FieldVector<double, 2> x0;
	x0[0] = 0.; x0[1] = 0.;

	Dune::FieldVector<double, 2> H;
	H[0] = 10.; H[1] = 10.;

	typedef PDL::GeometryBox<double, 2> Box;
	Box b (x0,H);

	Dune::FieldVector<double, 2> x;
	x[0] = 5.; x[1] = 5.;
	
	double dt = 0.01;

	// Gene
	double kon = 5.;
	double koff = 5.;
	double kmbasal = 5;
	double km = 25;

	// mRNA 
	double D = 0.1;
	double kdeg = 1.;

	typedef MinGenFactory<Box> Factory;
	Factory F (kon, koff, kmbasal, km, D, dt, kdeg);
	
	typedef MinGenReaction<Box> Rxn;
	PDL::System<Box, Factory, Rxn> system (b, F);

	Rxn rxn (F);
	system.addReaction (rxn);

	system.addParticle (x, MINGEN_GENE);

	std::cout << "Number of particles " << system.getNParticles() << std::endl;
	for (int i =0; i < 10000; i++)
	{
		system.evolve (dt);
		//system.printPositions();
		
		std::vector<Factory::Particle*> plist = system.particleList (); 
		std::cerr << std::endl;
		std::cerr << "*** t=" << system.time() << " (" << system.getNParticles() - 1 << " mRNA) ***" << std::endl;
		for (typename std::vector<Factory::Particle*>::iterator ps = plist.begin(); ps != plist.end(); ++ps)
		{
			std::cerr << (*ps)->position() << std::endl;
		}
	}

	// Check gene statistics
	MinGenParticle<Box> * p = system.getParticle (0);
	PDL::MinGen::Gene<Box> * g = p->getGene();
	g->printStat (&std::cerr);

	return 1;
}

