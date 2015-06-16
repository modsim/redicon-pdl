/* test MinGen class  2015-06-15
 *
 * Units: [D] = mu2/s, [k] = 1/s
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
	H[0] = 2.; H[1] = 2.;

	typedef PDL::GeometryBox<double, 2> Box;
	Box b (x0,H);

	Dune::FieldVector<double, 2> x;
	x[0] = 1.; x[1] = 1.;
	
	double dt = .1;

	// Gene
	// Kon = Koff = 0.16(6) 1/s (Nat Gen 6, 455, 2005)
	// for graded: half-life ~ 0.5/hour ~ 0.001 (1)/s
	double kon  = 0.0001;
	double koff = 0.0001;

	// Km = 50/min ~ 0.83(3) /s, Kmbasal = 5/min (Nat Gen 6, 455, 2005)
	double kmbasal = 0.02;
	double km = 0.2;

	// mRNA
	// D = 0.01 - 0.03 mum2/s (Nat Comm 4 2414, 2013)
	double D = 0.01;
	// Kdeg = 0.1/min ~ 0.0016(6)/s (Nat Gen 6, 455, 2005)
	double kdeg = 0.001;

	typedef MinGenFactory<Box> Factory;
	Factory F (kon, koff, kmbasal, km, D, dt, kdeg);
	
	typedef MinGenReaction<Box> Rxn;
	PDL::System<Box, Factory, Rxn> system (b, F);

	Rxn rxn (F);
	system.addReaction (rxn);

	system.addParticle (x, MINGEN_GENE);

	std::ofstream stream;
//	stream.open("mRNA.dat", std::ios_base::app);
	stream.open("mRNA.dat");

	int tfinal = 20*60*60;
	int N = rint (tfinal/dt);
	int skip = 60;
	int iskip = 1;
	for (int i =0; i < N; i++)
	{
		system.evolve (dt);

		stream << system.time() << "     " << system.getNParticles() - 1 << std::endl;
	
		if ( !(iskip < skip) )
		{
			std::cerr << "*** t=" << system.time()/60 << " (" << system.getNParticles() - 1 
				<< " mRNA) ***" << std::endl;
			iskip = 0;
		}
		iskip++;
	}

	// Check gene statistics
	MinGenParticle<Box> * p = system.getParticle (0);
	PDL::MinGen::Gene<Box> * g = p->getGene();
	g->printStat (&std::cerr);

	return 1;
}

