/* testBox.cc  2014-10-08 test GeometryBox class
 *
 * Copyright (C) 2012 Svyatoslav Kondrat (Valiska)
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

#include <pdlib/geom/box.hh>

int main (int argc, char ** argv) 
{
	Dune::FieldVector<double, 2> x0;
	x0[0] = 0.; x0[1] = 0.;

	Dune::FieldVector<double, 2> H;
	H[0] = 10.; H[1] = 10.;

	PDL::GeometryBox<double, 2> b (x0,H);

	Dune::FieldVector<double, 2> x;
	x[0] = 0.; x[1] = 1.;

	std::cerr << "Point (" << x << ") isnide? " << b.inside (x) << std::endl;

	Dune::FieldVector<double, 2> n;
	n[0] = 5.; n[1] = 5.;
	n.two_norm();

	std::cerr << "distance from (" << x << ") along (" << n << "): " << b.toBoundary(x,n) << std::endl;

	return 1;
}

