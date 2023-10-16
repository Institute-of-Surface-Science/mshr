// Copyright (C) 2017 Benjamin Kehlet
//
// This file is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//

#include <dolfin/mesh/Mesh.h>
#include <dolfin/generation/RectangleMesh.h>
#include <mshr/CSGGeometries3D.h>

int main(int /* argc */, char** /* argv */)
{
  std::shared_ptr<dolfin::Mesh> m =
    std::make_shared<dolfin::RectangleMesh>( dolfin::Point(1., 2.), dolfin::Point(2., 4.),
					     5, 6);

  std::shared_ptr<mshr::CSGGeometry> g = mshr::CSGGeometries::import_mesh(m);
      
  return 0;
}
