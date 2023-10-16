// Copyright (C) 2012 Anders Logg, Johannes Ring, 2012-2016 Benjamin Kehlet
//
// This file is part of mshr.
//
// mshr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mshr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mshr.  If not, see <http://www.gnu.org/licenses/>.

#include <mshr/MeshGenerator.h>
#include <mshr/CSGGeometry.h>
#include <mshr/CSGCGALMeshGenerator2D.h>
#include <mshr/CSGCGALMeshGenerator3D.h>
#include <mshr/TetgenMeshGenerator3D.h>

#include <mshr/CSGCGALDomain2D.h>

#include <dolfin/log/log.h>
#include <dolfin/mesh/BoundaryMesh.h>


namespace mshr
{

//-----------------------------------------------------------------------------
std::shared_ptr<dolfin::Mesh> generate_mesh(std::shared_ptr<const CSGGeometry> geometry,
                                            double resolution,
                                            std::string backend)
{
  if (resolution <= 0)
  {
    dolfin::dolfin_error("MeshGenerator.cpp",
                         "generate from CSG geometry",
                         "Resolution argument must be positive");
  }


  if (geometry->dim() == 2)
  {
    if (backend != "cgal")
    {
      const std::string e = "Unknown mesh generator backend: " + backend + ". The only supported 2D backend is cgal.";
      dolfin::dolfin_error("MeshGenerator.cpp",
                           "generate mesh of 2D geometry",
                           e);
    }

    CSGCGALMeshGenerator2D generator;

    // Estimate the bounding circle radius from the geometry
    const std::pair<dolfin::Point, double> bounding_circle = geometry->estimate_bounding_sphere();
    const double segment_granularity = 2*bounding_circle.second/resolution;
    const double cell_size = 2.0*bounding_circle.second/resolution;

    // Calculate surface representation
    log(dolfin::TRACE, "Request cell size: %f", cell_size);
    generator.parameters["mesh_resolution"] = -1.0;
    generator.parameters["cell_size"] = cell_size;

    std::shared_ptr<CSGCGALDomain2D> total_domain(new CSGCGALDomain2D(geometry, segment_granularity));
    std::vector<std::pair<std::size_t, std::shared_ptr<const CSGCGALDomain2D>>> subdomain_geometries;
    for (const std::pair<std::size_t, std::shared_ptr<const CSGGeometry>>& subdomain : geometry->get_subdomains())
    {
      subdomain_geometries.push_back(std::make_pair(subdomain.first,
                                                    std::shared_ptr<CSGCGALDomain2D>(new CSGCGALDomain2D(subdomain.second, segment_granularity))));
    }
    return generator.generate(total_domain, subdomain_geometries);
  }
  else if (geometry->dim() == 3)
  {
    std::shared_ptr<CSGCGALDomain3D> domain(new CSGCGALDomain3D(geometry));
    domain->ensure_meshing_preconditions();

    if (backend == "cgal")
    {
      CSGCGALMeshGenerator3D generator;
      generator.parameters["mesh_resolution"] = resolution;
      return generator.generate(std::move(domain));
    }
    else if (backend == "tetgen")
    {
      TetgenMeshGenerator3D generator;
      generator.parameters["mesh_resolution"] = resolution;
      return generator.generate(std::move(domain));
    }
    else
    {
      std::string e = "Unknown mesh generator backend: " + backend;
      dolfin::dolfin_error("MeshGenerator.cpp",
                           "Generator mesh of 3D geometry",
                           e);
      return std::shared_ptr<dolfin::Mesh>();
    }
  }
  else
  {
    dolfin::dolfin_error("MeshGenerator.cpp",
                         "create mesh from CSG geometry",
                         "Unhandled geometry dimension %d", geometry->dim());
    return std::shared_ptr<dolfin::Mesh>();
  }
}

}
