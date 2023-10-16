// Copyright (C) 2012-2016 Benjamin Kehlet
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
//
// Modified by Johannes Ring, 2012

#ifndef __MSHR_CGAL_MESH_GENERATOR2D_H
#define __MSHR_CGAL_MESH_GENERATOR2D_H

#include <dolfin/common/Variable.h>

#include <vector>
#include <memory>

// Forward declaration
namespace dolfin{ class Mesh; }

namespace mshr
{

  // Forward declaration
  class CSGCGALDomain2D;

  /// @brief Mesh generator for Constructive Solid Geometry (CSG)
  /// utilizing CGALs 2D Regularized Boolean Set-Operations
  class CSGCGALMeshGenerator2D : public dolfin::Variable
  {
  public :

    /// @brief Create mesh generator
    CSGCGALMeshGenerator2D();

    /// @brief Destructor
    ~CSGCGALMeshGenerator2D();

    /// @brief Generate mesh
    /// @param geometry The geometry to be meshed
    /// @param mesh The Dolfin mesh object to be returned. Will ble cleared.
    /// @param partition (experimental) If true then generate mesh on process 0 and partition and distribute
    ///   it. If false, then a full mesh is generated on all processes.(only
    ///   relevant when running in parallel)
    std::shared_ptr<dolfin::Mesh> generate(std::shared_ptr<const CSGCGALDomain2D> domain,
                                           const std::vector<std::pair<std::size_t, std::shared_ptr<const CSGCGALDomain2D>>>& subdomains
                                           = std::vector<std::pair<std::size_t, std::shared_ptr<const CSGCGALDomain2D>>>());

    /// Default parameter values
    static dolfin::Parameters default_parameters()
    {
      dolfin::Parameters p("csg_cgal_meshgenerator");

      // DEPRECATED! Use cell size
      p.add("mesh_resolution", 64.0);
      p.add("triangle_shape_bound", 0.125);
      p.add("cell_size", 0.25);
      p.add("lloyd_optimize", false);

      // shorter edges in the domain will be collapsed before meshing
      p.add("edge_truncate_tolerance", 1e-16);
      p.add("partition", true);

      return p;
    }
  };

}

#endif
