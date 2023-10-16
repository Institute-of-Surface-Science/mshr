// Copyright (C) 2014-2016 Benjamin Kehlet
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

#ifndef __MSHR_DOLFIN_MESH_UTILS_H
#define __MSHR_DOLFIN_MESH_UTILS_H

#include <utility>
#include <memory>

namespace dolfin { class Mesh; }

namespace mshr
{

class DolfinMeshUtils
{
 public:
  /// Compute the smallest and largest cell wrt. volume.
  /// @param m The mesh
  static std::pair<double, double> cell_volume_min_max(const dolfin::Mesh& m);

  /// Check that all vertices has at least one incident cell
  /// @param m The mesh
  static bool has_isolated_vertices(const dolfin::Mesh& m);

  /// Run all implemented checks mesh consistency
  /// @param m The mesh
  static bool check_mesh(const dolfin::Mesh& m);

  static std::shared_ptr<dolfin::Mesh>
    extract_subdomain(std::shared_ptr<const dolfin::Mesh>,
                      std::size_t cell_domain);

  /// @brief Glue together two meshes by matching boundary vertices.
  /// @param m1 Mesh 1 to merge
  /// @param m2 Mesh 2 to merge
  /// @param m1_boundary_marker Mark boundary facets the belongs to m1. A value < 0 disables marking.
  /// @param m1_boundary_marker Mark boundary facets the belongs to m2. A value < 0 disables marking.
  /// @param m1_boundary_marker Mark boundary facets the belongs to the interface between m1 an m2. A value < 0 disables marking.
  static std::shared_ptr<dolfin::Mesh>
    merge_meshes(std::shared_ptr<dolfin::Mesh> m1,
                 std::shared_ptr<dolfin::Mesh> m2,
                 int m1_marker=1,
                 int m2_marker=2,
                 int m1_boundary_marker=1,
                 int m2_boundary_marker=2,
                 int interface_marker=3);
};

}
#endif
