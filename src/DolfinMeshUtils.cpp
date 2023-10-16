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

#include <mshr/DolfinMeshUtils.h>
#include "FuzzyPointLocator.h"

#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/math/basic.h>

#include <limits>

namespace mshr
{

std::pair<double, double> DolfinMeshUtils::cell_volume_min_max(const dolfin::Mesh& m)
{
  std::pair<double, double> res(std::numeric_limits<double>::max(), 0.0);

  for (dolfin::CellIterator cell(m); !cell.end(); ++cell)
  {
    const double v = cell->volume();
    res.first = std::min(res.first, v);
    res.second = std::max(res.second, v);
  }

  return res;
}
//-----------------------------------------------------------------------------
bool DolfinMeshUtils::has_isolated_vertices(const dolfin::Mesh& m)
{
  std::set<std::size_t> vertices;
  for (dolfin::CellIterator cit(m); !cit.end(); ++cit)
  {
    const unsigned int* v = cit->entities(0);
    for (std::size_t i = 0; i < cit->num_global_entities(0); i++)
    {
      vertices.insert(v[i]);
    }
  }

  bool isolated_vertices = false;
  for (std::size_t i = 0; i < m.num_vertices(); i++)
  {
    if (vertices.count(i) < 1)
    {
      log(dolfin::DBG, "Vertex %u has no incident cells", i);
      isolated_vertices = true;
    }
  }

  return isolated_vertices;
}
//-----------------------------------------------------------------------------
bool DolfinMeshUtils::check_mesh(const dolfin::Mesh& m)
{
  return !has_isolated_vertices(m);
}
//-----------------------------------------------------------------------------
std::shared_ptr<dolfin::Mesh>
   DolfinMeshUtils::extract_subdomain(std::shared_ptr<const dolfin::Mesh> mesh,
                                      std::size_t cell_domain)
{
  dolfin_assert(mesh->geometry().dim() == 3);
  dolfin_assert(mesh->topology().dim() == 3);

  // Collect all vertices incident to all marked cells
  std::map<std::size_t, std::size_t> collected_vertices;
  std::size_t num_cells = 0;
  for (const std::pair<std::size_t, std::size_t>& marker : mesh->domains().markers(3))
  {
    if (marker.second == cell_domain)
    {
      num_cells++;
      dolfin::Cell c(*mesh, marker.second);
      for (std::size_t i = 0; i < 4; i++)
      {
        const std::size_t s = collected_vertices.size();
        collected_vertices.insert(std::make_pair(c.entities(0)[i], s));
      }
    }
  }

  std::shared_ptr<dolfin::Mesh> outmesh(new dolfin::Mesh);
  dolfin::MeshEditor editor;
  editor.open(*outmesh, dolfin::CellType::Type::tetrahedron, 3,3);

  editor.init_vertices(collected_vertices.size());
  for (std::pair<std::size_t, std::size_t> v : collected_vertices)
  {
    dolfin::Vertex existing_vertex(*mesh, v.first);
    editor.add_vertex(v.second, existing_vertex.point());
  }

  editor.init_cells(num_cells);
  std::size_t cell_counter = 0;
  for (const std::pair<std::size_t, std::size_t>& marker : mesh->domains().markers(3))
  {
    if (marker.second == cell_domain)
    {
      const dolfin::Cell c(*mesh, marker.second);
      const unsigned int* vertices = c.entities(0);
      editor.add_cell(cell_counter,
                      collected_vertices[vertices[0]],
                      collected_vertices[vertices[1]],
                      collected_vertices[vertices[2]],
                      collected_vertices[vertices[3]]);

    }
  }

  editor.close();
  return outmesh;
}


//-----------------------------------------------------------------------------
std::shared_ptr<dolfin::Mesh>
DolfinMeshUtils::merge_meshes(std::shared_ptr<dolfin::Mesh> m1,
                              std::shared_ptr<dolfin::Mesh> m2,
                              int m1_marker,
                              int m2_marker,
                              int m1_boundary_marker,
                              int m2_boundary_marker,
                              int interface_marker)
{
  log(dolfin::TRACE, "Merge meshes");

  if (m1->topology().dim() != m2->topology().dim() || m1->geometry().dim() != m2->geometry().dim())
    dolfin::dolfin_error("DolfinMeshUtils.cpp",
			 "merging meshes",
			 "Meshes dimensions must match");

  const std::size_t tdim = m1->topology().dim();

  // Map vertex (geometric) points to index for m1
  FuzzyPointMap vertex_map_inner(1e-14);
  std::vector<std::size_t> map_to_mesh_index(m1->num_vertices());;

  for (dolfin::VertexIterator v(*m1); !v.end(); ++v)
  {
    dolfin_assert(!vertex_map_inner.contains(v->point()));
    const std::size_t map_index = vertex_map_inner.insert_point(v->point());
    map_to_mesh_index[map_index] = v->global_index();
  }

  dolfin_assert(m1->num_vertices() == vertex_map_inner.size());

  // Count number of shared vertices
  // (We need this to initialize the mesh editor)
  // TOOO: Store these to avoid another lookup?
  std::size_t shared_vertices = 0;
  for (dolfin::VertexIterator v(*m2); !v.end(); ++v)
  {
    if (vertex_map_inner.contains(v->point()))
      shared_vertices++;
  }

  log(dolfin::TRACE, "Located %u shared vertices.", shared_vertices);

  std::shared_ptr<dolfin::Mesh> outmesh(new dolfin::Mesh);
  dolfin::MeshEditor editor;
  editor.open(*outmesh, m1->type().cell_type(), tdim, tdim);
  editor.init_vertices(m1->num_vertices() + m2->num_vertices() - shared_vertices);

  // Add vertices from m1
  for (dolfin::VertexIterator v(*m1); !v.end(); ++v)
  {
    editor.add_vertex(v->global_index(), v->point());
  }

  // Add vertices from m2
  std::vector<std::size_t> m2_vertex_map(m2->num_vertices());

  std::size_t m2_vertex_counter = m1->num_vertices();
  for (dolfin::VertexIterator v(*m2); !v.end(); ++v)
  {
    std::size_t index;
    if (vertex_map_inner.contains(v->point()))
    {
      const std::size_t map_index = vertex_map_inner.get_index(v->point());
      index = map_to_mesh_index[map_index];
    }
    else
    {
      index = m2_vertex_counter;
      editor.add_vertex(index, v->point());

      m2_vertex_counter++;
    }

    m2_vertex_map[v->global_index()] = index;
  }

  editor.init_cells(m1->num_cells() + m2->num_cells());

  std::size_t cell_counter = 0;

  std::vector<std::size_t> vertices(tdim+1);

  // Add cells from m1
  for (dolfin::CellIterator c(*m1); !c.end(); ++c)
  {
    for (std::size_t i = 0; i < tdim+1; i++)
      vertices[i] = c->entities(0)[i];

    editor.add_cell(cell_counter,vertices);
    cell_counter++;
  }

  for (dolfin::CellIterator c(*m2); !c.end(); ++c)
  {
    for (std::size_t i = 0; i < tdim+1; i++)
      vertices[i] = m2_vertex_map[c->entities(0)[i]];

    editor.add_cell(cell_counter,vertices);
    cell_counter++;
  }
  editor.close();

  return outmesh;

  // Mark cells from
  dolfin::MeshDomains &domain_markers = outmesh->domains();
  domain_markers.init(tdim-1);
  domain_markers.init(tdim);

  for (std::size_t i = 0; i < outmesh->num_cells(); ++i)
  {
    if (m1->num_cells())
    {
      if (m1_marker >= 0)
        domain_markers.set_marker(std::make_pair(i, m1_marker), tdim);
    }
    else
    {
      if (m2_marker >= 0)
        domain_markers.set_marker(std::make_pair(i, m2_marker), tdim);
    }
  }

  // Mark facets
  // Boundary facets are marked with 1 or 2 depending on the source meshes the
  // origined from. Facets shared by the two source meshes will be marked with
  // 3.

  // Compute facet - cell connectivity

  outmesh->init(tdim-1, tdim);
  for (dolfin::FacetIterator f(*outmesh); !f.end(); ++f)
  {
    if (f->exterior())
    {
      if (f->entities(tdim)[0] < m1->num_cells())
      {
        if (m1_boundary_marker >= 0)
          domain_markers.set_marker(std::make_pair(f->index(), m1_boundary_marker), tdim-1);
      }
      else
      {
        if (m2_boundary_marker >= 0)
          domain_markers.set_marker(std::make_pair(f->index(), m2_boundary_marker), tdim-1);
      }
    }
    else
    {
      if (interface_marker >= 0)
      {
        const std::size_t cell1 = f->entities(tdim)[0];
        const std::size_t cell2 = f->entities(tdim)[1];
        if ( (cell1 < m1->num_cells() && cell2 >= m1->num_cells()) ||
             (cell1 >= m1->num_cells() && cell2 < m1->num_cells()))
          domain_markers.set_marker(std::make_pair(f->index(), interface_marker), tdim-1);
      }
    }
  }

  return outmesh;
}

}
