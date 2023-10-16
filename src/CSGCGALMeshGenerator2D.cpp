// Copyright (C) 2012 Johannes Ring, 2012-2017 Benjamin Kehlet
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


#include <vector>
#include <cmath>
#include <limits>
#include <fstream>

#ifndef CGAL_HEADER_ONLY
#include <CGAL/compiler_config.h>
#endif
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>

#include <dolfin/common/constants.h>
#include <dolfin/common/MPI.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshDomains.h>
#include <dolfin/mesh/MeshValueCollection.h>
#include <dolfin/mesh/MeshPartitioning.h>
#include <dolfin/log/log.h>

#include <mshr/CSGCGALMeshGenerator2D.h>
#include <mshr/CSGGeometry.h>
#include <mshr/CSGOperators.h>
#include <mshr/CSGPrimitives2D.h>
#include <mshr/CSGCGALDomain2D.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel             Inexact_Kernel;
typedef CGAL::Delaunay_mesh_vertex_base_2<Inexact_Kernel>               Vertex_base;
typedef CGAL::Delaunay_mesh_face_base_2<Inexact_Kernel>                 Face_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base>    TDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<Inexact_Kernel, TDS> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>                        Mesh_criteria_2;
typedef CGAL::Delaunay_mesher_2<CDT, Mesh_criteria_2>                   CGAL_Mesher_2;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Face_handle Face_handle;
typedef Inexact_Kernel::Point_2 Point_2;

namespace mshr
{

//-----------------------------------------------------------------------------
CSGCGALMeshGenerator2D::CSGCGALMeshGenerator2D()
{
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
CSGCGALMeshGenerator2D::~CSGCGALMeshGenerator2D() {}
//-----------------------------------------------------------------------------
void explore_subdomain(const CDT& cdt,
                       Face_handle f,
                       std::size_t marker,
                       std::map<Face_handle, unsigned int>& subdomain_map,
                       std::set<Face_handle>& not_visited)
{
  std::deque<Face_handle> queue;
  queue.push_back(f);

  while (!queue.empty())
  {
    CDT::Face_handle face = queue.front();
    queue.pop_front();

    std::set<Face_handle>::iterator face_not_visited = not_visited.find(face);
    if (face_not_visited != not_visited.end())
    {
      dolfin_assert(subdomain_map.find(face) == subdomain_map.end());
      not_visited.erase(face_not_visited);
      subdomain_map[face] = marker;

      for(int i = 0; i < 3; i++)
      {
        const CDT::Edge e(face,i);
        if (!cdt.is_constrained(e))
          queue.push_back(face->neighbor(i));
      }
    }
  }
}
//-----------------------------------------------------------------------------
// Set the member in_domain and counter for all faces in the cdt
std::map<Face_handle, unsigned int>
  explore_subdomains(const CDT& cdt,
                     const CSGCGALDomain2D& total_domain,
                     const std::vector<std::pair<std::size_t, CSGCGALDomain2D>>&
                       subdomain_geometries)
{
  std::set<Face_handle> not_visited;

  for(CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
        fit != cdt.finite_faces_end(); ++fit)
  {
    if (fit->is_in_domain())
      not_visited.insert(fit);
  }

  std::map<Face_handle, unsigned int> subdomain_map;

  std::list<CDT::Face_handle> queue;
  queue.push_back(*not_visited.begin());

  while (!not_visited.empty())
  {
    const Face_handle f = *not_visited.begin();

    // set to default domain (0)
    unsigned int marker = 0;

    const Point_2 p0 = f->vertex(0)->point();
    const Point_2 p1 = f->vertex(1)->point();
    const Point_2 p2 = f->vertex(2)->point();

    const dolfin::Point p( (p0[0] + p1[0] + p2[0]) / 3.0,
                           (p0[1] + p1[1] + p2[1]) / 3.0 );

    for (const std::pair<std::size_t, CSGCGALDomain2D> subdomain : subdomain_geometries)
    {
      if (subdomain.second.point_in_domain(p))
      {
        marker = subdomain.first;
        break;
      }
    }

    explore_subdomain(cdt,
                      f,
                      marker,
                      subdomain_map,
                      not_visited);
  }

  return subdomain_map;
}
//-----------------------------------------------------------------------------
double shortest_constrained_edge(const CDT &cdt)
{
  double min_length = std::numeric_limits<double>::max();
  for (CDT::Finite_edges_iterator it = cdt.finite_edges_begin();
       it != cdt.finite_edges_end();
       it++)
  {
    if (!cdt.is_constrained(*it))
      continue;

    // An edge is an std::pair<Face_handle, int>
    // see CGAL/Triangulation_data_structure_2.h
    CDT::Face_handle f = it->first;
    const int i = it->second;

    CDT::Point p1 = f->vertex( (i+1)%3 )->point();
    CDT::Point p2 = f->vertex( (i+2)%3 )->point();

    min_length = std::min(CGAL::to_double((p1-p2).squared_length()),
                          min_length);
  }

  return min_length;
}
//-----------------------------------------------------------------------------
std::shared_ptr<dolfin::Mesh> CSGCGALMeshGenerator2D::generate(const std::shared_ptr<const CSGCGALDomain2D> total_domain,
                                                               const std::vector<std::pair<std::size_t,
                                                                                           std::shared_ptr<const CSGCGALDomain2D>>>& subdomains)
{
  const bool partition = parameters["partition"];

  std::shared_ptr<dolfin::Mesh> mesh(new dolfin::Mesh);

  // Note that if not in parallel (ie. size() == 0)
  // then both receiver and broadcaster will return false
  // If not partition, then generate the mesh on all processes
  if (!partition || !dolfin::MPI::is_receiver(mesh->mpi_comm()))
  {
    double cell_size;
    {
      const double mesh_resolution = parameters["mesh_resolution"];
      if (mesh_resolution > 0)
      {
        const double min_radius = total_domain->compute_boundingcircle_radius();
        cell_size = 2.0*min_radius/mesh_resolution;
      }
      else
        cell_size = parameters["cell_size"];
    }

    std::vector<std::pair<std::size_t, CSGCGALDomain2D> >
      subdomain_geometries;

    log(dolfin::TRACE, "Converting geometry to CGAL polygon");

    // Empty polygon, will be populated when traversing the subdomains
    CSGCGALDomain2D overlaying;

    if (!subdomains.empty())
      log(dolfin::TRACE, "Processing subdomains");

    // Add the subdomains to the PSLG. Traverse in reverse order to get the latest
    // added subdomain on top
    for (auto rit = subdomains.rbegin(); rit != subdomains.rend(); rit++)
    {
      const std::size_t current_index = rit->first;
      CSGCGALDomain2D cgal_geometry(*(rit->second));

      // Only the part inside the total domain
      cgal_geometry.intersect_inplace(*total_domain, 1e-15);

      // Only the part outside overlaying subdomains
      cgal_geometry.difference_inplace(overlaying);

      subdomain_geometries.push_back(std::make_pair(current_index,
                                                    cgal_geometry));

      overlaying.join_inplace(cgal_geometry);
    }

    CSGCGALDomain2D remaining(*total_domain);
    remaining.difference_inplace(overlaying);

    subdomain_geometries.push_back(std::make_pair(0, remaining));

    const auto pslg = CSGCGALDomain2D::compute_pslg(subdomain_geometries);


    // for (dolfin::Point p : pslg.first)
    // {
    //   std::cout << "\"Point " << p.x() << " " << p.y() << "\" ";
    // }
    // std::cout << std::endl;

    // for (std::pair<std::size_t, std::size_t> s : pslg.second)
    // {
    //   dolfin::Point p0 = pslg.first[s.first];
    //   dolfin::Point p1 = pslg.first[s.second];
    //   std::cout << "\"Segment " << p0.x() << " " << p0.y() << ", " << p1.x() << " " << p1.y() << "\" ";
    // }
    // std::cout << std::endl;

    // Create empty CGAL triangulation and copy data from the PSLG
    CDT cdt;

    {
      std::vector<CDT::Vertex_handle> vertices;
      {
        // Insert the vertices into the triangulation
        vertices.reserve(pslg.first.size());
        for (const dolfin::Point& vertex : pslg.first)
        {
          const Point_2 p(vertex.x(), vertex.y());
          vertices.push_back(cdt.insert(p));
        }
      }

      // Insert the edges as constraints
      for (const std::pair<std::size_t, size_t>& edge : pslg.second)
      {
        cdt.insert_constraint(vertices[edge.first], vertices[edge.second]);
      }
    }

    log(dolfin::TRACE, "Initializing mesh refinement");

    // Create mesher
    CGAL_Mesher_2 mesher(cdt);

    // Add seeds for all faces in the total domain
    std::vector<dolfin::Point> hole_points;
    total_domain->get_points_in_holes(hole_points);
    std::vector<Point_2> seed_points;
    for (dolfin::Point p : hole_points)
      seed_points.push_back(Point_2(p.x(), p.y()));

    log(dolfin::TRACE, "Added %d seed points", seed_points.size());
    mesher.set_seeds(seed_points.begin(), seed_points.end(), false);

    // Set shape and size criteria
    const double shape_bound = parameters["triangle_shape_bound"];
    mesher.set_criteria(Mesh_criteria_2(shape_bound, cell_size));

    // Refine CGAL mesh/triangulation
    mesher.refine_mesh();

    // Lloyd optimization
    if (parameters["lloyd_optimize"]) {
      log(dolfin::TRACE, "Optimizing mesh by Lloyd smoothing");
      CGAL::lloyd_optimize_mesh_2(cdt);
    }

    // Make sure triangulation is valid
    dolfin_assert(cdt.is_valid());

    // Mark the subdomains
    log(dolfin::TRACE, "Exploring subdomains in mesh");

    std::map<Face_handle, unsigned int> subdomain_map =
      explore_subdomains(cdt, *total_domain, subdomain_geometries);

    // Clear mesh
    // This function has been removed from dolfin.
    // mesh.clear();

    const std::size_t num_vertices = cdt.number_of_vertices();

    // Count valid cells
    std::size_t num_cells = 0;
    for (CDT::Finite_faces_iterator cgal_cell = cdt.finite_faces_begin();
         cgal_cell != cdt.finite_faces_end(); ++cgal_cell)
    {
      // Add cell if it is in the domain
      if (cgal_cell->is_in_domain())
      {
        num_cells++;
      }
    }

    log(dolfin::DBG, "Mesh with %d vertices and %d cells created", num_vertices, num_cells);

    // Create a MeshEditor and open
    dolfin::MeshEditor mesh_editor;
    mesh_editor.open(*mesh, dolfin::CellType::Type::triangle, 2, 2);
    mesh_editor.init_vertices(num_vertices);
    mesh_editor.init_cells(num_cells);

    // Add vertices to mesh
    unsigned int vertex_index = 0;
    std::map<Vertex_handle, unsigned int> vertex_map;
    for (CDT::Finite_vertices_iterator cgal_vertex = cdt.finite_vertices_begin();
         cgal_vertex != cdt.finite_vertices_end(); ++cgal_vertex)
    {
      // Get vertex coordinates and add vertex to the mesh
      dolfin::Point p;
      p[0] = cgal_vertex->point()[0];
      p[1] = cgal_vertex->point()[1];

      // Add mesh vertex
      mesh_editor.add_vertex(vertex_index, p);

      // Attach index to vertex and increment
      vertex_map[cgal_vertex] = vertex_index;
      vertex_index++;
    }

    dolfin_assert(vertex_index == num_vertices);

    // Add cells to mesh and build domain marker mesh function
    dolfin::MeshDomains &domain_markers = mesh->domains();
    std::size_t cell_index = 0;
    const bool mark_cells = !subdomains.empty();
    for (CDT::Finite_faces_iterator cgal_cell = cdt.finite_faces_begin();
         cgal_cell != cdt.finite_faces_end(); ++cgal_cell)
    {
      // Add cell if it is in the domain
      if (cgal_cell->is_in_domain())
      {
        mesh_editor.add_cell(cell_index,
                             vertex_map[cgal_cell->vertex(0)],
                             vertex_map[cgal_cell->vertex(1)],
                             vertex_map[cgal_cell->vertex(2)]);

        if (mark_cells)
        {
          domain_markers.set_marker(std::make_pair(cell_index,
                                                   subdomain_map[cgal_cell]), 2);
        }
        ++cell_index;
      }
    }
    dolfin_assert(cell_index == num_cells);

    // Close mesh editor
    mesh_editor.close();
  }

  // Distribute the mesh (if in parallel)
  if (partition)
    dolfin::MeshPartitioning::build_distributed_mesh(*mesh);

  return mesh;
}
}
//-----------------------------------------------------------------------------
