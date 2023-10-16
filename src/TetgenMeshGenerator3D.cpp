// Copyright (C) 2014 Benjamin Kehlet
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

#include <mshr/TetgenMeshGenerator3D.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshPartitioning.h>
#include <dolfin/log/log.h>

#include <tetgen.h>

// Bounding sphere computation
#include <CGAL/Cartesian.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>

#include <vector>
#include <array>

namespace
{
//-----------------------------------------------------------------------------
void build_dolfin_mesh(const tetgenio& tetgenmesh, dolfin::Mesh& dolfinmesh)
{
  // Clear mesh
  // dolfin::Mesh::clear has been removed
  // dolfinmesh.clear();

  // Create and initialize mesh editor
  dolfin::MeshEditor mesh_editor;
  mesh_editor.open(dolfinmesh, dolfin::CellType::Type::tetrahedron, 3, 3);
  mesh_editor.init_vertices(tetgenmesh.numberofpoints);

  const int offset = tetgenmesh.firstnumber;

  for (int i = 0; i < tetgenmesh.numberofpoints; i++)
  {

    dolfin::Point p(tetgenmesh.pointlist[i * 3],
                    tetgenmesh.pointlist[i * 3 + 1],
                    tetgenmesh.pointlist[i * 3 + 2]);
    mesh_editor.add_vertex(i, p);
  }


  mesh_editor.init_cells(tetgenmesh.numberoftetrahedra);
  dolfin_assert(tetgenmesh.numberofcorners == 4);

  for (int i = 0; i < tetgenmesh.numberoftetrahedra; i++)
  {
    mesh_editor.add_cell(i,
                         tetgenmesh.tetrahedronlist[i*4 + 0]-offset,
                         tetgenmesh.tetrahedronlist[i*4 + 1]-offset,
                         tetgenmesh.tetrahedronlist[i*4 + 2]-offset,
                         tetgenmesh.tetrahedronlist[i*4 + 3]-offset);
  }

  // Close mesh editor
  mesh_editor.close();
}
//-----------------------------------------------------------------------------
double bounding_sphere_radius(const std::vector<double>& vertices)
{
  typedef double FT;
  typedef CGAL::Cartesian<FT> K;
  typedef CGAL::Min_sphere_of_spheres_d_traits_3<K, FT> MinSphereTraits;
  typedef CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
  typedef MinSphereTraits::Sphere Sphere;

  std::vector<Sphere> S;

  for (std::size_t i = 0; i < vertices.size(); i += 3)
  {
    S.push_back(Sphere(K::Point_3( vertices[i],
                                   vertices[i+1],
                                   vertices[i+2]), 0.0));
  }

  Min_sphere ms(S.begin(), S.end());
  dolfin_assert(ms.is_valid());

  return ms.radius();
}

} // end anonymous namespace
//-----------------------------------------------------------------------------
namespace mshr
{

TetgenMeshGenerator3D::TetgenMeshGenerator3D()
{
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
TetgenMeshGenerator3D::~TetgenMeshGenerator3D()
{
}
//-----------------------------------------------------------------------------
std::shared_ptr<dolfin::Mesh>
TetgenMeshGenerator3D::generate(std::shared_ptr<const CSGCGALDomain3D> domain) const
{
  tetgenio in;

  double r;
  {
    // Copy the vertices to the tetgen structure
    dolfin::log(dolfin::TRACE, "Copying vertices");
    std::unique_ptr<const std::vector<double>> vertices = domain->get_vertices();
    r = bounding_sphere_radius(*vertices);

    in.numberofpoints = vertices->size()/3;
    in.pointlist = new REAL[in.numberofpoints * 3];

    for (std::size_t i = 0; i < vertices->size(); i++)
    {
      in.pointlist[i] = (*vertices)[i];
    }
  }

  // Copy the facets to the tetgen structure
  {
    dolfin::log(dolfin::TRACE, "Copying facets");
    std::unique_ptr<const std::vector<std::size_t>> facets = domain->get_facets();

    in.numberoffacets = facets->size()/3;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    //in.facetmarkerlist = new int[in.numberoffacets];

    for (std::size_t i = 0; i*3 < facets->size(); i++)
    {
      tetgenio::facet& f = in.facetlist[i];
      f.numberofpolygons = 1;
      f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
      f.numberofholes = 0;
      f.holelist = NULL;
      tetgenio::polygon& p = f.polygonlist[0];
      p.numberofvertices = 3;
      p.vertexlist = new int[p.numberofvertices];
      p.vertexlist[0] = (*facets)[i*3];
      p.vertexlist[1] = (*facets)[i*3+1];
      p.vertexlist[2] = (*facets)[i*3+2];
    }
  }

  // Mark holes in the domain
  {
    dolfin::log(dolfin::TRACE, "Marking holes");
    std::vector<dolfin::Point> holes;
    domain->get_points_in_holes(holes);

    in.numberofholes = holes.size();
    in.holelist = new REAL[in.numberofholes*3];
    std::size_t i = 0;
    for (std::vector<dolfin::Point>::const_iterator it = holes.begin();
         it != holes.end(); it++)
    {
      in.holelist[i*3]     = it->x();
      in.holelist[i*3 + 1] = it->y();
      in.holelist[i*3 + 2] = it->z();
    }

    i++;
  }

  // Release domain object, possibly deleting it
  domain.reset();

  // set tetgen parameters
  std::stringstream tetgenparams;
  tetgenparams << std::fixed << std::setprecision(16);

  // tetrahedralize a plc
  tetgenparams << "p";

  if (!parameters["disable_quality_improvement"])
  {
    // set quality constraints
    const double ratio = parameters["max_radius_edge_ratio"];
    const double angle = parameters["min_dihedral_angle"];
    tetgenparams << "q" << ratio << "/" << angle;

    tetgenparams << "a";
    if (double(parameters["max_tet_volume"]) > 0)
    {
      // set maximum cell volume
      tetgenparams << double(parameters["max_tet_volume"]);
    }
    else
    {
      const double resolution = parameters["mesh_resolution"];

      // try to compute reasonable parameters
      //const double r = bounding_sphere_radius(vertices);
      const double cell_size = r/static_cast<double>(resolution)*2.0;
      tetgenparams << cell_size;
    }

    if (parameters["preserve_surface"])
    {
      tetgenparams << "Y";
    }
  }

  if (dolfin::get_log_level() > dolfin::DBG)
  {
    // set verbosity level
    tetgenparams << "Q";
  }

  dolfin::log(dolfin::TRACE, "Calling tetgen with parameters: " + tetgenparams.str());

  // Tetgen requires a char[] (as opposed to a const char[])
  // so we need to copy of from the string
  const std::string str = tetgenparams.str();
  std::unique_ptr<char> writable(new char[str.size() + 1]);
  std::copy(str.begin(), str.end(), writable.get());
  writable.get()[str.size()] = '\0'; // terminating 0

  tetgenio out;
  tetrahedralize(writable.get(), &in, &out);

  std::shared_ptr<dolfin::Mesh> mesh(new dolfin::Mesh());
  build_dolfin_mesh(out, *mesh);

  // Distribute the mesh (if in parallel)
  dolfin::MeshPartitioning::build_distributed_mesh(*mesh);

  return mesh;
}
//-----------------------------------------------------------------------------
} // end namespace mshr
