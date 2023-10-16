#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <mshr/CSGGeometry.h>
#include <mshr/CSGPrimitives2D.h>
#include <mshr/CSGPrimitives3D.h>
#include <mshr/CSGCGALDomain2D.h>
#include <mshr/CSGCGALDomain3D.h>
#include <mshr/CSGOperators.h>
#include <mshr/MeshGenerator.h>
#include <mshr/CSGGeometries3D.h>
#include <mshr/CSGCGALMeshGenerator3D.h>
#include <mshr/CSGCGALMeshGenerator2D.h>
#include <mshr/TetgenMeshGenerator3D.h>
#include <mshr/Meshes.h>

namespace py = pybind11;


PYBIND11_MODULE(cpp, m)
{
  // Create module for C++ wrappers
  m.doc() ="mshr python interface";

  // CSGGeometry
  py::class_<mshr::CSGGeometry, dolfin::Variable,
             std::shared_ptr<mshr::CSGGeometry>> (m, "CSGGeometry")
    .def("dim", &mshr::CSGGeometry::dim)
    .def("set_subdomain",
         static_cast<void(mshr::CSGGeometry::*)(std::size_t i,
                                                std::shared_ptr<mshr::CSGGeometry> s)>(&mshr::CSGGeometry::set_subdomain))
    .def("has_subdomains", &mshr::CSGGeometry::has_subdomains)
    .def("inside", static_cast<bool (mshr::CSGGeometry::*)(dolfin::Point) const>(&mshr::CSGGeometry::inside))
    .def("__mul__",
         static_cast<std::shared_ptr<mshr::CSGIntersection>(*)(std::shared_ptr<mshr::CSGGeometry>,
                                                               std::shared_ptr<mshr::CSGGeometry>)>(&mshr::operator*),
         py::is_operator())

    .def("__mul__",
         static_cast<std::shared_ptr<mshr::CSGScaling>(*)(std::shared_ptr<mshr::CSGGeometry>,
                                                          double)>(&mshr::operator*),
         py::is_operator())

    .def("__add__",
         static_cast<std::shared_ptr<mshr::CSGUnion>(*)(std::shared_ptr<mshr::CSGGeometry>,
                                                        std::shared_ptr<mshr::CSGGeometry>)>(&mshr::operator+),
         py::is_operator())

    .def("__sub__",
         static_cast<std::shared_ptr<mshr::CSGDifference>(*)(std::shared_ptr<mshr::CSGGeometry>,
                                                             std::shared_ptr<mshr::CSGGeometry>)>(&mshr::operator-),
         py::is_operator());


  py::class_<mshr::CSGUnion, mshr::CSGGeometry,
             std::shared_ptr<mshr::CSGUnion>>(m, "CSGUnion");

  py::class_<mshr::CSGIntersection, mshr::CSGGeometry,
             std::shared_ptr<mshr::CSGIntersection>>(m, "CSGIntersection");

  py::class_<mshr::CSGDifference, mshr::CSGGeometry,
             std::shared_ptr<mshr::CSGDifference>>(m, "CSGDifference");

  py::class_<mshr::CSGScaling, mshr::CSGGeometry,
             std::shared_ptr<mshr::CSGScaling>>(m, "CSGScaling")
    .def(py::init<std::shared_ptr<mshr::CSGGeometry>, double>());

  py::class_<mshr::CSGTranslation, mshr::CSGGeometry,
             std::shared_ptr<mshr::CSGTranslation>>(m, "CSGTranslation")
    .def(py::init<std::shared_ptr<mshr::CSGGeometry>, dolfin::Point>());

  py::class_<mshr::CSGRotation, mshr::CSGGeometry,
             std::shared_ptr<mshr::CSGRotation>>(m, "CSGRotation")
    .def(py::init<std::shared_ptr<mshr::CSGGeometry>, dolfin::Point, double>())
   .def(py::init<std::shared_ptr<mshr::CSGGeometry>, double>());



  // Circle
  py::class_<mshr::Circle, mshr::CSGGeometry,
    std::shared_ptr<mshr::Circle>>(m, "Circle")
    .def(py::init<dolfin::Point, double, std::size_t>(),
         py::arg("c"), py::arg("r"), py::arg("segments")=0)
    .def("center", &mshr::Circle::center)
    .def("radius", &mshr::Circle::radius);

  // Ellipse
  py::class_<mshr::Ellipse, mshr::CSGGeometry,
             std::shared_ptr<mshr::Ellipse>>(m, "Ellipse")
    .def(py::init<dolfin::Point, double, double, std::size_t>(),
         py::arg("c"), py::arg("a"), py::arg("b"), py::arg("segments")=0)
    .def("center", &mshr::Ellipse::center)
    .def("a", &mshr::Ellipse::a)
    .def("b", &mshr::Ellipse::b);

  // Rectangle
  py::class_<mshr::Rectangle, mshr::CSGGeometry,
             std::shared_ptr<mshr::Rectangle>> (m, "Rectangle")
    .def(py::init<dolfin::Point, dolfin::Point>())
    .def("first_corner", &mshr::Rectangle::first_corner)
    .def("second_corne", &mshr::Rectangle::second_corner);

  // Polygon
  py::class_<mshr::Polygon, mshr::CSGGeometry,
             std::shared_ptr<mshr::Polygon>>(m, "Polygon")
    .def(py::init<std::vector<dolfin::Point>&>())
    .def("ccw", &mshr::Polygon::ccw)
    .def("vertices", &mshr::Polygon::vertices);

  // Sphere
  py::class_<mshr::Sphere, mshr::CSGGeometry,
             std::shared_ptr<mshr::Sphere>>(m, "Sphere")
    .def(py::init<dolfin::Point, double, std::size_t>(),
         py::arg("center"), py::arg("radius"), py::arg("segments")=10);

  // Box
  py::class_<mshr::Box, mshr::CSGGeometry,
             std::shared_ptr<mshr::Box>>(m, "Box")
    .def(py::init<dolfin::Point, dolfin::Point>());

  // Cylinder
  py::class_<mshr::Cylinder, mshr::CSGGeometry,
             std::shared_ptr<mshr::Cylinder>>(m, "Cylinder")
    .def(py::init<dolfin::Point, dolfin::Point, double, double, std::size_t>(),
         py::arg("top"),
         py::arg("bottom"),
         py::arg("top_radius"),
         py::arg("bottom_radius"),
         py::arg("segments")=32);

  // Cone
  py::class_<mshr::Cone, mshr::Cylinder,
             std::shared_ptr<mshr::Cone>>(m, "Cone")
    .def(py::init<dolfin::Point, dolfin::Point, double, std::size_t>(),
         py::arg("top"), py::arg("bottom"), py::arg("r"), py::arg("segments")=32);

  // Tetrahedron
  py::class_<mshr::Tetrahedron, mshr::CSGGeometry,
             std::shared_ptr<mshr::Tetrahedron>>(m, "Tetrahedron")
    .def(py::init<dolfin::Point, dolfin::Point, dolfin::Point, dolfin::Point>());

  // Surface3D (deprecated, use CSGCGALDomain3D)
  py::class_<mshr::Surface3D, mshr::CSGGeometry,
             std::shared_ptr<mshr::Surface3D>>(m, "Surface3D")
    .def(py::init<std::string>());

  // Ellipsoid
  py::class_<mshr::Ellipsoid, mshr::CSGGeometry,
             std::shared_ptr<mshr::Ellipsoid>>(m, "Ellipsoid")
    .def(py::init<dolfin::Point, double, double, double, std::size_t>(),
         py::arg("center"),
         py::arg("a"), py::arg("b"), py::arg("c"), py::arg("segments")=15);

  // Extrude2D
  py::class_<mshr::Extrude2D, mshr::CSGGeometry,
             std::shared_ptr<mshr::Extrude2D>>(m, "Extrude2D")
    .def(py::init<std::shared_ptr<mshr::CSGGeometry>, double>());

  // CSGCGALDomain3D
  py::class_<mshr::CSGCGALDomain3D, mshr::CSGGeometry,
             std::shared_ptr<mshr::CSGCGALDomain3D>>(m, "CSGCGALDomain3D")
    .def(py::init<std::shared_ptr<const mshr::CSGGeometry>>())
    .def("num_vertices", &mshr::CSGCGALDomain3D::num_vertices)
    .def("num_facets", &mshr::CSGCGALDomain3D::num_facets)
    .def("num_halfedges", &mshr::CSGCGALDomain3D::num_halfedges)
    .def("num_degenerate_facets", &mshr::CSGCGALDomain3D::num_degenerate_facets)
    .def("facet_area_minmax", &mshr::CSGCGALDomain3D::facet_area_minmax)
    .def("edge_length_minmax", &mshr::CSGCGALDomain3D::edge_length_range)
    .def("num_short_edges", &mshr::CSGCGALDomain3D::num_short_edges)
    .def("volume", &mshr::CSGCGALDomain3D::volume)
    .def("ensure_meshing_preconditions", &mshr::CSGCGALDomain3D::ensure_meshing_preconditions)
    .def("is_selfintersecting", &mshr::CSGCGALDomain3D::is_selfintersecting)
    .def("remove_selfintersectione", &mshr::CSGCGALDomain3D::remove_selfintersections)
    .def("save", &mshr::CSGCGALDomain3D::save)
    .def("num_disconnected_components", &mshr::CSGCGALDomain3D::num_disconnected_components)
    .def("remesh_surface", &mshr::CSGCGALDomain3D::remesh_surface)
    .def("remove_degenerate_facets", &mshr::CSGCGALDomain3D::remove_degenerate_facets)
    .def("convex_hull", static_cast<std::shared_ptr<mshr::CSGCGALDomain3D>(mshr::CSGCGALDomain3D::*)() const>(&mshr::CSGCGALDomain3D::convex_hull))
    ;

  py::class_<mshr::CSGGeometries>(m, "CSGGeometries")
    .def_static("lego", &mshr::CSGGeometries::lego)
    .def_static("propeller", &mshr::CSGGeometries::propeller);

  py::class_<mshr::UnitSphereMesh>(m, "UnitSphereMesh")
    .def(py::init<std::size_t>());

  py::class_<mshr::CSGCGALMeshGenerator3D, dolfin::Variable,
	     std::shared_ptr<mshr::CSGCGALMeshGenerator3D>>(m, "CSGCGALMeshGenerator3D")
    .def(py::init<>())
    .def("generate", static_cast<std::shared_ptr<dolfin::Mesh>(mshr::CSGCGALMeshGenerator3D::*)(std::shared_ptr<const mshr::CSGCGALDomain3D>) const>(&mshr::CSGCGALMeshGenerator3D::generate));

  py::class_<mshr::CSGCGALDomain2D, dolfin::Variable,
	     std::shared_ptr<mshr::CSGCGALDomain2D>>(m, "CSGCGALDomain2D")
    .def(py::init<std::shared_ptr<const mshr::CSGGeometry>, double>());

  py::class_<mshr::CSGCGALMeshGenerator2D, dolfin::Variable,
	     std::shared_ptr<mshr::CSGCGALMeshGenerator2D>>(m, "CSGCGALMeshGenerator2D")
    .def(py::init<>())
    .def("generate", &mshr::CSGCGALMeshGenerator2D::generate);


  py::class_<mshr::TetgenMeshGenerator3D, dolfin::Variable,
	     std::shared_ptr<mshr::TetgenMeshGenerator3D>>(m, "TetgenMeshGenerator3D")
    .def(py::init<>())
    .def("generate", &mshr::TetgenMeshGenerator3D::generate);

  // generate_mesh
  m.def("generate_mesh",
	&mshr::generate_mesh,
	py::arg("geoemtry"), py::arg("resolution"), py::arg("backend")="cgal");
}
