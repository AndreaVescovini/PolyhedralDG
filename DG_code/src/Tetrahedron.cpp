#include "Tetrahedron.hpp"
#include <Eigen/Core>
#include <cmath>

namespace geom {

Tetrahedron::Tetrahedron(Vertex& v1,  Vertex& v2, Vertex& v3, Vertex& v4,
                         Polyhedron* poly)
  : id_{counter_}, vertices_{{v1, v2, v3, v4}}, poly_{poly}
{
  counter_++;

  map_ = (Eigen::Matrix4d() << (v2.getCoords() - v1.getCoords()),
                               (v3.getCoords() - v1.getCoords()),
                               (v4.getCoords() - v1.getCoords()),
                                v1.getCoords(), Eigen::RowVector4d::Zero()).finished();
  // std::cout << map_.linear().matrix() << std::endl;
  absDetJacobian_ = std::abs(map_.linear().matrix().determinant());
}

Tetrahedron::Tetrahedron(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4,
                         Polyhedron& poly)
  : Tetrahedron(v1, v2, v3, v4, &poly) {}

void Tetrahedron::resetCounter(unsigned counter)
{
  counter_ = counter;
}

std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra)
{
  out << tetra.id_ << " " << "V: " << tetra.vertices_[0].get().getId() << " "
                                   << tetra.vertices_[1].get().getId() << " "
                                   << tetra.vertices_[2].get().getId() << " "
                                   << tetra.vertices_[3].get().getId() << ", P: "
                                   << tetra.poly_->getId();
  return out;
}

unsigned Tetrahedron::counter_ = 0;

} // namespace geom
