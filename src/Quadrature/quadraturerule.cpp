#include "quadraturerule.hpp"

#include <quad2d.hpp>
#include <quad1d.hpp>
#include <vertex.hpp>

namespace HArDCore2D
{

  QuadratureRule generate_quadrature_rule(
					  const Cell & T,
					  const int doe,
					  const bool force_split
					  )
  { 
    QuadRuleTriangle quadCell(std::max(doe,0), true);
    QuadratureRule quad;
    size_t nedges = T.n_edges();

    if ( (nedges == 3) && (!force_split) ) {
        // Triangle
        auto x0 = T.vertex(0)->coords();
        auto x1 = T.vertex(1)->coords();
        auto x2 = T.vertex(2)->coords();

        double xT[] = {x0.x(), x1.x(), x2.x()};
        double yT[] = {x0.y(), x1.y(), x2.y()};

        quadCell.setup(xT, yT);
        for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
            quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.wq(iqn));
        }

    } else if ( (nedges == 4) && (!force_split) ) {  
        // quadrilateral split into two triangles
        auto x0 = T.vertex(0)->coords();
        for (size_t isplit = 0; isplit < 2; isplit++) {
            auto x1 = T.vertex(1 + isplit)->coords();
            auto x2 = T.vertex(2 + isplit)->coords();
            double xTr[] = {x0.x(), x1.x(), x2.x()};
            double yTr[] = {x0.y(), x1.y(), x2.y()};
            quadCell.setup(xTr, yTr);

            for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
                quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.wq(iqn));
            }
        }


    } else {
        // split at barycentre into triangles
        auto xT = T.center_mass();

        for (size_t ilF = 0; ilF < nedges; ilF++) {
            const Edge* e = T.edge(ilF);
            auto x0 = e->vertex(0)->coords();
            auto x1 = e->vertex(1)->coords();

            double xTr[] = {xT.x(), x0.x(), x1.x()};
            double yTr[] = {xT.y(), x0.y(), x1.y()};
            quadCell.setup(xTr, yTr);

            for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
                quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.wq(iqn));
            }
        }
    }

    return quad;

  }

  //------------------------------------------------------------------------------

  QuadratureRule generate_quadrature_rule(
					  const Edge & E,
					  const int doe
					  )
  {
    QuadRuleEdge quadEdge(std::max(doe,0), true);
    QuadratureRule quad;

    auto x0 = E.vertex(0)->coords();
    auto x1 = E.vertex(1)->coords();
    double xT[] = {x0.x(), x1.x()};
    double yT[] = {x0.y(), x1.y()};
  
    quadEdge.setup(xT, yT);
    for (size_t iqn = 0; iqn < quadEdge.nq(); iqn++) {
        quad.emplace_back(quadEdge.xq(iqn), quadEdge.yq(iqn), quadEdge.wq(iqn));
    }
    return quad;
};

} // End namespace HArDCore2D
