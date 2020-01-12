// Class to create and store values of cell and edge basis functions on quadrature points
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "elementquad.hpp"
#include <string>
#include <mesh.hpp>
//#include <vector>
//#include <cell.hpp>
//#include <edge.hpp>

using namespace HArDCore2D;

// Creation class

ElementQuad::ElementQuad(const HybridCore& hho, const size_t iT, const size_t doeT, const size_t doeF )
  : _hho(hho),
    _iT(iT),
    _doeT(doeT),
    _doeF(doeF) {
      // Create cell elements
      _quadT = _hho.cell_qrule(_iT, _doeT);
      _phiT_quadT = _hho.basis_quad("cell", _iT, _quadT, _hho.K()+1);
      _dphiT_quadT = _hho.grad_basis_quad(_iT, _quadT, _hho.K()+1);

      // Create edge elements
      const Mesh* mesh = _hho.get_mesh_ptr();
      Cell* cell = mesh->cell(_iT);
      size_t nedges = cell->n_edges();
      _quadF.resize(nedges);
      _phiT_quadF.resize(nedges);
      _phiF_quadF.resize(nedges);
      _dphiT_quadF.resize(nedges);
      for (size_t ilF = 0; ilF < nedges ; ilF++){
        const size_t iF = cell->edge(ilF)->global_index();
        _quadF[ilF] = _hho.edge_qrule(iF, _doeF);
        _phiT_quadF[ilF] = _hho.basis_quad("cell", _iT, _quadF[ilF], _hho.K()+1);
        _phiF_quadF[ilF] = _hho.basis_quad("edge", iF, _quadF[ilF], _hho.K());
        _dphiT_quadF[ilF] = _hho.grad_basis_quad(_iT, _quadF[ilF], _hho.K()+1);
      }

  }

//ElementQuad::~ElementQuad() {}

