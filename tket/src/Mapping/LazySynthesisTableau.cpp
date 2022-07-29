// Copyright 2019-2022 Cambridge Quantum Computing
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "Mapping/LazySynthesisAAS.hpp"

namespace tket {


void LazySynthesisTableau::add_qubit(const Qubit& qubit, const Edge& edge) 
    // merge tableaux
    // if there are overlapping qubits this method will throw an error, so don't check first
    this->tableau_ = UnitaryTableau::add_disjoint_tableaux(this->tableau_, UnitaryTableau({qb}));
    // add new entry to dependencies
    this->dependencies_.insert({qb, {}});
    // add new entry to qubit edge map
    this->q_in_hole.insert({qubit, edge});
    // no change to held pauli gadgets or vertices
}

void LazySynthesisTableau::merge_tableaux(const LazySynthesisTableau& merge){
  // merge tableaux
  // if there are overlapping qubits this method will throw an error, so don't check first
  this->tableau_ = UnitaryTableau::add_disjoint_tableaux(this->tableau_, merge.tableau_);
  // merge dependencies maps
  this->dependencies_.insert(
      merge.dependencies_.begin(), merge.dependencies_.end());
  // merge qubit edge map
  this->q_in_hole.insert(merge.q_in_hole.begin(), merge.q_in_hole.end());
  // merge verts
  this->verts.insert(merge.verts.begin(), merge.verts.end());
  // merge pauli gadgets
  this->pauli_gadgets_.insert(this->pauli_gadgets_.end(), merge.pauli_gadgets_.begin(), merge.pauli_gadgets_.end());
}

void LazySynthesisTableau::update_gadgets(
    const MappingFrontier::Operation& operation) {
  // rz_pauli and rx_pauli updated depending on input OpType
  // Only accepts Clifford operations
  OpType type = operation.op->get_type();
  switch (type) {
    case OpType::S: {
      Qubit q(operation.nodes[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.")
      }
      this->tableau_.apply_S_at_end(q);
      break;
    }
    case OpType::V: {
      Qubit q(operation.nodes[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->tableau_.apply_V_at_end(q);
      break;
    }
    case OpType::Z: {
      Qubit q(operation.nodes[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->tableau_.apply_pauli_at_end(QubitPauliTensor(q, Pauli::Z), 2);
      break;
    }
    case OpType::X: {
      Qubit q(operation.nodes[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->tableau_.apply_pauli_at_end(QubitPauliTensor(q, Pauli::X), 2);
      break;
    }
    case OpType::Y {
      Qubit q(operation.nodes[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->tableau_.apply_pauli_at_end(QubitPauliTensor(q, Pauli::Y), 2);
      break;
    }
    case OpType::H: {
      Qubit q(operation.nodes[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.")
      }
      this->tableau_.apply_S_at_end(q);
      this->tableau_.apply_V_at_end(q);
      this->tableau_.apply_S_at_end(q);
      break;
    }
    case OpType::CX: {
      Qubit q_ctrl(operation.nodes[0]);
      if (this->dependencies_.count(q_ctrl) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      Qubit q_trgt(operation.nodes[1]);
      if (this->dependencies_.count(q_trgt) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->tableau_.apply_CX_at_end(q_ctrl, q_trgt);
      this->dependencies_[q_ctrl].insert(q_trgt);
      this->dependencies_[q_trgt].insert(q_ctrl);
      break;
    }
    case OpType::CZ: {
      Qubit q_ctrl(operation.nodes[0]);
      if (this->dependencies_.count(q_ctrl) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      Qubit q_trgt(operation.nodes[1]);
      if (this->dependencies_.count(q_trgt) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->tableau_.apply_S_at_end(q_trgt);
      this->tableau_.apply_V_at_end(q_trgt);
      this->tableau_.apply_S_at_end(q_trgt);
      this->tableau_.apply_CX_at_end(q_ctrl, q_trgt);
      this->tableau_.apply_S_at_end(q_trgt);
      this->tableau_.apply_V_at_end(q_trgt);
      this->tableau_.apply_S_at_end(q_trgt);

      this->dependencies_[q_ctrl].insert(q_trgt);
      this->dependencies_[q_trgt].insert(q_ctrl);
      break;
    }
    default: {
      throw BadOpType(type);
    }
  }
  this->verts.insert(operation.vertex);
}
}  // namespace tket