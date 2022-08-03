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

bool LazyCliffordRoutingMethod::update_from_passed_operations(
    MappingFrontierPtr& mapping_frontier, const ArchitecturePtr& architecture) {
  // i.e. we have at least one tableau
  if (this->tableau_qubits_.size() > 0) {
    for (const MappingFrontier::Operation& operation :
         mapping_frontier->passed_operations) {
      OpType type = operation.second->get_type();
      if (type != OpType::Barrier) {
        // for a first implementation, only support:
        // S, V, Z, X, Sdg, Vdg, CX, CZ, H, Rz, Rx
        // And throw error if this isn't supported
        // Can extend this range of operations once we've benchmarked etc
        switch (type) {
          case OpType::S:
          case OpType::V:
          case OpType::Z:
          case OpType::X:
          case OpType::Sdg:
          case OpType::Vdg:
          case OpType::H: {
            for (const LazySynthesisTableau& t : this->tableau_) {
              if (t->dependencies_.count(operation.nodes[0]) != 0) {
                t->update_gadgets(operation);
                modified = true;
                break;
              }
            }
          }
          case OpType::Rz:
          case OpType::Rx: {
            for (const LazySynthesisTableau& t : this->tableau_) {
              if (t->dependencies_.count(operation.nodes[0]) != 0) {
                t->pauli_gadgets_.push_back(t->extract_qubit_pauli_tensor(
                    operation.nodes[0], operation.op));
                modified = true;
                break;
              }
            }
          }
          case OpType::CX:
          case OpType::CZ: {
            // TODO: currently we always amend a tableau if we hit a CX/CZ gate
            // change this condition in future
            modified = true;

            Qubit qubit0(operation.nodes[0]);
            Qubit qubit1(operation.nodes[1]);

            // see if gate is in tableau
            // any operation encountered should be architecturally valid
            auto kt = this->tableau_.begin();
            while (kt != this->tableau_.end() &&
                   kt->dependencies_.count(q0) != 1) {
              ++kt;
            }

            auto lt = this->tableau_.end();
            while (lt != this->tableau_.end() &&
                   lt->dependencies_.count(q1) != 1) {
              ++lt;
            }

            if (kt == this->tableau_.end()) {
              // n.b if lt == end => neither qubit in a tableau
              // therefore just pass
              // this condition is:
              // => q1 is in tableau, q0 isn't
              if (lt != this->tableau_.end()) {
                // TODO: think about several options:
                // remove q1 from tableau via gadgets and make a new tableau
                // remove q1 from tableau via gadgets and route gate as usual
                // if q1 is control, codiagonalise circuit s.t. gate can commute
                // through to start and route as usual
                // For now:
                // Add new qubit, update tableau
                lt->add_qubit(q0);
                kt->update_gadgets(
                    {{qubit0, qubit1}, operation.op, operation.vertex});
              }
            } else {
              // => q0 is in tableau, q1 isn't
              if (lt == this->tableau_.end()) {
                // TODO: think about several options:
                // remove q0 from tableau via gadgets and make a new tableau
                // remove q0 from tableau via gadgets and route gate as usual
                // if q0 is control, codiagonalise circuit s.t. gate can commute
                // through to start and route as usual
                // For now:
                // Add new qubit, update tableau
                kt->add_qubit(q1);
                kt->update_gadgets(
                    {{qubit0, qubit1}, operation.op, operation.vertex});
              } else {
                // => both qubits in same tableau, so just add gate
                if (lt == kt) {
                  kt->update_gadgets(
                      {{qubit0, qubit1}, operation.op, operation.vertex});
                }
                // => qubits are in different tableau
                else {
                  // TODO: think about several options:
                  // remove one qubit from tableau with gadgets and add to other
                  // remove both qubits from tableau with gadgets and make new
                  // tableau
                  // remove both qubits from tableau with gadgets and route as
                  // normal
                  // For now:
                  // merge tableau and add gate
                  kt->merge_tableau(*lt);
                  kt->update_gadgets(
                      {{qubit0, qubit1}, operation.op, operation.vertex});
                  // remove merged tableau
                  this->tableau_.erase(lt);
                }
              }
            }
            break;
          }
          default: {
            throw BadOpType(type);
          }
        }
      }
    }
  }
  return modified;
}

bool LazyCliffordRoutingMethod::update_from_mapping_frontier(
    MappingFrontierPtr& mapping_frontier, const ArchitecturePtr& architecture) {
  bool modified = false;
  // if tableau are made, move boundary on

  std::map<UnitID, VertPort> replacement_edges;
  // iterate through edges
  for (auto it = mapping_frontier_->linear_boundary->get<TagKey>().begin();
       it != mapping_frontier_->linear_boundary->get<TagKey>().end(); ++it) {
    Edge e0 = this->mapping_frontier_->circuit_.get_nth_out_edge(
        it->second.first, it->second.second);
    Vertex v0 = this->mapping_frontier_->circuit_.target(e0);
    Op_ptr op = this->mapping_frontier_->circuit_.get_Op_ptr_from_Vertex(v0);
    OpType type = op->get_type();

    // We update given tableau by looking at "passed operations" after advancing
    // frontier We make new tableau by looking at gates at frontier boundary We
    // manually move the boundary past a gate if it's added
    if ((type == OpType::CX || type == OpType::CZ) &&
        architecture->node_exists(Node(it->first))) {
      auto jt = it;
      ++jt;
      // iterate through edges
      while (jt != mapping_frontier_->linear_boundary->get<TagKey>().end()) {
        Edge e1 = this->mapping_frontier_->circuit_.get_nth_out_edge(
            jt->second.first, jt->second.second);
        Vertex v1 = this->mapping_frontier_->circuit_.target(e1);
        if (v0 == v1 && architecture->node_exists(Node(jt->first))) {
          // TODO: currently we always amend a tableau if we hit a CX/CZ gate
          // change this condition in future
          modified = true;

          Qubit qubit0(it->first);
          Qubit qubit1(jt->first);

          unsigned port0 = it->second.second;
          unsigned port1 = jt->second.second;

          auto kt = this->tableau_.begin();
          while (kt != this->tableau_.end() && kt->dependencies_.count(q0) != 1)
            ++kt;

          auto lt = this->tableau_.end();
          while (lt != this->tableau_.end() && lt->dependencies_.count(q1) != 1)
            ++lt;

          auto update_2q_gadget = [&kt, &op, &qubit0, &qubit1, &port0,
                                   &port1]() {
            // TODO: Doesn't matter for CZ case, though there might be logic for
            // improving still keep this way for now as easy
            // => qubit0 is control
            if (port0 < port1) {
              kt->update_gadgets(
                  {{qubit0, qubit1}, operation.op, operation.vertex});
            }
            // => qubit1 is control
            else {
              kt->update_gadgets(
                  {{qubit1, qubit0}, operation.op, operation.vertex});
            }
          };

          if (kt == this->tableau_.end()) {
            // => neither qubit in a tableau
            if (lt == this->tableau_.end()) {
              // TODO: think about several options:
              // ignore gate and route as usual
              // For now:
              // make a new tableau
              LazySynthesisTableau lst;
              lst.add_qubit(qubit0);
              lst.add_qubit(qubit1);
              lst.q_ kt = this->tableau_.insert(lst);
              update_2q_gadget();

              replacement_edges.insert(
                  {qubit0,
                   this->mapping_frontier_->circuit_.get_next_edge(v0, e0)});
              replacement_edges.insert(
                  {qubit1,
                   this->mapping_frontier_->circuit_.get_next_edge(v1, e1)});

            }
            // => q1 is in tableau, q0 isn't
            else {
              // TODO: think about several options:
              // remove q1 from tableau via gadgets and make a new tableau
              // remove q1 from tableau via gadgets and route gate as usual
              // if q1 is control, codiagonalise circuit s.t. gate can commute
              // through to start and route as usual
              // For now:
              // Add new qubit, update tableau
              lt->add_qubit(q0);
              update_2q_gadget();
              replacement_edges.insert(
                  {qubit0,
                   this->mapping_frontier_->circuit_.get_next_edge(v0, e0)});
              replacement_edges.insert(
                  {qubit1,
                   this->mapping_frontier_->circuit_.get_next_edge(v1, e1)});
            }
          } else {
            // => q0 is in tableau, q1 isn't
            if (lt == this->tableau_.end()) {
              // TODO: think about several options:
              // remove q0 from tableau via gadgets and make a new tableau
              // remove q0 from tableau via gadgets and route gate as usual
              // if q0 is control, codiagonalise circuit s.t. gate can commute
              // through to start and route as usual
              // For now:
              // Add new qubit, update tableau
              kt->add_qubit(q1);
              update_2q_gadget();
              replacement_edges.insert(
                  {qubit0,
                   this->mapping_frontier_->circuit_.get_next_edge(v0, e0)});
              replacement_edges.insert(
                  {qubit1,
                   this->mapping_frontier_->circuit_.get_next_edge(v1, e1)});
            } else {
              // => both qubits in same tableau, so just add gate
              if (lt == kt) {
                update_2q_gadget();
                replacement_edges.insert(
                    {qubit0,
                     this->mapping_frontier_->circuit_.get_next_edge(v0, e0)});
                replacement_edges.insert(
                    {qubit1,
                     this->mapping_frontier_->circuit_.get_next_edge(v1, e1)});
              }
              // => qubits are in different tableau
              else {
                // TODO: think about several options:
                // remove one qubit from tableau with gadgets and add to other
                // remove both qubits from tableau with gadgets and make new
                // tableau
                // remove both qubits from tableau with gadgets and route as
                // normal
                // For now:
                // merge tableau and add gate
                kt->merge_tableau(*lt);
                update_2q_gadget();
                // remove merged tableau
                this->tableau_.erase(lt);
                replacement_edges.insert(
                    {qubit0,
                     this->mapping_frontier_->circuit_.get_next_edge(v0, e0)});
                replacement_edges.insert(
                    {qubit1,
                     this->mapping_frontier_->circuit_.get_next_edge(v1, e1)});
              }
            }
          }
        }
        ++jt;
      }
    }
  }
  // for each edge added to a tableau, advance mapping_frontier past them
  // two reasons:
  // 1) advance_frontier_boundary would consider them invalid and not pass them
  // 2) avoid adding twice from passed_operations
  for (const std::pair<UnitID, VertPort>& uid_vp : replacement_edges) {
    mapping_frontier_->linear_boundary->replace(
        this->linear_boundary->get<TagKey>().find(uid_vp.first),
        {uid_vp.first, uid_vp.second});
  }
  return modified;
}

LazyCliffordRoutingMethod::LazyCliffordRoutingMethod() {}

std::pair<bool, unit_map_t> LazyCliffordRoutingMethod::routing_method(
    MappingFrontierPtr& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  bool modified_passed_operation =
      this->update_from_passed_operations(mapping_frontier, architecture);
  bool modified_mapping_frontier =
      this->update_from_mapping_frontier(mapping_frontier, architecture);

  auto check_finish = [&mapping_frontier]() {
    for (const std::pair<UnitID, VertPort>& pair :
         mapping_frontier->linear_boundary->get<TagKey>()) {
      Edge e = mapping_frontier->circuit_.get_nth_out_edge(
          pair.second.first, pair.second.second);
      Vertex v = mapping_frontier->circuit_.target(e);
      OpType ot = mapping_frontier->circuit_.get_OpType_from_Vertex(v);
      if (!is_final_q_type(ot) && ot != OpType::ClOutput) {
        return false;
      }
    }
    return true;
  };

  if (check_finish) {
    // decompose_all_gadgets
  }

  return {modified_passed_operation || modified_mapping_frontier, {}};
}

nlohmann::json LazyCliffordRoutingMethod::serialize() const {
  nlohmann::json j;
  j["name"] = "LazyCliffordRoutingMethod";
  return j;
}

LazyCliffordRoutingMethod LazyCliffordRoutingMethod::deserialize(
    const nlohmann::json& j) {
  return LazyCliffordRoutingMethod();
}

}  // namespace tket