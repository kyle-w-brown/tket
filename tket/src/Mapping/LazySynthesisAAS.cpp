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

// Methods for "LazySynthesisTableau"
void LazySynthesisTableau::add_qubit(const Qubit& qubit) {
  this->clifford_.add_qubit(qb);
  this->rx_pauli_.insert({qb, QubitPauliTensor(qb, Pauli::X)});
  this->rz_pauli_.insert({qb, QubitPauliTensor(qb, Pauli::Z)});
  this->dependencies_.insert({qb, {}});
}

void LazySynthesisTableau::update_gadgets(
    const unit_vector_t& unitids, const Op_ptr& op_ptr) {
  // rz_pauli and rx_pauli updated depending on input OpType
  // Only accepts Clifford operations
  OpType type = op_ptr->get_type();
  switch (type) {
    case OpType::S: {
      Qubit q(args[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.")
      }
      this->rx_pauli_[q] = i_ * this->rz_pauli__[q] * this->rx_pauli_[q];
      break;
    }
    case OpType::V: {
      Qubit q(args[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->rz_pauli_[q] = i_ * this->rx_pauli_[q] * this->rz_pauli_[q];
      break;
    }
    case OpType::Z: {
      Qubit q(args[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->rx_pauli_[q] = -1. * this->rx_pauli_[q];
      break;
    }
    case OpType::X: {
      Qubit q(args[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->rz_pauli_[q] = -1. * this->rz_pauli_[q];
      break;
    }
    case OpType::Sdg: {
      Qubit q(args[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->rx_pauli_[q] = -i_ * this->rz_pauli_[q] * this->rx_pauli_[q];
      break;
    }
    case OpType::Vdg: {
      Qubit q(args[0]);
      if (this->dependencies_.count(q) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->rz_pauli_[q] = -i_ * this->rx_pauli_[q] * this->rz_pauli_[q];
      break;
    }
    case OpType::CX: {
      Qubit q_ctrl(args[0]);
      if (this->dependencies_.count(q_ctrl) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      Qubit q_trgt(args[1]);
      if (this->dependencies_.count(q_trgt) == 0) {
        throw LazySynthesisError(
            "LazySynthesisTableau object does not contain given Qubit.");
      }
      this->rx_pauli_[q_ctrl] =
          this->rx_pauli_[q_ctrl] * this->rx_pauli_[q_trgt];
      this->rz_pauli_[q_trgt] =
          this->rz_pauli_[q_ctrl] * this->rz_pauli_[q_trgt];
      break;
    }
    default: {
      throw BadOpType(type);
    }
  }
  this->clifford_.add_op(op_ptr, args);
}

// TODO overload + operator to do this, purely for fun
void LazySynthesisTableau::merge_tableau(const LazySynthesisTableau& merge) {
  // first check there are no overlapping keys in dependencies which would
  // prevent a merge
  // TODO: overload some intersection operator
  for (const Qubit& q : merge.depencies) {
    if (std::find(q, this->dependencies_.begin(), this->dependencies_.end()) ==
        this->dependencies_.end()) {
      throw LazySynthesisError("Tableau to be merged has overlapping qubits.");
    }
  }
  // merge dependencies maps
  this->dependencies_.insert(
      merge.dependencies_.begin(), merge.dependencies_.end());
  this->pauli_gadgets_.insert(
      this->pauli_gadgets_.end(), merge.pauli_gadgets_.begin(),
      merge.pauli_gadgets_.end());
  this->q_in_hole.insert(merge.q_in_hole.begin(), merge.q_in_hole.end());

  // TOOD: Update with tableau
  // merge rz_paulis
  this->this->rz_pauli_.insert(merge.rz_pauli_.begin(), merge.rz_pauli_.end());
  // merge rz_paulis
  this->rx_pauli_.insert(merge.rx_pauli_.begin(), merge.rx_pauli_.end());
  // append Clifford circuit
  // TODO: We don't want to relabel qubits, and all qubit should be unique,
  // double check this
  this->clifford_.append(merge.clifford_);
}

std::pair<QubitPauliTensor, Expr>
LazySynthesisTableau::extract_qubit_pauli_tensor(
    const Qubit& qubit, Op_ptr op_ptr) {
  OpType type = op_ptr->get_type();
  switch (type) {
    case OpType::Rz {
      return {this->rz_pauli_[q], op_ptr->get_params()[0]};
    } case OpType::Rx {
      return {this->rx_pauli_[q], op_ptr->get_params()[0]};
    } default: {
      throw BadOpType(type);
    }
  }
}

// Circuit LazySynthesisTableau::expel_clifford() {
//   this->rz_pauli_ = {};
//   this->rx_pauli_ = {};
//   this->dependencies_ = {};
//   return this->clifford_;
// }

LazyCliffordRoutingMethod::LazyCliffordRoutingMethod() {}

bool LazyCliffordRoutingMethod::update_from_mapping_frontier(
    MappingFrontierPtr& mapping_frontier, const ArchitecturePtr& architecture,
    bool allow_non_adjacent_cx) {
  bool modified = false;
  for (auto it = mapping_frontier_->linear_boundary->get<TagKey>().begin();
       it != mapping_frontier_->linear_boundary->get<TagKey>().end(); ++it) {
    Edge e0 = this->mapping_frontier_->circuit_.get_nth_out_edge(
        it->second.first, it->second.second);
    Vertex v0 = this->mapping_frontier_->circuit_.target(e0);
    Op_ptr op = this->mapping_frontier_->circuit_.get_Op_ptr_from_Vertex(v0);
    OpType type = op->get_type();
    // for a first implementation, only support:
    // S, V, Z, X, Sdg, Vdg, CX, Rz, Rx
    // And throw error if this isn't supported
    // Can extend this range of operations once we've benchmarked etc
    switch (type) {
      case OpType::S:
      case OpType::V:
      case OpType::Z:
      case OpType::X:
      case OpType::Sdg:
      case OpType::Vdg: {
        for (const LazySynthesisTableau& t : this->tableau_) {
          if (t->dependencies_.count(it->first) != 0) {
            t->update_gadgets({it->first}, op);
            modified = true;
            break;
          }
        }
      }
      case OpType::Rz:
      case OpType::Rx: {
        for (const LazySynthesisTableau& t : this->tableau_) {
          if (t->dependencies_.count(it->first) != 0) {
            std::pair<QubitPauliTensor, Expr> qpt_e =
                t->extract_qubit_pauli_tensor(it->first, op);
            modified = true;

            /*
            do some things with the qubit pauli tensor i.e. the routing bit...
            */
            break;
          }
        }
      }
      case OpType::CX {
        auto jt = it; ++jt;
        while (jt != mapping_frontier_->linear_boundary->get<TagKey>().end()) {
          Edge e1 = this->mapping_frontier_->circuit_.get_nth_out_edge(
              jt->second.first, jt->second.second);
          Vertex v1 = this->mapping_frontier_->circuit_.target(e1);
          if (v0 == v1) {
            // TODO: currently we always amend a tableau if we hit a CX gate
            // change this condition in future
            modified = true;

            Qubit qubit0(it->first);
            Qubit qubit1(jt->first);

            // We assume that in outer mapping loop, non architecture permitted
            // gates will always be present at the frontier
            // We update the tableau for operations that would be missed in this
            // model
            // TODO: this means reiterating over the same vertices often - can
            // we update RV3 structure to advance and update data structures
            // held by internal Routing Method?
            // probably yeah though likely some solid plumbing - leave for
            // future
            if (!architecture->valid_operation({Node(qubit0), Node(qubit1)}) &&
                !allow_non_adjacent_cx) break;

            unsigned port0 = it->second.second;
            unsigned port1 = jt->second.second;

            auto kt = this->tableau_.begin();
            while (kt != this->tableau_.end() &&
                   kt->dependencies_.count(q0) != 1)++ kt;

            auto lt = this->tableau_.end();
            while (lt != this->tableau_.end() &&
                   lt->dependencies_.count(q1) != 1)++ lt;

            auto update_cx_gadget =
                [&kt, &op, &qubit0, &qubit1, &port0, &port1 ]() {
                  // => qubit0 is control
                  if (port0 < port1) {
                    kt->update_gadgets({qubit0, qubit1}, op);
                  }
                  // => qubit1 is control
                  else {
                    kt->update_gadgets({qubit1, qubit0}, op);
                  }
                };

            if (kt == this->tableau_.end()) {
              // => neither qubit in a tableau
              if (lt == this->tableau_.end()) {
                // TODO: think about several options:
                // ignore gate and route as usual
                // For now:
                // make a new tableau
                LazySynthesisTableau lst; lst.add_qubit(qubit0);
                lst.add_qubit(qubit1);
                kt = this->tableau_.insert(lst);
                update_cx_gadget();
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
                lt->add_qubit(q0); update_cx_gadget();
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
                kt->add_qubit(q1); update_cx_gadget();
              } else {
                // => both qubits in same tableau, so just add gate
                if (lt == kt) {
                  update_cx_gadget();
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
                  kt->merge_tableau(*lt); update_cx_gadget();
                  // remove merged tableau
                  this->tableau_.erase(lt);
                }
              }
            }
          } ++jt;
        }
      } default: {
        throw BadOpType(type);
      }
    }
  }
}

// TODO:
// Helper method, code taken from MappingFrontier.cpp
// make accessible without copying
std::shared_ptr<unit_frontier_t> frontier_convert_vertport_to_edge(
    const Circuit& circuit,
    const std::shared_ptr<unit_vertport_frontier_t>& u_frontier) {
  // make empty unit_frontier_t object
  std::shared_ptr<unit_frontier_t> output_frontier =
      std::make_shared<unit_frontier_t>();
  // iterate through u_frontier, convert VertPort to Edge and insert
  for (const std::pair<UnitID, VertPort>& pair : u_frontier->get<TagKey>()) {
    output_frontier->insert(
        {pair.first,
         circuit.get_nth_out_edge(pair.second.first, pair.second.second)});
  }
  return output_frontier;
}
/*

n.b.  in the below the term "tableau" is used as it presents a common data
structure for which the below might be done. In practice, probably implement
something similar to "pairwise_pauli_gadgets"

1) mapping_frontier_->linear_boundary (probably) has some CX gates (or 2-qubit)
which are not architecture permitted For each of these CX gate:
* if a "tableau" doesn't exist for either Qubit -> optionally start a new
2-qubit "tableau" based on upcoming gates and connectivity of other "tableau"
* if a "tableau" exists for one Qubit in pair -> optionally extend "tableau"
with new qubit based on upcoming gates and connectivty of other "tableau", or
expel internal Clifford, or if tableau qubit is control then codiagonalise
Clifford until the Qubit in tableau is diagonal
* if the same "tableau" exists for both Qubits -> update "tableau"
* if each Qubit has a different "tableau" -> merge "tableau" and update, or
expel Clifford of both

2) If there is at least one tableau, start a new process of slice iteration
(probably need to copy some logic from advance_next_2qb_slice)

3) For each slice:
* If there is a single qubit gate on a qubit not in a "tableau", skip
* If there is a single qubit Clifford on a qubit in a "tableau", update
"tableau", skip edge
* If there is a single qubit non-Clifford on a qubit in a "tableau", perform fan
out operation, optionally synthesise full "tableau" (see 5),skip edge
* If there is an architecture permitted CX gate between a pair of qubits not in
a "tableau", skip edge
* If there is an architecture permitted CX gate between a pair of qubits both in
the same "tableau", update "tableau" and skip edge
* If there is an architecture permitted CX gate between a pair of qubits, one in
a tableau, one not, update "tableau" with new qubit and CX operation, skip edge,
or, if tableau is control potentially codiagonalise up to that qubit and
commute, or expel Qubit from Tableau (ask Will i.e. Pauli Graph stuff)
* If there is an architecture permitted CX gate between a pair of qubits in
different "tableau", merge tableau, update "tableau", skip edge
* If there is an architecture non-permitted CX gate, freeze edges
(advance_frontier_boundary will find them later and utilise step 1)

4) Once every edge is "frozen", exit method and allow Routing to continue

5) For synthesis:
There is a (metaphorically kind of?) continuous set of choices that can be made.
Things to consider:
* The set of other interacting physical qubits is known by asking other held
"tableau"
* "fan out" gates for a single "tableau" will likely utilise Qubits in other
"tableau". We can avoid a "fan in" if we can leave these in arbitrary states, so
probably some decision making is required as to what physical qubit path to "fan
out" along
* At any point in a "fan out" we can replace the "fan out" operation with a SWAP
gate, based on the interaction graph of all other "tableau"
* Note that a SWAP operation could also be utilised to avoid "fan out" gates on
Qubit in other "tableau" Or to summarise, the decision made here should attempt
to be "optima l" for all known "tableau"

6)
It is (highly) possible that the end of the Circuit is reached without some
Tableau synthesised i.e. we have some final Clifford operator tagged on the end
of the Circuit. In conjunction with this work we should in some capacity (maybe
by a new Circuit attribute?) update TKET to deal with how to classically emulate
a final non-trival Clifford operator (Section A:
https://quantum-journal.org/papers/q-2022-06-07-729/pdf/). Given this,
(paraphrased) a co-diagonalization circuit can be added to the circuit with a
maximum of three layers of two-qubit gates which then leaves this final layer as
diagonal operators with an effect we can fix when sampling bit strings.
*/

std::pair<bool, unit_map_t> LazyCliffordRoutingMethod::routing_method(
    MappingFrontierPtr& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  bool modified =
      this->update_from_mapping_frontier(mapping_frontier, architecture, true);
  // => that there is at least one tableau so we need to start checking

  if (this->tableau_.size() > 0) {
    // std::shared_ptr<unit_vertport_frontier_t> linear_boundary =
    // mapping_frontier_->linear_boundary; std::shared_ptr<unit_frontier_t>
    // l_frontier_edges =
    // frontier_convert_vertport_to_edge(mapping_frontier->circuit_,);
    // CutFrontier next_cut =
    //     this->circuit_.next_cut(l_frontier_edges, this->boolean_boundary);

    // get next slice of gates from edges
    // update tableau
    //
    //
    // repeat until
  }
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