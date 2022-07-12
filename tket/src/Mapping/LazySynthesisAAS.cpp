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
void LazySynthesisTableau::add_qubit(const Qubit& qubit){
  this->clifford.add_qubit(qb);
  this->rx_pauli.insert({qb, QubitPauliTensor(qb, Pauli::X)});
  this->rz_pauli.insert({qb, QubitPauliTensor(qb, Pauli::Z)});
  this->dependencies.insert({qb, {}});
}

void LazySynthesisTableau::update_gadgets(const unit_vector_t& unitids, const Op_ptr& op_ptr){
  // rz_pauli and rx_pauli updated depending on input OpType
  // Only accepts Clifford operations
  OpType type = op_ptr->get_type();
  switch (type) {
        case OpType::S: {
          Qubit q(args[0]);
          if(this->dependencies.count(q) == 0){
            throw LazySynthesisError("LazySynthesisTableau object does not contain given Qubit.")
          }
          rx_pauli[q] = i_ * rz_pauli[q] * rx_pauli[q];
          break;
        }
        case OpType::V: {
          Qubit q(args[0]);
          if(this->dependencies.count(q) == 0){
            throw LazySynthesisError("LazySynthesisTableau object does not contain given Qubit.");
          }
          rz_pauli[q] = i_ * rx_pauli[q] * rz_pauli[q];
          break;
        }
        case OpType::Z: {
          Qubit q(args[0]);
          if(this->dependencies.count(q) == 0){
            throw LazySynthesisError("LazySynthesisTableau object does not contain given Qubit.");
          }
          rx_pauli[q] = -1. * rx_pauli[q];
          break;
        }
        case OpType::X: {
          Qubit q(args[0]);
          if(this->dependencies.count(q) == 0){
            throw LazySynthesisError("LazySynthesisTableau object does not contain given Qubit.");
          }
          rz_pauli[q] = -1. * rz_pauli[q];
          break;
        }
        case OpType::Sdg: {
          Qubit q(args[0]);
          if(this->dependencies.count(q) == 0){
            throw LazySynthesisError("LazySynthesisTableau object does not contain given Qubit.");
          }
          rx_pauli[q] = -i_ * rz_pauli[q] * rx_pauli[q];
          break;
        }
        case OpType::Vdg: {
          Qubit q(args[0]);
          if(this->dependencies.count(q) == 0){
            throw LazySynthesisError("LazySynthesisTableau object does not contain given Qubit.");
          }
          rz_pauli[q] = -i_ * rx_pauli[q] * rz_pauli[q];
          break;
        }
        case OpType::CX: {
          Qubit q_ctrl(args[0]);
          if(this->dependencies.count(q_ctrl) == 0){
            throw LazySynthesisError("LazySynthesisTableau object does not contain given Qubit.");
          }
          Qubit q_trgt(args[1]);
          if(this->dependencies.count(q_trgt) == 0){
            throw LazySynthesisError("LazySynthesisTableau object does not contain given Qubit.");
          }
          rx_pauli[q_ctrl] = rx_pauli[q_ctrl] * rx_pauli[q_trgt];
          rz_pauli[q_trgt] = rz_pauli[q_ctrl] * rz_pauli[q_trgt];
          break;
        }
        default: {
          throw BadOpType(type);
        }
  }
  clifford.add_op(op_ptr, args);
}

void LazySynthesisTableau::merge_tableau(const LazySynthesisTableau& merge){
  // first check there are no overlapping keys in dependencies which would prevent a merge
  for(const Qubit& q : merge.depencies){
    if(std::find(q, this->depencies.begin(), this->dependencie.end()) == this->dependencies.end()){
      throw LazySynthesisError("Tableau to be merged has overlapping qubits.");
    }
  }
  // merge dependencies maps
  this->dependencies.insert(merge.dependencies.begin(), merge.dependencies.end());
  // merge rz_paulis
  this->rz_paulis.insert(merge.rz_paulis.begin(), merge.rz_paulis.end());
  // merge rz_paulis
  this->rx_paulis.insert(merge.rx_paulis.begin(), merge.rx_paulis.end());
  // append Clifford circuit
  // TODO: We don't want to relabel qubits, and all qubit should be unique, double check this 
  this->clifford.append(merge.clifford);
}


std::pair<QubitPauliTensor, Expr> LazySynthesisTableau::extract_qubit_pauli_tensor(const Qubit& qubit, Op_ptr op_ptr){
  OpType type = op_ptr->get_type();
  switch (type) {
    case OpType::Rz {
      return {rz_pauli[q], op_ptr->get_params()[0]};
    }
    case OpType::Rx {
      return {rx_pauli[q], op_ptr->get_params()[0]};
    }
    default: {
      throw BadOpType(type);
    }
  }
}
 
// std::vector<std::pair<QubitPauliTensor, Expr>> LazySynthesisTableau::extract_qubit(const Qubit& qubit){
//   // Chatting to Will, it's possible to remove a Qubit from a Tableau with 3/4 Pauli gadgets.
//   // We may find that this is a useful result at times
//   // Worth adding and experimenting with
// }


Circuit LazySynthesisTableau::expel_clifford(){
  this->rz_pauli = {};
  this->rx_pauli = {};
  this->dependencies = {};
  return this->clifford;
}

LazyAASRoutingMethod::LazyAASRoutingMethod(){};

std::pair<bool, unit_map_t> LazyAASRoutingMethod::routing_method(
    MappingFrontierPtr& mapping_frontier,
    const ArchitecturePtr& architecture) const {

/*

n.b.  in the below the term "tableau" is used as it presents a common data structure for which the below might be done.
In practice, probably implement something similar to "pairwise_pauli_gadgets"

1) mapping_frontier_->linear_boundary (probably) has some CX gates (or 2-qubit)  which are not architecture permitted
For each of these CX gate:
* if a "tableau" doesn't exist for either Qubit -> optionally start a new 2-qubit "tableau" based on upcoming gates and connectivity of other "tableau"
* if a "tableau" exists for one Qubit in pair -> optionally extend "tableau" with new qubit based on upcoming gates and connectivty of other "tableau", or expel internal Clifford, or if tableau qubit is control then codiagonalise Clifford until the Qubit in tableau is diagonal
* if the same "tableau" exists for both Qubits -> update "tableau"
* if each Qubit has a different "tableau" -> merge "tableau" and update, or expel Clifford of both

2) If there is at least one tableau, start a new process of slice iteration (probably need to copy some logic from advance_next_2qb_slice)

3) For each slice:
* If there is a single qubit gate on a qubit not in a "tableau", skip
* If there is a single qubit Clifford on a qubit in a "tableau", update "tableau", skip edge
* If there is a single qubit non-Clifford on a qubit in a "tableau", perform fan out operation, optionally synthesise full "tableau" (see 5),skip edge
* If there is an architecture permitted CX gate between a pair of qubits not in a "tableau", skip edge
* If there is an architecture permitted CX gate between a pair of qubits both in the same "tableau", update "tableau" and skip edge
* If there is an architecture permitted CX gate between a pair of qubits, one in a tableau, one not, update "tableau" with new qubit and CX operation, skip edge, or, if tableau is control potentially codiagonalise up to that qubit and commute, or expel Qubit from Tableau (ask Will i.e. Pauli Graph stuff)
* If there is an architecture permitted CX gate between a pair of qubits in different "tableau", merge tableau, update "tableau", skip edge
* If there is an architecture non-permitted CX gate, freeze edges (advance_frontier_boundary will find them later and utilise step 1)

4) Once every edge is "frozen", exit method and allow Routing to continue

5) For synthesis:
There is a (metaphorically kind of?) continuous set of choices that can be made. Things to consider:
* The set of other interacting physical qubits is known by asking other held "tableau"
* "fan out" gates for a single "tableau" will likely utilise Qubits in other "tableau". We can avoid a "fan in" if we can leave these in arbitrary states, 
so probably some decision making is required as to what physical qubit path to "fan out" along
* At any point in a "fan out" we can replace the "fan out" operation with a SWAP gate, based on the interaction graph of all other "tableau"
* Note that a SWAP operation could also be utilised to avoid "fan out" gates on Qubit in other "tableau"
Or to summarise, the decision made here should attempt to be "optimal" for all known "tableau"

6) 
It is (highly) possible that the end of the Circuit is reached without some Tableau synthesised i.e. we have some final Clifford operator tagged on the end of the Circuit.
In conjunction with this work we should in some capacity (maybe by a new Circuit attribute?) update TKET to deal with how to classically emulate a final non-trival Clifford operator
(Section A: https://quantum-journal.org/papers/q-2022-06-07-729/pdf/).
Given this, (paraphrased) a co-diagonalization circuit can be added to the circuit with a maximum of three layers of two-qubit gates which then leaves this final layer as diagonal 
operators with an effect we can fix when sampling bit strings.
*/



        
  for (auto it =
           mapping_frontier_->linear_boundary->get<TagKey>().begin();
       it != mapping_frontier_->linear_boundary->get<TagKey>().end();
       ++it) {


       }
    }

nlohmann::json LazyAASRoutingMethod::serialize() const {
  nlohmann::json j;
  j["name"] = " LazyAASRoutingMethod";
  return j;
}

LazyAASRoutingMethod LazyAASRoutingMethod::deserialize(
    const nlohmann::json& j) {
  return LazyAASRoutingMethod();
}

}  // namespace tket