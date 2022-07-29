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

#pragma once

#include "Mapping/RoutingMethod.hpp"

namespace tket {

class LazyCliffordSynthesisError : public std::logic_error {
 public:
  explicit LazyCliffordSynthesisError(const std::string& message)
      : std::logic_error(message) {}
};

class LazySynthesisTableau {
 public:
  LazySynthesisTableau(){};

  void add_qubit(const Qubit& qubit, const Edge& edge);

  void update_gadgets(const unit_vector_t& unitids, const Op_ptr& op_ptr);
  /*
  TODO: Will has some ideas for an optimal number (I believe Will said 3/4
  but WIP) of pauli gadgets that need to be inserted to allow a qubit removal
  from a tableau This is a degree of freedom that could provide gains in the
  future remove_qubit(const Qubtit& qubit);
  */

 protected:
  UnitaryTableau tableau_;
  std::map<Qubit, std::set<Qubit>> dependencies_;

  // LST consumes vertices as it's passed them, store and replace later
  std::map<Qubit, Edge> q_in_hole;
  VertexSet verts;
}

class LazyCliffordRoutingMethod : public RoutingMethod {
 public:
  /**
   *
   *
   */
  LazyCliffordRoutingMethod();

  /**
   *
   */
  std::pair<bool, unit_map_t> routing_method(
      MappingFrontierPtr& mapping_frontier,
      const ArchitecturePtr& architecture) const override;

  nlohmann::json serialize() const override;

  static LexiRouteRoutingMethod deserialize(const nlohmann::json& j);

 private:
  bool update_from_mapping_frontier(
      MappingFrontierPtr& mapping_frontier,
      const ArchitecturePtr& architecture);
  std::set<LazySynthesisTableau> tableau_;
  std::set<Qubit> tableau_qubits_;
};

JSON_DECL(LazyCliffordRoutingMethod);

}  // namespace tket
