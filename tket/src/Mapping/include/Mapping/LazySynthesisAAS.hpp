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


class LazySynthesisTableau {
  public:
  
  LazySynthesisTableau() {};

  add_qubit(const Qubit& qubit);
  




  protected;

  Circuit clifford;
  std::map<Qubit, QubitPauliTensor> rx_pauli;
  std::map<Qubit, QubitPauliTensor> rz_pauli;
  std::map<Qubit, std::set<Qubit>> dependencies;
}



class LazyArchitectureAwareSynthesis {
  public: 

  LazyArchitectureAwareSynthesis();



}



class LazyAASRoutingMethod : public RoutingMethod {
 public:
  /**
   *
   *
   */
  LazyAASRoutingMethod();

  /**
   *
   */
  std::pair<bool, unit_map_t> routing_method(
      MappingFrontierPtr& mapping_frontier,
      const ArchitecturePtr& architecture) const override;

  nlohmann::json serialize() const override;

  static LexiRouteRoutingMethod deserialize(const nlohmann::json& j);
};

JSON_DECL(LazyAASRoutingMethod);

}  // namespace tket
