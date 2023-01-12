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

#include "Command.hpp"
#include "Utils/Json.hpp"

namespace tket {

void to_json(nlohmann::json& j, const Command& com) {
  const Op_ptr op = com.get_op_ptr();
  const std::optional<std::string> opgroup = com.get_opgroup();
  j["op"] = op;
  if (opgroup) {
    j["opgroup"] = opgroup.value();
  }

  const op_signature_t& sig = op->get_signature();
  const unit_vector_t& args = com.get_args();

  nlohmann::json j_args;
  for (size_t i = 0; i < sig.size(); i++) {
    if (sig[i] == EdgeType::Quantum) {
      const Qubit& qb = static_cast<const Qubit&>(args[i]);
      j_args.push_back(qb);
    } else if (sig[i] == EdgeType::Classical || sig[i] == EdgeType::Boolean) {
      const Bit& cb = static_cast<const Bit&>(args[i]);
      j_args.push_back(cb);
    } else {
      std::cout << "found wasm edge\n";
    }
  }

  j["args"] = j_args;
}
void from_json(const nlohmann::json& j, Command& com) {
  std::cout << "command from json - 0\n";
  const auto op = j.at("op").get<Op_ptr>();
  std::cout << "command from json - 1\n";
  std::optional<std::string> opgroup;
  std::cout << "command from json - 2\n";
  if (j.contains("opgroup")) {
    opgroup = j.at("opgroup").get<std::string>();
  }
  std::cout << "command from json - 3\n";
  const op_signature_t& sig = op->get_signature();
  std::cout << "command from json - 4\n";
  const nlohmann::json& j_args = j.at("args");
  std::cout << "command from json - 5\n";
  if (sig.size() != j_args.size()) {
    JsonError("Number of args does not match signature of op.");
  }
  std::cout << "command from json - 6\n";
  unit_vector_t args;
  std::cout << "command from json - 7\n";
  for (size_t i = 0; i < sig.size(); i++) {
    std::cout << "command from json - 8 - x\n";
    if (sig[i] == EdgeType::Quantum) {
      args.push_back(j_args[i].get<Qubit>());
    } else if (sig[i] == EdgeType::Classical || sig[i] == EdgeType::Boolean) {
      args.push_back(j_args[i].get<Bit>());
    } else {
      std::cout << "found wasm edge\n";
    }
  }

  com = Command(op, args, opgroup);
}

}  // namespace tket
