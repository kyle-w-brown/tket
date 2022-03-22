# Copyright 2019-2022 Cambridge Quantum Computing
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from enum import Enum
from pytket.mapping import get_token_swapping_network  # type: ignore
from pytket.circuit import Circuit, Bit, Qubit, OpType  # type: ignore
from pytket.predicates import CompilationUnit  # type: ignore
from pytket.passes import RoutingPass  # type: ignore
from pytket.architecture import Architecture, FullyConnected  # type: ignore
from typing import Dict, Tuple, List, Union


class Knitter(Enum):
    @staticmethod
    def separate(
        pattern_list: List[Tuple[Circuit, Dict[str, Dict[Qubit, int]]]],
        arch: Architecture,
    ) -> Tuple[Circuit, Dict[str, Dict[Qubit, int]]]:
        """
        Takes a list of tuples containing pytket circuits and their corresponding
        i/o maps as well as an architecture. It then routes each circuit on the
        architecture separately. Finally it connects the outputs of one segment to
        the inputs of the next by generating a network of SWAP gates. Then returns
        a tuple containing the new circuit object and its corresponding i/o map.
        
        :param pattern_list:  List of tuples of circuits with their i/o maps.
        :param type:          List[Tuple[Circuit,Dict[str, Dict[Qubit, int]]]]
        
        :param arch:          The architecture onto which to fit the circuits.
        :param type:          Architecture
        
        :returns:        A tuple containing a circuit and an i/o map dictionary.
        :rtype:          Tuple[Circuit,Dict[str, Dict[Qubit, int]]]
        """
        new_c = Circuit()
        for q in arch.nodes:
            new_c.add_qubit(q)
        for p in range(len(pattern_list)):
            (curr_circ, curr_map) = pattern_list[p]
            route_map = {}
            if p == 0:
                for k in curr_map["i"].keys():
                    if k in arch.nodes:
                        temp_val = curr_map["i"][k]
                        route_map[curr_circ.qubits[temp_val]] = k
                        curr_map["i"][k] = k
                        if curr_map["o"][k] == temp_val:
                            curr_map["o"][k] = k
                        else:
                            curr_map["o"][k] = curr_circ.qubits[curr_map["o"][k]]
                    else:
                        curr_map["i"][k] = curr_circ.qubits[curr_map["i"][k]]
                        curr_map["o"][k] = curr_circ.qubits[curr_map["o"][k]]
                curr_circ.rename_units(route_map)
            else:
                for k in curr_map["i"].keys():
                    curr_map["i"][k] = curr_circ.qubits[curr_map["i"][k]]
                    curr_map["o"][k] = curr_circ.qubits[curr_map["o"][k]]
            cu = CompilationUnit(curr_circ)
            if arch == FullyConnected:
                raise TypeError(
                    "This method doesn't work with a FullyConnected architecture."
                )
            RoutingPass(arch).apply(cu)
            used_nodes = set()
            all_nodes = set(arch.nodes)
            unassigned_qubits = []
            for k in curr_map["i"].keys():
                if cu.initial_map[curr_map["i"][k]] in all_nodes:
                    used_nodes |= {cu.initial_map[curr_map["i"][k]]}
                    used_nodes |= {cu.initial_map[curr_map["o"][k]]}
                    curr_map["i"][k] = cu.initial_map[curr_map["i"][k]]
                    curr_map["o"][k] = cu.final_map[curr_map["o"][k]]
                else:
                    unassigned_qubits.append(k)
            for command in cu.circuit.get_commands():
                if command.op.get_name() == "CZ":
                    qubits = command.qubits
                    for q in qubits:
                        if q in used_nodes:
                            used_nodes |= set(qubits)
                            break
                elif command.op.get_name() == "SWAP":
                    qubits = command.qubits
                    if (
                        (qubits[0] in used_nodes) and not (qubits[1] in used_nodes)
                    ) or (not (qubits[0] in used_nodes) and (qubits[1] in used_nodes)):
                        used_nodes ^= set(qubits)
                elif command.op.get_name() == "BRIDGE":
                    qubits = command.qubits
                    for q in {qubits[0], qubits[2]}:
                        if q in used_nodes:
                            used_nodes |= {qubits[0], qubits[2]}
                            break
            permutation = {x: x for x in all_nodes}
            for com in cu.circuit.commands_of_type(OpType.SWAP):
                permutation[com.qubits[0]] = com.qubits[1]
                permutation[com.qubits[1]] = com.qubits[0]
            unused_nodes = all_nodes - used_nodes
            segment_circuit = cu.circuit.copy()
            for uq in unassigned_qubits:
                temp = unused_nodes.pop()
                curr_map["o"][uq] = temp
                for k in permutation.keys():
                    if permutation[k] == temp:
                        segment_circuit.rename_units({curr_map["i"][uq]: k})
                        curr_map["i"][uq] = k
                        break
            new_tuple = (segment_circuit, curr_map)
            pattern_list[p] = new_tuple
            for q in arch.nodes:
                if not q in segment_circuit.qubits:
                    segment_circuit.add_qubit(q)
            if p > 0:
                prev_map = pattern_list[p - 1][1]
                matching_dict = {}
                for k in curr_map["i"].keys():
                    matching_dict[prev_map["o"][k]] = curr_map["i"][k]
                swaps_as_pairs = get_token_swapping_network(arch, matching_dict)
                for pair in swaps_as_pairs:
                    new_c.SWAP(qubit_0=pair[0], qubit_1=pair[1])
            prev_bits = len(new_c.bits)
            b_map = {}
            for b in range(len(segment_circuit.bits)):
                b_map[segment_circuit.bits[b]] = Bit(prev_bits + b)
            segment_circuit.rename_units(b_map)
            new_c.add_circuit(segment_circuit, [])
        final_map = {"i": pattern_list[0][1]["i"], "o": pattern_list[-1][1]["o"]}
        return (new_c, final_map)

    @staticmethod
    def sequential(
        pattern_list: List[Tuple[Circuit, Dict[str, Dict[Qubit, int]]]],
        arch: Architecture,
    ) -> Tuple[Circuit, Dict[str, Dict[Qubit, int]]]:
        """
        Takes a list of tuples containing pytket circuits and their corresponding
        i/o maps as well as an architecture. It then routes the first circuit onto
        the architecture. For every subsequent circuit, the inputs are forcefully
        assigned to the location of the corresponding outputs of the previous segment
        on the architecture. The current subcircuit is then routed as best as possible
        given those initial conditions. One all subcircuits have been placed, the
        method returns a tuple containing the new circuit object and its corresponding
        i/o map.
        
        :param pattern_list:  List of tuples of circuits with their i/o maps.
        :param type:          List[Tuple[Circuit,Dict[str, Dict[Qubit, int]]]]
        
        :param arch:          The architecture onto which to fit the circuits.
        :param type:          Architecture
        
        :returns:        A tuple containing a circuit and an i/o map dictionary.
        :rtype:          Tuple[Circuit,Dict[str, Dict[Qubit, int]]]
        """
        new_c = Circuit()
        for q in arch.nodes:
            new_c.add_qubit(q)
        for p in range(len(pattern_list)):
            (curr_circ, curr_map) = pattern_list[p]
            if p > 0:
                (prev_circ, prev_map) = pattern_list[p - 1]
                route_map = {}
                for k in curr_map["i"].keys():
                    if curr_map["o"][k] == curr_map["i"][k]:
                        curr_map["o"][k] = prev_map["o"][k]
                    else:
                        curr_map["o"][k] = curr_circ.qubits[curr_map["o"][k]]
                    route_map[curr_circ.qubits[curr_map["i"][k]]] = prev_map["o"][k]
                    curr_map["i"][k] = prev_map["o"][k]
                curr_circ.rename_units(route_map)
            else:
                route_map = {}
                for k in curr_map["i"].keys():
                    if k in arch.nodes:
                        temp_val = curr_map["i"][k]
                        route_map[curr_circ.qubits[temp_val]] = k
                        curr_map["i"][k] = k
                        if curr_map["o"][k] == temp_val:
                            curr_map["o"][k] = k
                        else:
                            curr_map["o"][k] = curr_circ.qubits[curr_map["o"][k]]
                    else:
                        curr_map["i"][k] = curr_circ.qubits[curr_map["i"][k]]
                        curr_map["o"][k] = curr_circ.qubits[curr_map["o"][k]]
                curr_circ.rename_units(route_map)
            cu = CompilationUnit(curr_circ)
            if arch == FullyConnected:
                raise TypeError(
                    "This method doesn't work with a FullyConnected architecture."
                )
            RoutingPass(arch).apply(cu)
            used_nodes = set()
            all_nodes = set(arch.nodes)
            unassigned_qubits = []
            for k in curr_map["i"].keys():
                if cu.initial_map[curr_map["i"][k]] in all_nodes:
                    used_nodes |= {cu.initial_map[curr_map["i"][k]]}
                    used_nodes |= {cu.initial_map[curr_map["o"][k]]}
                    curr_map["i"][k] = cu.initial_map[curr_map["i"][k]]
                    curr_map["o"][k] = cu.final_map[curr_map["o"][k]]
                else:
                    unassigned_qubits.append(k)
            for command in cu.circuit.get_commands():
                if command.op.get_name() == "CZ":
                    qubits = command.qubits
                    for q in qubits:
                        if q in used_nodes:
                            used_nodes |= set(qubits)
                            break
                elif command.op.get_name() == "SWAP":
                    qubits = command.qubits
                    if (
                        (qubits[0] in used_nodes) and not (qubits[1] in used_nodes)
                    ) or (not (qubits[0] in used_nodes) and (qubits[1] in used_nodes)):
                        used_nodes ^= set(qubits)
                elif command.op.get_name() == "BRIDGE":
                    qubits = command.qubits
                    for q in {qubits[0], qubits[2]}:
                        if q in used_nodes:
                            used_nodes |= {qubits[0], qubits[2]}
                            break
            permutation = {x: x for x in all_nodes}
            for com in cu.circuit.commands_of_type(OpType.SWAP):
                permutation[com.qubits[0]] = com.qubits[1]
                permutation[com.qubits[1]] = com.qubits[0]
            unused_nodes = all_nodes - used_nodes
            segment_circuit = cu.circuit.copy()
            for uq in unassigned_qubits:
                temp = unused_nodes.pop()
                curr_map["o"][uq] = temp
                for k in permutation.keys():
                    if permutation[k] == temp:
                        segment_circuit.rename_units({curr_map["i"][uq]: k})
                        curr_map["i"][uq] = k
                        break
            new_tuple = (segment_circuit, curr_map)
            pattern_list[p] = new_tuple
            for q in arch.nodes:
                if not q in segment_circuit.qubits:
                    segment_circuit.add_qubit(q)
            prev_bits = len(new_c.bits)
            b_map = {}
            for b in range(len(segment_circuit.bits)):
                b_map[segment_circuit.bits[b]] = Bit(prev_bits + b)
            segment_circuit.rename_units(b_map)
            new_c.add_circuit(segment_circuit, [])
        final_map = {"i": pattern_list[0][1]["i"], "o": pattern_list[-1][1]["o"]}
        return (new_c, final_map)

    @staticmethod
    def unrouted(
        pattern_list: List[Tuple[Circuit, Dict[str, Dict[Qubit, int]]]],
        arch: Architecture = FullyConnected,
    ) -> Tuple[Circuit, Dict[str, Dict[Qubit, int]]]:
        """
        Takes a list of tuples containing pytket circuits and their corresponding
        i/o maps and joins them together into a new pytket circuit object assuming
        full connectivity. Then returns a tuple containing the new circuit object
        and its corresponding i/o map.
        
        :param pattern_list:  List of tuples of circuits with their i/o maps.
        :param type:          List[Tuple[Circuit,Dict[str, Dict[Qubit, int]]]]
        
        :returns:        A tuple containing a circuit and an i/o map dictionary.
        :rtype:          Tuple[Circuit,Dict[str, Dict[Qubit, int]]]
        """
        new_c = Circuit()
        prev_map: Dict[str, Dict[Qubit, Union[Qubit, int]]] = {}
        for p in range(len(pattern_list)):
            (curr_circ, curr_map) = pattern_list[p]
            if len(prev_map) == 0:
                new_c.add_circuit(curr_circ, [])
                prev_map = curr_map.copy()
            else:
                q_map = {}
                for k in prev_map["o"].keys():
                    q_map[curr_circ.qubits[curr_map["i"][k]]] = new_c.qubits[
                        prev_map["o"][k]
                    ]
                prev_ancillas = []
                curr_ancillas = []
                for q in new_c.qubits:
                    if not q in list(q_map.values()):
                        prev_ancillas.append(q)
                for q in curr_circ.qubits:
                    if not q in q_map.keys():
                        curr_ancillas.append(q)
                while len(prev_ancillas) > 0 and len(curr_ancillas) > 0:
                    q_map[curr_ancillas.pop()] = prev_ancillas.pop()
                if len(curr_ancillas) > 0:
                    unused_id_pool = []
                    for q in curr_circ.qubits:
                        if not q in list(q_map.values()):
                            unused_id_pool.append(q)
                    for q in curr_circ.qubits:
                        if not q in q_map.keys():
                            q_map[q] = unused_id_pool.pop()
                for k in curr_map["o"].keys():
                    curr_map["o"][k] = q_map[curr_circ.qubits[curr_map["o"][k]]]
                curr_circ.rename_units(q_map)
                prev_bits = len(new_c.bits)
                b_map = {}
                for b in range(len(curr_circ.bits)):
                    b_map[curr_circ.bits[b]] = Bit(prev_bits + b)
                curr_circ.rename_units(b_map)
                new_c.add_circuit(curr_circ, [])
                for k in curr_map["o"].keys():
                    curr_map["o"][k] = new_c.qubits.index(curr_map["o"][k])
                prev_map = curr_map.copy()
        final_map = {"i": pattern_list[0][1]["i"], "o": prev_map["o"]}
        for io in final_map.keys():
            for q in final_map[io].keys():
                final_map[io][q] = new_c.qubits[final_map[io][q]]
        return (new_c, final_map)
