# Copyright 2019-2023 Cambridge Quantum Computing
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

from pytket._tket.circuit_library import (
    _BRIDGE_using_CX_0,
    _BRIDGE_using_CX_1,
    _CX_using_TK2,
    _TK2_using_CX,
    _TK2_using_CX_and_swap,
    _approx_TK2_using_1xCX,
    _approx_TK2_using_2xCX,
    _TK2_using_3xCX,
    _CX_using_flipped_CX,
    _CX_using_ECR,
    _CX_using_ZZMax,
    _CX_using_ZZPhase,
    _CX_using_XXPhase_0,
    _CX_using_XXPhase_1,
    _CX_VS_CX_reduced,
    _CX_V_CX_reduced,
    _CX_S_CX_reduced,
    _CX_V_S_XC_reduced,
    _CX_S_V_XC_reduced,
    _CX_XC_reduced,
    _SWAP_using_CX_0,
    _SWAP_using_CX_1,
    _two_Rz1,
    _X1_CX,
    _Z0_CX,
    _CCX_modulo_phase_shift,
    _CCX_normal_decomp,
    _C3X_normal_decomp,
    _C4X_normal_decomp,
    _ladder_down,
    _ladder_down_2,
    _ladder_up,
    _X,
    _CX,
    _CCX,
    _BRIDGE,
    _H_CZ_H,
    _CZ_using_CX,
    _CY_using_CX,
    _CH_using_CX,
    _CV_using_CX,
    _CVdg_using_CX,
    _CSX_using_CX,
    _CSXdg_using_CX,
    _CSWAP_using_CX,
    _ECR_using_CX,
    _ZZMax_using_CX,
    _CRz_using_TK2,
    _CRz_using_CX,
    _CRx_using_TK2,
    _CRx_using_CX,
    _CRy_using_TK2,
    _CRy_using_CX,
    _CU1_using_TK2,
    _CU1_using_CX,
    _CU3_using_CX,
    _ISWAP_using_TK2,
    _ISWAP_using_CX,
    _XXPhase_using_TK2,
    _XXPhase_using_CX,
    _YYPhase_using_TK2,
    _YYPhase_using_CX,
    _ZZPhase_using_TK2,
    _ZZPhase_using_CX,
    _TK2_using_ZZPhase,
    _TK2_using_ZZPhase_and_swap,
    _TK2_using_TK2_or_swap,
    _approx_TK2_using_1xZZPhase,
    _approx_TK2_using_2xZZPhase,
    _TK2_using_ZZMax,
    _TK2_using_ZZMax_and_swap,
    _XXPhase3_using_TK2,
    _XXPhase3_using_CX,
    _ESWAP_using_TK2,
    _ESWAP_using_CX,
    _FSim_using_TK2,
    _FSim_using_CX,
    _PhasedISWAP_using_TK2,
    _PhasedISWAP_using_CX,
    _NPhasedX_using_PhasedX,
    _TK2_using_normalised_TK2,
    _TK1_to_PhasedXRz,
    _TK1_to_RzRx,
    _TK1_to_RzH,
    _TK1_to_RzSX,
    _TK1_to_TK1,
)