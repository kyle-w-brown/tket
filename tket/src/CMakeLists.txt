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

cmake_minimum_required(VERSION 3.24.1)

project(tket LANGUAGES CXX)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE) # Release/Debug
    set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT TKLOG_SHARED) # True/False
    set(TKLOG_SHARED False)
endif()

if(NOT TKET_SHARED) # True/False
    set(TKET_SHARED False)
endif()

if(NOT TKET_COVERAGE) # True/False
    set(TKET_COVERAGE False)
endif()

if(CONAN_EXPORTED)
    if(TKET_WITH_TESTS OR TKET_WITH_PROPTESTS OR TKET_COVERAGE)
        message(FATAL_ERROR
            "When building with conan, these variables should be unset (or set
            to False): TKET_WITH_TESTS; TKET_WITH_PROPTESTS; TKET_COVERAGE."
            )
    endif()
endif()

if(NOT CONAN_EXPORTED)
    # not building with conan: generate conanbuildinfo.cmake

    include(${CMAKE_CURRENT_SOURCE_DIR}/../../conan.cmake)

    conan_check(VERSION 1.52.0 REQUIRED)

    conan_add_remote(
        NAME tket-libs
        URL https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
        VERIFY_SSL True)
    conan_cmake_autodetect(settings)
    conan_cmake_install(
        PATH_OR_REFERENCE ${CMAKE_CURRENT_SOURCE_DIR}/../../recipes/tket
        SETTINGS ${settings}
        BUILD missing
        OPTIONS
            tklog:shared=${TKLOG_SHARED}
            tket:shared=${TKET_SHARED})

    if(TKET_WITH_TESTS)
        conan_cmake_install(
            PATH_OR_REFERENCE catch2/3.1.0@
            SETTINGS ${settings}
            BUILD missing
            GENERATOR deploy)
    endif()

    if(TKET_WITH_PROPTESTS)
        conan_cmake_install(
            PATH_OR_REFERENCE rapidcheck/cci.20220514@
            SETTINGS ${settings}
            BUILD missing
            GENERATOR deploy)
    endif()
endif()

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

file(WRITE ${CMAKE_BINARY_DIR}/tk_ver.txt ${CONAN_PACKAGE_VERSION})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /EHsc /MD")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wunreachable-code -Wunused")
ENDIF()
if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winconsistent-missing-override -Wloop-analysis")
endif()

IF (TKET_COVERAGE)
    IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
        # Bug in gcc 10: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95353
        IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-stringop-overflow")
        ENDIF()
    ELSE()
        message(FATAL_ERROR "TKET_COVERAGE only supported with GNU compiler.")
    ENDIF()
ENDIF()

IF (WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
ELSEIF(APPLE)
    # set correct install_name
    set(CMAKE_INSTALL_NAME_DIR "@loader_path")
    set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
ELSEIF(UNIX)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH "\${ORIGIN}")
ENDIF()

if(TKET_SHARED STREQUAL "True")
    set(BUILD_SHARED_LIBS "TRUE")
endif()

# if you add new modules here make sure that it is added at the right position
# this list corresponds to a topological sorting of the dependency graph of the different modules
list(APPEND TKET_COMPS
    Utils
    OpType
    Clifford
    Ops
    Graphs
    Gate
    PauliGraph
    Circuit
    Architecture
    Simulation
    Diagonalisation
    Characterisation
    ZX
    Converters
    Placement
    ArchAwareSynth
    Mapping
    MeasurementSetup
    Transformations
    Predicates
    )

foreach(COMP ${TKET_COMPS})
    set(TKET_${COMP}_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${COMP}/include)
    add_subdirectory(${COMP})
    file(GLOB TKET_${COMP}_PUBLIC_HEADERS ${COMP}/include/${COMP}/*.hpp)
    set_target_properties(tket-${COMP} PROPERTIES PUBLIC_HEADER "${TKET_${COMP}_PUBLIC_HEADERS}")
    install(TARGETS tket-${COMP}
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/${COMP})
endforeach()

if(TKLOG_SHARED)
    find_library(TKLOG_LIB NAMES tklog PATHS CONAN_PKG_LIBS_TKLOG)
    install(FILES ${TKLOG_LIB} DESTINATION lib)
endif()

if(TKET_SHARED)
    list(APPEND TKET_DEPS
        ${CONAN_LIBS_SYMENGINE}
        ${CONAN_LIBS_TKLOG}
        ${CONAN_LIBS_TKASSERT}
        ${CONAN_LIBS_TKRNG}
        ${CONAN_LIBS_TKTOKENSWAP}
        ${CONAN_LIBS_TKWSM}
        )
endif()

if(TKET_WITH_TESTS)
    set(TKET_TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../tests)

    include(${TKET_TESTS_DIR}/tkettestutilsfiles.cmake)
    include(${TKET_TESTS_DIR}/tkettestsfiles.cmake)

    add_executable(test_tket ${TESTUTILS_SOURCES} ${TEST_SOURCES})

    foreach(COMP ${TKET_COMPS})
        target_include_directories(test_tket PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${COMP}/include)
        target_link_libraries(test_tket PRIVATE tket-${COMP})
    endforeach()

    if(TKET_SHARED)
        target_link_libraries(test_tket PRIVATE ${TKET_DEPS})
    endif()

    target_include_directories(test_tket PRIVATE ${CMAKE_BINARY_DIR}/catch2/include)
    target_link_directories(test_tket PRIVATE ${CMAKE_BINARY_DIR}/catch2/lib)
    if(${CMAKE_BUILD_TYPE} MATCHES "Release")
        set(CATCH2MAIN_LIB "Catch2Main")
        set(CATCH2_LIB "Catch2")
    else()
        set(CATCH2MAIN_LIB "Catch2Maind")
        set(CATCH2_LIB "Catch2d")
    endif()
    target_link_libraries(test_tket
        PRIVATE ${CATCH2MAIN_LIB}
        PRIVATE ${CATCH2_LIB})
    install(TARGETS test_tket DESTINATION bin)
    install(FILES ${TKET_TESTS_DIR}/test_circuits/lexiroute_circuit.json DESTINATION bin)
endif()

if(TKET_WITH_PROPTESTS)
    set(TKET_PROPTESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../proptests)
    add_executable(proptest ${TKET_PROPTESTS_DIR}/ComparisonFunctions.cpp ${TKET_PROPTESTS_DIR}/proptest.cpp)
    foreach(COMP ${TKET_COMPS})
        target_include_directories(proptest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${COMP}/include)
        target_link_libraries(proptest PRIVATE tket-${COMP})
    endforeach()

    if(TKET_SHARED)
        target_link_libraries(proptest PRIVATE ${TKET_DEPS})
    endif()

    target_include_directories(proptest PRIVATE ${CMAKE_BINARY_DIR}/rapidcheck/include)
    target_link_directories(proptest PRIVATE ${CMAKE_BINARY_DIR}/rapidcheck/lib)
    target_link_libraries(proptest PRIVATE rapidcheck)
    install(TARGETS proptest DESTINATION bin)
endif()
