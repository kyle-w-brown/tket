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

if(NOT TKET_SHARED) # True/False
    set(TKET_SHARED True)
endif()

if(NOT CONAN_EXPORTED)
    # not building with conan: generate conanbuildinfo.cmake
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
      message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
      file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
                    "${CMAKE_BINARY_DIR}/conan.cmake"
                    TLS_VERIFY ON)
    endif()
    include(${CMAKE_BINARY_DIR}/conan.cmake)

    conan_check(VERSION 1.52.0 REQUIRED)

    conan_add_remote(
        NAME tket-libs
        URL https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
        VERIFY_SSL True)

    # Make sure to use conanfile.py to define dependencies, to stay consistent
    conan_cmake_configure(
        REQUIRES
            boost/1.80.0
            symengine/0.9.0@tket/stable
            eigen/3.4.0
            nlohmann_json/3.11.2
            tklog/0.1.2@tket/stable
            tkassert/0.1.1@tket/stable
            tkrng/0.1.2@tket/stable
            tktokenswap/0.1.2@tket/stable
            tkwsm/0.2.1@tket/stable
        GENERATORS cmake_find_package)
    conan_cmake_autodetect(settings)
    conan_cmake_install(
        PATH_OR_REFERENCE ../../../recipes/tket
        BUILD missing
        SETTINGS ${settings}
        OPTIONS
            tklog:shared=True
            tket:shared=${TKET_SHARED})
endif()

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /EHsc")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wunreachable-code -Wunused")
ENDIF()
if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winconsistent-missing-override -Wloop-analysis")
endif()

set(PROFILE_COVERAGE no CACHE BOOL "Build library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (PROFILE_COVERAGE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
        # Bug in gcc 10: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95353
        IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-stringop-overflow")
        ENDIF()
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
