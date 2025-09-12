# vim:ts=2:et
#=============================================================================#
#                        "CMakeLists.txt" (Top-Level):                        #
#=============================================================================#
#=============================================================================#
# CMake Settings:                                                             #
#=============================================================================#
# Export flags for convenience when using Vim's YCM plugin:
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.27)
CMAKE_POLICY(          VERSION 3.27)
CMAKE_POLICY(SET       CMP0167 NEW)

IF (NOT UNIX)
  MESSAGE(FATAL_ERROR "Non-UNIX platforms are not supported")
ENDIF()

#=============================================================================#
# Project Settings:                                                           #
#=============================================================================#
IF (NOT CMAKE_CXX_COMPILER OR NOT TOOL_CHAIN OR
    NOT CMAKE_BUILD_TYPE   OR NOT PROJ_NAME  OR
    NOT LIB_DIR            OR NOT BIN_DIR    OR NOT DEFINED UNCHECKED_MODE)
  MESSAGE(FATAL_ERROR "Required Definition(s) Missing")
ENDIF()
# NB: Also, CUDA_VER, CUDA_CC and CUDA_INCL may be passed from the CMake cmdl,
# but it is not an error if no CUDA is provided...

# NB: In this case, PROJECT_SOURCE_DIR = CMAKE_SOURCE_DIR
#              and  PROJECT_BINARY_DIR = CMAKE_BINARY_DIR ,
# so we use the "CMAKE_"-prefixed vars:
#
PROJECT(${PROJ_NAME} LANGUAGES CXX)

SET (CMAKE_EXPORT_COMPILE_COMMANDS ON)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${LIB_DIR}")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}")

#=============================================================================#
# Check for Dependencies:                                                     #
#=============================================================================#
# First of all, we need "DimType" (which is a Headers-Only library):
# 
SET(DIMTYPES_TOP_DIR  "${CMAKE_SOURCE_DIR}/../DimTypes")
IF (NOT IS_DIRECTORY  "${DIMTYPES_TOP_DIR}")
  MESSAGE(FATAL_ERROR "${DYMTYPES_TOP_DIR} not found")
ENDIF()

# Extend the PkgConfig search path (for *.pc files) -- useful for those modules
# which can be found via PkgConfig:
SET(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${ENV_PREFIX}/lib/pkgconfig")

# GSL:
FIND_PACKAGE(GSL REQUIRED)
SET(GSL_LIBS "gsl" "openblas")

# Boost:
# Boost (no specific version checks):
FIND_PACKAGE(Boost  REQUIRED)

# NLOpt:
FIND_PACKAGE(NLOPT  REQUIRED)

# XXX: For NOMAD, we currently need an explicit location:
SET(NOMAD_DIR  "/opt/NOMAD/4.5.0")
SET(NOMAD_LIBS "nomadAlgos" "nomadUtils" "nomadEval")

#=============================================================================#
# Compiler Settings:                                                          #
#=============================================================================#
MESSAGE(STATUS "Configuring for the ${TOOL_CHAIN} ToolChain...")

#-----------------------------------------------------------------------------#
IF("${TOOL_CHAIN}" STREQUAL "GCC")
#-----------------------------------------------------------------------------#
  # Optimisation and debugging flags:
  IF ("${CMAKE_BUILD_TYPE}"     STREQUAL "Release")
    # NB: Do NOT use:
    # "-ffast-math" option -- it will result in incorrect treatment of "NaN"s!
    # -flto"               -- it intereacts in a strange way with header-only
    #                         code!
    SET(CXXFLAGS "-O3 -mavx2 -funroll-loops -fomit-frame-pointer")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    SET(CXXFLAGS "-O2 -g -funroll-loops")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    SET(CXXFLAGS "-O0 -g")
  ELSE()
    MESSAGE(FATAL_ERROR "UnSupported BuildType: ${CMAKE_BUILD_TYPE}")
  ENDIF()

  SET(OPENMP_FLAGS "-fopenmp")
  SET(OPENMP_LIBS  "gomp")

  SET(CXXFLAGS "${CXXFLAGS} -march=native -mtune=native -std=gnu++23 \
                -Wall                     -Wextra            \
                -Wctor-dtor-privacy       -Wnon-virtual-dtor \
                -Woverloaded-virtual      -Wcast-qual        \
                -Winit-self               -Wwrite-strings    \
                -Wmissing-braces          -Wdisabled-optimization  \
                -Wconversion              -Wno-stringop-truncation \
                -Wno-class-memaccess      -Wno-double-promotion    \
                -Wno-c++20-compat         -Wno-ctor-dtor-privacy   \
                ${OPENMP_FLAGS}")

#-----------------------------------------------------------------------------#
ELSEIF("${TOOL_CHAIN}" STREQUAL "CLang")
#-----------------------------------------------------------------------------#
  # Optimisation and debugging flags:
  IF ("${CMAKE_BUILD_TYPE}"     STREQUAL "Release")
    SET(CXXFLAGS "-O3 -mavx2 -funroll-loops -fomit-frame-pointer")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    SET(CXXFLAGS "-O2 -g -funroll-loops")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"  )
    SET(CXXFLAGS "-O0 -g")
  ELSE()
    MESSAGE(FATAL_ERROR "UnSupported BuildType: ${CMAKE_BUILD_TYPE}")
  ENDIF()

  SET(OPENMP_FLAGS "-fopenmp=libomp")
  SET(OPENMP_LIBS  "omp")

  SET(CXXFLAGS "${CXXFLAGS}  -march=native    -mtune=native   -std=gnu++23  \
                -Weverything                  -Wno-unsafe-buffer-usage      \
                -Wno-c++98-compat             -Wno-c++98-compat-pedantic    \
                -Wno-weak-vtables             -Wno-double-promotion         \
                -Wno-weak-template-vtables    -Wno-covered-switch-default   \
                -Wno-switch-enum              -Wno-exit-time-destructors    \
                -Wno-global-constructors      -Wno-reserved-id-macro        \
                -Wno-padded                   -Wno-unused-macros            \
                -Wno-gnu-zero-variadic-macro-arguments                      \
                -Wno-disabled-macro-expansion                               \
                -Wno-float-equal              -Wno-missing-prototypes       \
                -Wno-vla                      -Wno-vla-extension            \
                -Wno-c99-extensions           -Wno-zero-length-array        \
                -Wno-packed                   -Wno-undefined-func-template  \
                -Wno-gnu-statement-expression -Wctor-dtor-privacy           \
                -Wno-c++20-compat             ${OPENMP_FLAGS}")

#-----------------------------------------------------------------------------#
ELSE()
  MESSAGE(FATAL_ERROR  "UnSupported ToolChain: ${TOOL_CHAIN}")
ENDIF()

#=============================================================================#
# Compiler-Independent Compilation Env:                                       #
#=============================================================================#
# Reset the standard flags for the selected mode, to prevent a mix-up:
SET(CMAKE_CXX_FLAGS_RELEASE        "")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")
SET(CMAKE_CXX_FLAGS_DEBUG          "")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS}")

ADD_DEFINITIONS(-DUNCHECKED_MODE=${UNCHECKED_MODE}
                -DBOOST_SYSTEM_NO_DEPRECATED
                -UBOOST_CB_ENABLE_DEBUG
                -DBOOST_CB_DISABLE_DEBUG)

IF (UNCHECKED_MODE)
  ADD_DEFINITIONS(-DNDEBUG)
ENDIF()

#-----------------------------------------------------------------------------#
# Common Configs:                                                             #
#-----------------------------------------------------------------------------#
# Use CCache if available:
OPTION(USE_CCACHE "Use ccache" OFF)
IF (USE_CCACHE)
  FIND_PROGRAM(CCACHE_FOUND ccache)
  IF (CCACHE_FOUND)
    SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_LINK    ccache)
  ENDIF()
ENDIF()

# During build, link executables with build rpath, not install rpath; they are
# then automatically re-linked on install. This is the default behaviour, any-
# way:
SET(CMAKE_SKIP_BUILD_RPATH         FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/Include"
                    "${DIMTYPES_TOP_DIR}/Include"
                     ${NLOPT_INCLUDE_DIR}
                    "${NOMAD_DIR}/include/nomad")

LINK_DIRECTORIES   ("${CMAKE_BINARY_DIR}"
                    "${LIB_DIR}"
                    "${NOMAD_DIR}/lib64")

#=============================================================================#
# The Target Library:                                                         #
#=============================================================================#
SET(SRCS
  Src/CoOrds/TimeScales.cpp
  Src/LVSC/Soyuz-2.1b/Stage1.cpp
  Src/LVSC/Soyuz-2.1b/Stage2.cpp
  Src/LVSC/Soyuz-2.1b/Stage3.cpp
  Src/PhysEffects/GravityFld-Earth.cpp
  Src/PhysEffects/GravityFld-Moon.cpp
  Src/PhysEffects/DE440T.cpp
  Src/PhysEffects/DE440T-Data.cpp
  Src/PhysEffects/EarthRotationModel.cpp
  Src/Missions/Ascent2.cpp)

# With CLang, we need some more .CPP files (to circumvent the unavailability of
# some "constexpr" functions):
IF ("${TOOL_CHAIN}" STREQUAL "CLang")
  LIST(APPEND SRCS
    Src/PhysEffects/EarthAtmosphereModel.cpp)
ENDIF()

ADD_LIBRARY          (${PROJECT_NAME} STATIC ${SRCS})

TARGET_LINK_LIBRARIES(${PROJECT_NAME}
                      ${NLOPT_LIBRARY}
                      ${NOMAD_LIBS}
                      ${OPENMP_LIBS})

#=============================================================================#
# Tests:                                                                      #
#=============================================================================#
SET(SB_TESTS
  ME_Test
  Soyuz21b_Stage1_Booster_Test
  Soyuz21b_Stage2_Test
  Soyuz21b_Stage3_Test
  LocationsTest
  AzimuthTest
  LagrangeNormTest
  LunarOrbiterTest
  ChebyshevTest
  TT_UTC_Test
  DE440T_Test0
  DE440T_Test1
  DE440T_Test2
  ERM_Test
  PrecNutTest
  TransLunarTest
  PlateAeroDynTest
  MkZhukovskyMesh
  MkSoyuz21bMesh
  LinAlgTuplesTest
  RKF5Test
  EarthAtmosphereTest
  AscentTest)

FOREACH(SB_TEST ${SB_TESTS})
  ADD_EXECUTABLE         (${SB_TEST} Tests/${SB_TEST}.cpp)
  TARGET_LINK_LIBRARIES  (${SB_TEST} ${PROJECT_NAME} ${GSL_LIBS})

  IF (${SB_TEST} STREQUAL "MkSoyuz21bMesh")
    TARGET_LINK_LIBRARIES(${SB_TEST} "gmsh")
  ENDIF()
ENDFOREACH(SB_TEST)
