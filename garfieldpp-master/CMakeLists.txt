# Build Garfield++ library and examples 
# Authors: Klaus Zenker, Pere Mato, Stefano Caiazza 
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# Setting the default build type to be Release.
# Otherwise the default mode is Debug.
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "${default_build_type}"
      CACHE STRING "Choose the type of build.")
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                               "RelWithDebInfo")
endif()

if (USEGPU)
  ADD_DEFINITIONS(-DUSEGPU)
  project(Garfield
    LANGUAGES CXX Fortran CUDA
    VERSION 0.3.0)
  set(CMAKE_CUDA_ARCHITECTURES native)
else()
  project(
    Garfield
    LANGUAGES CXX Fortran
    VERSION 0.3.0)
endif()

#--- Set default installation prefix ------------------------------------------
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(DEFINED ENV{GARFIELD_INSTALL} )
    set(CMAKE_INSTALL_PREFIX 
        $ENV{GARFIELD_INSTALL}
        CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
  elseif(DEFINED ENV{GARFIELD_HOME} )
    set(CMAKE_INSTALL_PREFIX
        $ENV{GARFIELD_HOME}/install
        CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
  else()
    set(CMAKE_INSTALL_PREFIX
        ${PROJECT_SOURCE_DIR}/install
        CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
  endif()
else()
  message(
    STATUS "The installation folder is set to its default value ${CMAKE_INSTALL_PREFIX}")
endif()
message(STATUS "Install path: ${CMAKE_INSTALL_PREFIX}")
message(
  STATUS
    "If you want to change this path call cmake -DCMAKE_INSTALL_PREFIX=my_install_path ..."
)

# Adding the CMake folder to the search path for modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)

#--- Use GNU-style hierarchy for installing build products --------------------
include(GNUInstallDirs)
include(BuildUtils)
# force_color_output()

#--- Dependencies -------------------------------------------------------------
if(DEFINED ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
elseif(DEFINED ROOTSYS)
  list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})
elseif(DEFINED ROOT_CMAKE_DIR)
  list(APPEND CMAKE_PREFIX_PATH ${ROOT_CMAKE_DIR})
# else()
#   message(
#     STATUS
#       "Please consider to set ROOTSYS or use -DROOTSYS=..."
#       "If still there is a problem, point to the directory which includes FindROOT.cmake using -DROOT_CMAKE_DIR=..."
# )
endif(DEFINED ENV{ROOTSYS})

find_package(ROOT 6.0 REQUIRED COMPONENTS Geom Gdml MathCore)
#--- Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY) --
if(${ROOT_VERSION} VERSION_LESS "6.20.0")
  # It only makes sense for ROOT < 6.20. 
  # Afterwards it is loaded automatically by the find_package
  get_filename_component(Root_Use_File_Dir ${ROOT_USE_FILE} DIRECTORY)
  include("${Root_Use_File_Dir}/modules/RootNewMacros.cmake")
endif()

#--- Set the C++ standard (depending on the ROOT version) ---------------------
if(NOT DEFINED CMAKE_CXX_STANDARD)
  if(DEFINED ROOT_CXX_STANDARD)
    set(GARFIELD_CXX_STANDARD ${ROOT_CXX_STANDARD})
  elseif(${ROOT_CXX_FLAGS} MATCHES ".*-std=c\\+\\+([0-9][0-9])")
    set(GARFIELD_CXX_STANDARD ${CMAKE_MATCH_1})
  elseif(${ROOT_CXX_FLAGS} MATCHES ".*-std=c\\+\\+1z")
    set(GARFIELD_CXX_STANDARD 17)
  endif()
else()
  set(GARFIELD_CXX_STANDARD ${CMAKE_CXX_STANDARD})
endif()

find_package(GSL REQUIRED)
find_package(OpenMP)

# Add Magboltz
add_subdirectory(Magboltz)

# Add Degrade
option(WITH_DEGRADE "Build Degrade." ON)
if(WITH_DEGRADE)
  add_subdirectory(Degrade)
endif()

# Create GarfieldRandom library
add_library(GarfieldRandom SHARED
  Source/Random.cc 
  Source/RandomEngineRoot.cc 
  Source/RandomEngineSTL.cc
)
target_include_directories(GarfieldRandom PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>)
target_link_libraries(GarfieldRandom PRIVATE ROOT::MathCore)
set_property(TARGET GarfieldRandom PROPERTY CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(TARGET GarfieldRandom PROPERTY CXX_EXTENSIONS OFF)
set_property(TARGET GarfieldRandom PROPERTY CXX_STANDARD ${GARFIELD_CXX_STANDARD})

# Creating the main Garfield library
add_library(Garfield SHARED "")
target_link_libraries(Garfield PUBLIC GarfieldRandom)
set_property(TARGET Garfield PROPERTY CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(TARGET Garfield PROPERTY CXX_EXTENSIONS OFF)
set_property(TARGET Garfield PROPERTY CXX_STANDARD ${GARFIELD_CXX_STANDARD})
add_library(Garfield::Garfield ALIAS Garfield)
target_sources(
  Garfield
  PRIVATE Source/AvalancheGrid.cc
          Source/AvalancheGridSpaceCharge.cc
          Source/AvalancheMC.cc
          Source/AvalancheMicroscopic.cc
          Source/Component.cc
          Source/ComponentAnalyticField.cc
          Source/ComponentAnsys121.cc
          Source/ComponentAnsys123.cc
          Source/ComponentCST.cc
          Source/ComponentComsol.cc
          Source/ComponentConstant.cc
          Source/ComponentElmer.cc
          Source/ComponentElmer2d.cc
          Source/ComponentFieldMap.cc
          Source/ComponentGrid.cc
          Source/ComponentNeBem2d.cc
          Source/ComponentNeBem3d.cc
          Source/ComponentNeBem3dMap.cc
          Source/ComponentParallelPlate.cc
          Source/ComponentTcad2d.cc
          Source/ComponentTcad3d.cc
          Source/ComponentTcadBase.cc
          Source/ComponentUser.cc
          Source/ComponentVoxel.cc
          Source/DriftLineRKF.cc
          Source/GeometryRoot.cc
          Source/GeometrySimple.cc
          Source/KDTree.cc
          Source/Medium.cc
          Source/MediumCdTe.cc
          Source/MediumDiamond.cc
          Source/MediumGaAs.cc
          Source/MediumGaN.cc
          Source/MediumGas.cc
          Source/MediumMagboltz.cc
          Source/MediumSilicon.cc
          Source/NeBemInterface.cpp
          Source/Numerics.cc
          Source/OpticalData.cc
          Source/PlottingEngine.cc
          Source/Polygon.cc
          Source/QuadTree.cc
          Source/Sensor.cc
          Source/Shaper.cc
          Source/Solid.cc
          Source/SolidBox.cc
          Source/SolidExtrusion.cc
          Source/SolidHole.cc
          Source/SolidRidge.cc
          Source/SolidSphere.cc
          Source/SolidTube.cc
          Source/SolidWire.cc
          Source/TetrahedralTree.cc
          Source/TGeoTet.cc
          Source/Track.cc
          Source/TrackBichsel.cc
          Source/TrackElectron.cc
          Source/TrackHeed.cc
          Source/TrackPAI.cc
          Source/TrackSimple.cc
          Source/TrackSrim.cc
          Source/TrackTrim.cc
          Source/ViewBase.cc
          Source/ViewCell.cc
          Source/ViewDrift.cc
          Source/ViewFEMesh.cc
          Source/ViewField.cc
          Source/ViewGeometry.cc
          Source/ViewIsochrons.cc
          Source/ViewMedium.cc
          Source/ViewSignal.cc)
if(WITH_DEGRADE)
  target_sources(
    Garfield
    PRIVATE Source/TrackDegrade.cc
  )
endif()

#--- Add GPU support if requested ---------------------------------------------
if (USEGPU)
  set(gpu_sources
    GPU/AvalancheMicroscopicGPU.cu
    GPU/ComponentGPU.cu
    GPU/GPUFunctions.cu
    GPU/MediumGPU.cu
    GPU/RandomEngineGPU.cu
    GPU/SensorGPU.cu
    GPU/TetrahedralTreeGPU.cu)

  # for cuda builds, need to remove -pedantic and -Wshadow to avoid *a lot* of errors/warnings
  set_source_files_properties(${gpu_sources} PROPERTIES POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 17
    CMAKE_CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    COMPILE_OPTIONS "-Wall;-Wextra;-Wfatal-errors;-Wno-missing-braces"
    INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/Source)

  target_sources( Garfield PRIVATE ${gpu_sources})
  set_property(TARGET Garfield PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

#--- Create ROOT dictionary ---------------------------------------------------
if(ROOT_VERSION VERSION_LESS 6.17)
  file(
    GLOB garfield_header
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/Include
    Include/Garfield/*.hh)
  root_generate_dictionary(
    GarfieldDict
    ${garfield_header}
    LINKDEF
    Include/Garfield/LinkDef.h
    OPTIONS
    -I${CMAKE_CURRENT_SOURCE_DIR}/Include
    -I${CMAKE_CURRENT_SOURCE_DIR}/GPU
    -I${PROJECT_SOURCE_DIR}/Heed)
else()
  file(
    GLOB garfield_header
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/Include
    Include/Garfield/*.hh)
  include_directories(Include Heed GPU)
  root_generate_dictionary(GarfieldDict ${garfield_header} LINKDEF
                           Include/Garfield/LinkDef.h)
endif()
target_sources(Garfield PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/GarfieldDict.cxx)

#--- Install pcm file in case of ROOT 6 ---------------------------------------
if(ROOT_VERSION VERSION_GREATER 6.06)
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/libGarfieldDict_rdict.pcm
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT dict)
endif()

#--- If you cannot create symbolic links turn FALSE the library version -------
option(BUILD_WITH_VERSION_NUMBER
       "If you cannot create symbolic links turn FALSE the library version" TRUE)
if(BUILD_WITH_VERSION_NUMBER)
  set_target_properties(Garfield PROPERTIES SOVERSION ${Garfield_VERSION})
endif()
mark_as_advanced(BUILD_WITH_VERSION_NUMBER)

# Add neBEM 
# TODO Decouple Nebem interface from the Garfield components
# The Nebem interface is somewhat coupled with Garfield itself.
# We should decouple them to keep them independent, but until then 
# they have to be compiled in the same library to make it work
set(nebem_sources
  NeBem/ComputeProperties.c
  NeBem/Isles.c
  NeBem/luc.c
  NeBem/neBEM.c
  NeBem/neBEMInterface.c
  NeBem/nrutil.c
  NeBem/ReTriM.c
  NeBem/svdcmp.c
  NeBem/Vector.c)

# Until CMake 3.18 the source file properties are only visible to targets 
# added in the same directory therefore we cannot add Nebem through the 
# add_subdirectory mechanism unless we change the cpp file extension so that
# CMake can recognize them automatically as CXX and not C
set_source_files_properties(${nebem_sources} PROPERTIES LANGUAGE CXX)
target_sources(Garfield PRIVATE ${nebem_sources})

# Add Heed
add_subdirectory(Heed)

target_link_libraries(Garfield PUBLIC ROOT::Geom
                                      ROOT::Gdml
                                      ROOT::Graf3d 
                                      GSL::gsl
                                      Garfield::Heed
                                      Magboltz::Magboltz
                                      $<$<BOOL:${WITH_DEGRADE}>:Degrade::Degrade>
                                      $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
                                      )

target_include_directories(
  Garfield
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/NeBem>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GPU>)

#--- Install the header files -------------------------------------------------
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  REGEX LinkDef EXCLUDE)

install(
  TARGETS Garfield GarfieldRandom
  EXPORT "${PROJECT_NAME}Targets"
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Data
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Garfield/)

#--- Target for Doxygen documentation -----------------------------------------
option(WITH_DOC "Whether or not to create doxygen doc target." OFF)
if(WITH_DOC)
  add_subdirectory(Doc)
endif()

#--- Build the examples -------------------------------------------------------
option(WITH_EXAMPLES "Build Garfield++ examples" ON)
if(WITH_EXAMPLES)
  message(
    STATUS
      "Garfield++ examples will be built. In order disable this option use -DWITH_EXAMPLES=OFF"
  )
  add_subdirectory(Examples)
else()
  message(STATUS "Building of examples is switched off.")
endif()

add_subdirectory(CMake)

option(WITH_TESTS "Build Garfield++ tests" ON)
if(WITH_TESTS)
  include(CTest)
  add_subdirectory(Tests)
endif()
