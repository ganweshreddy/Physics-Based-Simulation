cmake_minimum_required(VERSION 3.18)

project(vislab
        LANGUAGES CXX
        DESCRIPTION "VisLab framework by Professorship for Visual Computing at FAU Erlangen-Nuernberg."
        HOMEPAGE_URL "https://vc.tf.fau.de/")

# ------------------------------------
# ------- Helper Functions -----------
# ------------------------------------
include(FetchContent)
include(cmake/enableIdeSupport.cmake)

# ------------------------------------
# ------------ Options ---------------
# ------------------------------------

option(VISLAB_BUILD_TESTS "Built vislab tests?" ON)

# The following options' default value depends on whether we are using a multi config generator such as Visual Studio or
# or XCode. We first have to load that global property into a local variable to be able to read it.
get_property(IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

if (IS_MULTI_CONFIG)
    option(VISLAB_IDE_SUPPORT
            "Enable support for IDEs such as Visual Studio or XCode (targets grouped into folders and source groups specified)." ON)
else()
    option(VISLAB_IDE_SUPPORT
            "Enable support for IDEs such as Visual Studio or XCode (targets grouped into folders and source groups specified)." OFF)
endif()

if (IS_MULTI_CONFIG OR CMAKE_BUILD_TYPE MATCHES Release)
    option(VISLAB_USE_OPENMP "Enable OpenMP support?" ON)
else()
    option(VISLAB_USE_OPENMP "Enable OpenMP support?" OFF)
endif()

# ------------------------------------
# ----- Compiler settings ------------
# ------------------------------------
# C++17
set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_definitions(-DNOMINMAX)
  add_compile_options(-bigobj)
  add_compile_options(-MP)
endif()

set(BUILD_SHARED_LIBS OFF)

# ------------------------------------
# ------- CMake Settings -------------
# ------------------------------------

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (VISLAB_IDE_SUPPORT)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()


# ------------------------------------
# ------------ Tests -----------------
# ------------------------------------
if(VISLAB_BUILD_TESTS)
    include(cmake/googletest.cmake)
    enable_testing()
endif()

# ------------------------------------
# ----- Compile vislab libraries -----
# ------------------------------------
add_subdirectory(vislab)

# ------------------------------------
# ----- Compile vislab apps ----------
# ------------------------------------
add_subdirectory(pbr)
add_subdirectory(physsim)