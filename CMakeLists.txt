cmake_minimum_required(VERSION 3.24)
project(MicrofluidicDataCenter VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(MF_BUILD_VIEWER "Build the raylib 3D viewer" ON)
option(MF_BUILD_TESTS  "Build the unit tests"       ON)

# Some dependencies declare very old cmake_minimum_required() values that CMake 4 rejects.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "Minimum policy version for dependencies")

include(FetchContent)
set(FETCHCONTENT_QUIET ON)

# ---- Dependencies (pinned by version and SHA-256) --------------------------------------------
FetchContent_Declare(nlohmann_json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
  URL_HASH SHA256=d6c65aca6b1ed68e7a182f4757257b107ae403032760ed6ef121c9d55e81757d
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_MakeAvailable(nlohmann_json)

# ---- Headless physics core ------------------------------------------------------------------
add_library(mfcore STATIC
  src/core/fluid/WaterIAPWS.cpp
  src/core/hydraulics/RectChannel.cpp
  src/core/hydraulics/Edges.cpp
  src/core/hydraulics/HydraulicNetwork.cpp
  src/core/thermal/ChipGrid.cpp
  src/core/thermal/PowerMap.cpp
  src/core/thermal/AirGroup.cpp
  src/core/system/Server.cpp
  src/core/system/Rack.cpp
  src/core/system/CoolingPlant.cpp
  src/core/system/DataCenter.cpp
  src/core/workload/Workload.cpp
  src/core/config/Config.cpp
)
target_include_directories(mfcore PUBLIC src/core)
target_link_libraries(mfcore PUBLIC nlohmann_json::nlohmann_json)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(mfcore PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ---- Headless experiments CLI ---------------------------------------------------------------
add_executable(mfsim apps/mfsim/main.cpp)
target_link_libraries(mfsim PRIVATE mfcore)

# ---- 3D viewer --------------------------------------------------------------------------------
if(MF_BUILD_VIEWER)
  set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(BUILD_GAMES    OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(raylib
    URL https://github.com/raysan5/raylib/archive/refs/tags/5.5.tar.gz
    URL_HASH SHA256=aea98ecf5bc5c5e0b789a76de0083a21a70457050ea4cc2aec7566935f5e258e
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(raylib)

  # raygui is header-only; SOURCE_SUBDIR points at a folder with no CMakeLists so it is only downloaded.
  FetchContent_Declare(raygui
    URL https://github.com/raysan5/raygui/archive/refs/tags/4.0.tar.gz
    URL_HASH SHA256=299c8fcabda68309a60dc858741b76c32d7d0fc533cdc2539a55988cee236812
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR header-only)
  FetchContent_MakeAvailable(raygui)

  add_executable(mfsim_viewer
    apps/viewer/main.cpp
    apps/viewer/Scenes.cpp
    apps/viewer/UiPanel.cpp
    apps/viewer/raygui_impl.c
  )
  target_include_directories(mfsim_viewer PRIVATE apps/viewer ${raygui_SOURCE_DIR}/src ${raygui_SOURCE_DIR})
  target_link_libraries(mfsim_viewer PRIVATE mfcore raylib)
  if(APPLE)
    target_link_libraries(mfsim_viewer PRIVATE "-framework IOKit" "-framework Cocoa" "-framework OpenGL")
  endif()
endif()

# ---- Tests ------------------------------------------------------------------------------------
if(MF_BUILD_TESTS)
  FetchContent_Declare(doctest
    URL https://github.com/doctest/doctest/archive/refs/tags/v2.4.11.tar.gz
    URL_HASH SHA256=632ed2c05a7f53fa961381497bf8069093f0d6628c5f26286161fbd32a560186
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR header-only)
  FetchContent_MakeAvailable(doctest)

  enable_testing()
  add_executable(mfcore_tests
    tests/test_main.cpp
    tests/test_fluid.cpp
    tests/test_hydraulics.cpp
    tests/test_chip.cpp
    tests/test_system.cpp
    tests/test_plant.cpp
  )
  target_include_directories(mfcore_tests PRIVATE ${doctest_SOURCE_DIR})
  target_link_libraries(mfcore_tests PRIVATE mfcore)
  add_test(NAME mfcore_tests COMMAND mfcore_tests)
endif()
