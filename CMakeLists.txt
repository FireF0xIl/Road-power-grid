cmake_minimum_required(VERSION 3.17)
project(Road-power-Grid)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)

add_link_options(-static -static-libgcc -static-libstdc++)

set(SOURCE_FILES
        GeoJSONHandler.cpp
        GeoJSONUtility.cpp
        Geometry.cpp
        Graph.cpp
        Map.cpp
        main.cpp
        GraphElement.h
        Config.cpp)

add_executable(Road-power-Grid
        ${SOURCE_FILES})
target_link_libraries(Road-power-Grid PRIVATE nlohmann_json::nlohmann_json)
