cmake_minimum_required(VERSION 3.17)
project(HDD)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)

add_link_options(-static -static-libgcc -static-libstdc++)

set(SOURCE_FILES
        GeoJSONUtility.cpp
        Geometry.cpp
        Graph.cpp
        main.cpp
        Config.cpp)

add_executable(HDD
        ${SOURCE_FILES})
target_link_libraries(HDD PRIVATE nlohmann_json::nlohmann_json)
