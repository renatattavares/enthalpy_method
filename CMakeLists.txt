cmake_minimum_required(VERSION 3.21)

project(VCM2D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fontes
set(SOURCES
    VCM2D.cpp
    Functions/Model.cpp
    Functions/Plots.cpp
    Functions/Solver.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

# Includes
target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/Headers
        ${CMAKE_CURRENT_SOURCE_DIR}/Functions
        ${CMAKE_CURRENT_SOURCE_DIR}/Results
        ${CMAKE_CURRENT_SOURCE_DIR}/Validation
)

# yaml-cpp (via vcpkg)
find_package(yaml-cpp CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE yaml-cpp::yaml-cpp)

# Flags MSVC
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4 /permissive-)
endif()