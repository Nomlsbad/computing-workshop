cmake_minimum_required(VERSION 3.26)
project(computing_workshop LANGUAGES CXX)

#C++ compiler
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

#Eigen
FetchContent_Declare(
        eigen
        GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
        GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

#Task1
set(TASK_1_NAME ${PROJECT_NAME}_task1)

set (TASK_1_PATH "${PROJECT_SOURCE_DIR}/app/task1")
set (TASK_1_SRC_PATH  "${TASK_1_PATH}/src" )
set (TASK_1_INCLUDE_PATH  "${TASK_1_PATH}/include" )
set (TASK_1_MAIN "${TASK_1_SRC_PATH}/main.cpp")

file (GLOB_RECURSE TASK_1_SOURCES CONFIGURE_DEPENDS "${TASK_1_SRC_PATH}/*.cpp")
list(FILTER TASK_1_SOURCES EXCLUDE REGEX "${TASK_1_SRC_PATH}/*/main.cpp")
message(STATUS ${TASK_1_SOURCES} )

add_subdirectory(app/task1)