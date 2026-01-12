cmake_minimum_required(VERSION 3.20)
project(OPT)

set(CMAKE_CXX_STANDARD 17)

file(GLOB SRC CONFIGURE_DEPENDS "*.cpp" "*.h")

add_executable(${PROJECT_NAME} ${SRC})