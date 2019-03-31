cmake_minimum_required(VERSION 3.13)

project(numopt)

add_executable(${PROJECT_NAME} "")

add_subdirectory(external)
add_subdirectory(src)
