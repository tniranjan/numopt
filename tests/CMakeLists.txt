add_executable("${PROJECT_NAME}UnitTests" "")

target_compile_features("${PROJECT_NAME}UnitTests" PRIVATE cxx_std_17)
FILE(GLOB_RECURSE SOURCES *.cpp)
target_sources("${PROJECT_NAME}UnitTests" PUBLIC ${SOURCES})
target_link_libraries("${PROJECT_NAME}UnitTests" PUBLIC ${PROJECT_NAME} gtest gtest_main)
target_include_directories("${PROJECT_NAME}UnitTests" PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
add_test( "${PROJECT_NAME}UnitTests" "${PROJECT_NAME}UnitTests" )
