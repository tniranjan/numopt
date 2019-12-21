target_include_directories("${PROJECT_NAME}UnitTests" PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR})
FILE(GLOB_RECURSE SOURCES *.cpp *.h)
target_sources("${PROJECT_NAME}UnitTests" PUBLIC ${SOURCES})