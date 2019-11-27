target_include_directories(numopt PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR})
FILE(GLOB_RECURSE SOURCES *.cpp *.h)
target_sources(numopt PUBLIC ${SOURCES})