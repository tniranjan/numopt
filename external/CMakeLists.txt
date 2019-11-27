SET(EIGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen-git-mirror)
target_include_directories(numopt PUBLIC ${EIGEN_INCLUDE_DIR})