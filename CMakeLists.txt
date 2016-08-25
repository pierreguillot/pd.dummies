cmake_minimum_required(VERSION 2.8)

project(pd-dummy)
include(pd.cmake)

include_directories(${PROJECT_SOURCE_DIR}/pd/src/)
if(${WIN32})
  link_directories(${PROJECT_SOURCE_DIR}/pdwinlib)
endif()

set_external_path("${PROJECT_SOURCE_DIR}/dummies/")
add_external(pak_project pak ${PROJECT_SOURCE_DIR}/src/pak.c)
add_external(crash_project crash ${PROJECT_SOURCE_DIR}/src/crash.c)
add_external(multisig_tilde_project multisig~ ${PROJECT_SOURCE_DIR}/src/multisig_tilde.c)
