cmake_minimum_required(VERSION 2.8)

include(pd.build/pd.cmake)
project(pd-dummy)

set_pd_sources(${PROJECT_SOURCE_DIR}/pd/src/)
set_external_path("${PROJECT_SOURCE_DIR}/dummies/")

add_external(pak_project pak ${PROJECT_SOURCE_DIR}/src/pak.c)
add_external(crash_project crash ${PROJECT_SOURCE_DIR}/src/crash.c)
add_external(multisig_tilde_project multisig~ ${PROJECT_SOURCE_DIR}/src/multisig_tilde.c)
add_external(kgen_tilde_project kgen~ ${PROJECT_SOURCE_DIR}/src/kgen_tilde.c)
add_external(kop_tilde_project kop~ ${PROJECT_SOURCE_DIR}/src/kop_tilde.c)
