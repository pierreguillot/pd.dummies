cmake_minimum_required(VERSION 2.8)

set(CMAKE_SUPPRESS_REGENERATION true)
set(CMAKE_MACOSX_RPATH Off)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.4)
set(CMAKE_OSX_ARCHITECTURES "ppc;i386;x86_64")

include(pd.build/pd.cmake)

project(pd-dummy)
set_pd_sources(${PROJECT_SOURCE_DIR}/pd/src/)
set_pd_external_path("${PROJECT_SOURCE_DIR}/dummies/")

add_pd_external(pak_project pak ${PROJECT_SOURCE_DIR}/src/pak.c)
add_pd_external(crash_project crash ${PROJECT_SOURCE_DIR}/src/crash.c)
add_pd_external(multisig_tilde_project multisig~ ${PROJECT_SOURCE_DIR}/src/multisig_tilde.c)
add_pd_external(kgen_tilde_project kgen~ ${PROJECT_SOURCE_DIR}/src/kgen_tilde.c)
add_pd_external(kop_tilde_project kop~ ${PROJECT_SOURCE_DIR}/src/kop_tilde.c)
