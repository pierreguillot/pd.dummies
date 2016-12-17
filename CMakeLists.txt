cmake_minimum_required(VERSION 2.8)

set(CMAKE_SUPPRESS_REGENERATION true)
set(CMAKE_MACOSX_RPATH Off)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.4)
set(CMAKE_OSX_ARCHITECTURES "i386;x86_64")

include(pd.build/pd.cmake)

project(dummies)
set_pd_sources(${PROJECT_SOURCE_DIR}/pd/src/)
set_pd_external_path("${PROJECT_SOURCE_DIR}/dummies/")

add_pd_external(pak_project pak ${PROJECT_SOURCE_DIR}/src/pak.c)
add_pd_external(crash_project crash ${PROJECT_SOURCE_DIR}/src/crash.c)
add_pd_external(leeloo_tilde_project leeloo~ ${PROJECT_SOURCE_DIR}/src/leeloo_tilde.c)
add_pd_external(connected_tilde_project connected~ ${PROJECT_SOURCE_DIR}/src/connected_tilde.c)
add_pd_external(tabosco_tilde_project tabosco~ ${PROJECT_SOURCE_DIR}/src/tabosco_tilde.c)
add_pd_external(maxicatch_tilde_project maxicatch~ ${PROJECT_SOURCE_DIR}/src/maxicatch_tilde.c)
add_pd_external(meanblock_tilde_project meanblock~ ${PROJECT_SOURCE_DIR}/src/meanblock_tilde.c)
#add_pd_external(vbaf_project vbaf ${PROJECT_SOURCE_DIR}/src/vbaf.c)
