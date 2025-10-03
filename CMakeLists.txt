cmake_minimum_required (VERSION 3.20.1) 

project(HYBIRT VERSION 0.1 LANGUAGES CXX)

# Release mode default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()


set(HYBIRT_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_program(Git git)

function(hybirt_git_get_or_update name dir url branch)
  if (NOT EXISTS ${dir})
    message("cloning ${url} ${branch}  into ${dir}")
    execute_process(
      COMMAND ${Git} clone ${url} ${dir} -b ${branch} --recursive --depth 1 --shallow-submodules
    )
endif()
endfunction(hybirt_git_get_or_update)

function(hybirt_github_get_or_update name dir repo branch)
  hybirt_git_get_or_update(${name} ${dir} https://github.com/${repo} ${branch})
endfunction(hybirt_github_get_or_update)


if(NOT DEFINED HYBIRT_HIGHFIVE_VERSION)
    SET(HYBIRT_HIGHFIVE_VERSION "main")
endif()

set (HIGHFIVE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/highfive)

hybirt_github_get_or_update(HighFive ${HIGHFIVE_SRC} highfive-devs/highfive ${HYBIRT_HIGHFIVE_VERSION})

set(HIGHFIVE_UNIT_TESTS OFF) # silence warning
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_BUILD_DOCS OFF) # conflicts with phare doc target
set(HIGHFIVE_EXAMPLES OFF)

add_subdirectory(${HIGHFIVE_SRC})

message("HighFive enabled - checking HDF5")



set(SOURCE_INC
   src/ampere.hpp
   src/boundary_condition.hpp
   src/diagnostics.hpp
   src/faraday.hpp
   src/field.hpp
   src/gridlayout.hpp
   src/moments.hpp
   src/ohm.hpp
   src/particle.hpp
   src/population.hpp
   src/pusher.hpp
   src/utils.hpp
   src/vecfield.hpp
)

set (SOURCE_CPP
    src/hybirt.cpp
)

add_executable(hybirt ${SOURCE_INC} ${SOURCE_CPP})

target_link_libraries(hybirt PRIVATE HighFive)


add_subdirectory(tests/boris)

add_subdirectory(tests/ampere)

add_subdirectory(tests/faraday)


