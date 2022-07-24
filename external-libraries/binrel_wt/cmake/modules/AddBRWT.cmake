function(brwt_enable_all_warnings TARGET_ID)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${TARGET_ID} PRIVATE
      -Weverything
      -Wno-c++98-compat
      -Wno-c++98-compat-pedantic
      -Wno-padded
    )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(${TARGET_ID} PRIVATE
      -Wall
      -Wextra
    )
  endif()
endfunction()

function(brwt_make_strict_modern_cpp TARGET_ID)
  set_property(TARGET ${TARGET_ID} PROPERTY CXX_EXTENSIONS       OFF)
  set_property(TARGET ${TARGET_ID} PROPERTY CXX_STANDARD          17)
  set_property(TARGET ${TARGET_ID} PROPERTY CXX_STANDARD_REQUIRED ON)
endfunction()

# Adds flags to the LINK_FLAGS property of the given target.
function(brwt_target_link_flags TARGET_ID)
  foreach(FLAG ${ARGN})
    set_property(TARGET ${TARGET_ID} APPEND_STRING
                 PROPERTY LINK_FLAGS " ${FLAG}")
  endforeach()
endfunction()

function(brwt_target_enable_lto TARGET_ID)
  set_property(TARGET ${TARGET_ID} PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)

  # The current CMake version does nothing when AppleClang is used, so...
  if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(LTO_FLAG "-flto") # todo: Try -flto=thin when available.
    brwt_target_link_flags(${TARGET_ID} ${LTO_FLAG})
    target_compile_options(${TARGET_ID} PRIVATE ${LTO_FLAG})
  endif()
endfunction()

function(brwt_configure TARGET_ID)
  brwt_enable_all_warnings(${TARGET_ID})
  brwt_make_strict_modern_cpp(${TARGET_ID})

  if (BRWT_ENABLE_LTO)
    brwt_target_enable_lto(${TARGET_ID})
  endif()
endfunction()

function(add_brwt_library TARGET_ID)
  add_library(${TARGET_ID} ${ARGN})
  brwt_configure(${TARGET_ID})
endfunction()

function(add_brwt_executable TARGET_ID)
  add_executable(${TARGET_ID} ${ARGN})
  brwt_configure(${TARGET_ID})
endfunction()
