find_path(OPTIMIST_INCLUDE_DIR NAMES optimist/Optimist.h PATHS
  ${CMAKE_PREFIX_PATH}
  /usr/local/include
  /usr/include
)

find_library(OPTIMIST_LIBRARY NAMES Optimist PATHS
  ${CMAKE_PREFIX_PATH}
  /usr/local/lib
  /usr/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Optimist DEFAULT_MSG OPTIMIST_LIBRARY OPTIMIST_INCLUDE_DIR)

if(OPTIMIST_FOUND)
  add_library(Optimist::Optimist UNKNOWN IMPORTED)
  set_target_properties(Optimist::Optimist PROPERTIES
    IMPORTED_LOCATION "${OPTIMIST_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${OPTIMIST_INCLUDE_DIR}"
  )
endif()