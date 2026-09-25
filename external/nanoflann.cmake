if(NOT TARGET nanoflann)
  include(FetchContent)

  FetchContent_Declare(
    nanoflann
    GIT_REPOSITORY https://github.com/jlblancoc/nanoflann.git
    GIT_TAG v1.9.0
    GIT_SHALLOW TRUE)

  FetchContent_GetProperties(nanoflann)
  if(NOT nanoflann_POPULATED)
    FetchContent_Populate(nanoflann)

    add_library(nanoflann INTERFACE)
    add_library(nanoflann::nanoflann ALIAS nanoflann)
    target_include_directories(nanoflann
                               INTERFACE ${nanoflann_SOURCE_DIR}/include)
  endif()
endif()
