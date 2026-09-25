include(FetchContent)

set(TINYUSDZ_BUILD_TESTS
    OFF
    CACHE BOOL "" FORCE)
set(TINYUSDZ_BUILD_EXAMPLES
    OFF
    CACHE BOOL "" FORCE)
set(TINYUSDZ_WITH_TOOL_USDA_PARSER
    OFF
    CACHE BOOL "" FORCE)
set(TINYUSDZ_WITH_TOOL_USDC_PARSER
    OFF
    CACHE BOOL "" FORCE)
set(TINYUSDZ_WITH_PYTHON
    OFF
    CACHE BOOL "" FORCE)
set(LIGHTUSD_BUILD_TESTS
    OFF
    CACHE BOOL "" FORCE)
set(LIGHTUSD_BUILD_EXAMPLES
    OFF
    CACHE BOOL "" FORCE)
set(LIGHTUSD_WITH_PYTHON
    OFF
    CACHE BOOL "" FORCE)
set(LIGHTUSD_NO_WERROR
    ON
    CACHE BOOL "" FORCE)

FetchContent_Declare(
  LightUSD
  GIT_REPOSITORY https://github.com/lighttransport/LightUSD.git
  GIT_TAG main
  GIT_SHALLOW TRUE)

FetchContent_MakeAvailable(LightUSD)

if(TARGET lightusd_static)
  target_include_directories(lightusd_static PUBLIC ${lightusd_SOURCE_DIR}/src)
endif()
