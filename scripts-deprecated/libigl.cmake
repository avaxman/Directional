if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/avaxman/libigl.git
    GIT_TAG v2.5.1
)
FetchContent_MakeAvailable(libigl)
