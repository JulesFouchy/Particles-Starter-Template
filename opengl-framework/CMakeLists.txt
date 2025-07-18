cmake_minimum_required(VERSION 3.20)

add_library(opengl_framework)
add_library(opengl_framework::opengl_framework ALIAS opengl_framework)
target_compile_features(opengl_framework PUBLIC cxx_std_20)

set(WARNINGS_AS_ERRORS_FOR_OPENGL_FRAMEWORK OFF CACHE BOOL "ON iff you want to treat warnings as errors for opengl_framework")

# Set warning level
if(MSVC)
    target_compile_options(opengl_framework PRIVATE /W4)
else()
    target_compile_options(opengl_framework PRIVATE -Wall -Wextra -Wpedantic -pedantic-errors -Wconversion -Wsign-conversion -Wimplicit-fallthrough)
endif()

# Maybe enable warnings as errors
if(WARNINGS_AS_ERRORS_FOR_OPENGL_FRAMEWORK)
    if(MSVC)
        target_compile_options(opengl_framework PRIVATE /WX)
    else()
        target_compile_options(opengl_framework PRIVATE -Werror)
    endif()
endif()

# ---Add glad---
add_library(glad lib/glad/src/gl.c)
target_include_directories(glad SYSTEM PUBLIC lib/glad/include)
target_link_libraries(opengl_framework PUBLIC glad)

# ---Add glm---
add_subdirectory(lib/glm)
target_link_libraries(opengl_framework PUBLIC glm)

# # ---Add img---
add_subdirectory(lib/img)
target_link_libraries(opengl_framework PUBLIC img::img)

# ---Add exe_path---
add_subdirectory(lib/exe_path)
target_link_libraries(opengl_framework PUBLIC exe_path::exe_path)

# ---Add tinyobjloader---
target_include_directories(opengl_framework PUBLIC lib/tinyobjloader)

# ---Add glfw---
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
add_subdirectory(lib/glfw)
target_link_libraries(opengl_framework PRIVATE glfw)

# ---Add Dear ImGui---
# add_library(ImGui)
# target_compile_features(ImGui PUBLIC cxx_std_11)

# if(WIN32 AND NOT MSVC)
# target_link_libraries(ImGui PRIVATE imm32) # Required only when using mingw, not Visual Studio
# endif()

# target_include_directories(ImGui SYSTEM PUBLIC lib/imgui)
# target_sources(ImGui PRIVATE
# lib/imgui/imgui_demo.cpp
# lib/imgui/imgui_draw.cpp
# lib/imgui/imgui_tables.cpp
# lib/imgui/imgui_widgets.cpp
# lib/imgui/imgui.cpp
# lib/imgui/misc/cpp/imgui_stdlib.cpp
# lib/imgui/backends/imgui_impl_glfw.cpp # Build ImGui with glfw
# lib/imgui/backends/imgui_impl_opengl3.cpp # Build ImGui with OpenGL
# )
# target_link_libraries(ImGui PRIVATE glfw)
# target_link_libraries(opengl_framework PUBLIC ImGui)

# ---Add source files---
file(GLOB_RECURSE OPENGL_FRAMEWORK_SOURCES CONFIGURE_DEPENDS src/*)
target_sources(opengl_framework PRIVATE ${OPENGL_FRAMEWORK_SOURCES})

if(WARNINGS_AS_ERRORS_FOR_OPENGL_FRAMEWORK)
    set(MAYBE_SYSTEM "")
else()
    set(MAYBE_SYSTEM "SYSTEM")
endif()

target_include_directories(opengl_framework ${MAYBE_SYSTEM} PUBLIC
    "include"

    # "lib" # For Dear ImGui
)

# ---Functions to copy files and folders---
include("CMakeUtils/files_and_folders.cmake")

function(gl_target_copy_file TARGET_NAME FILENAME)
    Cool__target_copy_file(${TARGET_NAME} ${FILENAME})
endfunction()

function(gl_target_copy_folder TARGET_NAME FOLDERNAME)
    Cool__target_copy_folder(${TARGET_NAME} ${FOLDERNAME})
endfunction()
