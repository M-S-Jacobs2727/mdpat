-- Premake test

workspace "MDPAT"
    configurations { "debug", "release" }

project "MDPAT"
    architecture "x64"
    kind "ConsoleApp"
    language "C++"
    location "build"
    links { "mpi" }
    libdirs { os.findlib("mpi", "${HOME}/.local") }

    files { "src/*.hpp", "src/*.cpp" }

    includedirs { "${HOME}/.local/include" }

    filter "action:gmake2"
        buildoptions {"-std=c++17"}

    filter "configurations:debug"
        defines {"DEBUG"}
        symbols "On"

    filter "configurations:release"
        defines {"NDEBUG"}
        optimize "Speed"