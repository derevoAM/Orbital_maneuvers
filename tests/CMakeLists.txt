cmake_minimum_required(VERSION 3.24)
project(Orbital_maneuvers)

find_package(GTest REQUIRED)
file(GLOB files "*.cpp")

foreach (file ${files})
    get_filename_component(TName ${file} NAME)
    add_executable("${TName}" ${file})
    target_link_libraries(${TName}
            PRIVATE
            GTest::GTest)
endforeach ()

