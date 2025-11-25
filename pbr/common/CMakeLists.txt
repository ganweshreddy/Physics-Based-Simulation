set(LIBRARY_NAME pbr_common)

# Add your project files
file(GLOB SRCFILES src/*.cpp)
file(GLOB HFILES include/*.hpp)

add_library(${LIBRARY_NAME} ${SRCFILES} ${HFILES})
target_link_libraries(${LIBRARY_NAME} PUBLIC vislab_graphics)
target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(${LIBRARY_NAME} PROPERTIES FOLDER "pbr")
