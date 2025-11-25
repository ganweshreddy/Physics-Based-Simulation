# Get the name of the folder
get_filename_component(FolderName ${CMAKE_CURRENT_LIST_DIR} NAME)
string(REPLACE " " "_" FolderName ${FolderName})

set(EXECUTABLE_NAME pbr_${FolderName})

add_definitions(-DNOMINMAX)
add_executable(${EXECUTABLE_NAME} main.cpp)
target_link_libraries(${EXECUTABLE_NAME} PRIVATE 
	vislab_graphics 
	pbr_common
)

# Place in a folder
set_target_properties(${EXECUTABLE_NAME} PROPERTIES FOLDER "pbr")