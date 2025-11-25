set(TEST_NAME ${LIBRARY_NAME}_test)

# Find source files
file(GLOB SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Create test executable
add_executable(${TEST_NAME} ${SOURCES})

# link
target_link_libraries(${TEST_NAME} PRIVATE vislab_analytic gtest gtest_main)

# add test
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})

if (VISLAB_IDE_SUPPORT)
	# Place in a folder
	set_target_properties(${TEST_NAME} PROPERTIES FOLDER "test")
endif()