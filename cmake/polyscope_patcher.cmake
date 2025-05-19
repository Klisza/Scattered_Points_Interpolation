function(patch_poly TARGET_NAME NEW_NAME SOURCE_DIR)
    set(_cm "${SOURCE_DIR}/CMakeLists.txt")
    if(EXISTS "${_cm}")
        file(READ "${_cm}" _content)
        foreach(cmd IN ITEMS add_library target_link_libraries target_compile_definitions target_include_directories target_compile_options set_target_properties install)
            string(REGEX REPLACE "${cmd}\\([\t\r\n ]*${TARGET_NAME}" "${cmd}(${NEW_NAME}" _content "${_content}")
        endforeach()
        # Special case for glad
        string(REPLACE "get_target_property(library_type ${TARGET_NAME} TYPE)" "get_target_property(library_type ${NEW_NAME} TYPE)" _content "${_content}")
        string(REGEX REPLACE "add_library(\n ${TARGET_NAME})" "add_library(\n ${NEW_NAME})" _content "${_content}")
        file(WRITE "${_cm}" "${_content}")
        message(STATUS "Polyscope patched CMakeLists.txt for ${TARGET_NAME}")
    else()
        message(WARNING "CMakeLists.txt not found in ${SOURCE_DIR}!")
    endif()
endfunction()



