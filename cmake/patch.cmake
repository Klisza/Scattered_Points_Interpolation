function(patch_dep TARGET_NAME NEW_NAME SOURCE_DIR)
    set(_cm "${SOURCE_DIR}/CMakeLists.txt")
    if(EXISTS "${_cm}")
        file(READ "${_cm}" _content)
        foreach(cmd IN ITEMS add_library target_link_libraries target_compile_definitions target_include_directories target_compile_options set_target_properties install)
            string(REGEX REPLACE "${cmd}\\s*\\(\\s*${TARGET_NAME}" "${cmd}(${NEW_NAME}" _content "${_content}")
        endforeach()
        # Extra: also replace any occurrence of just the bare target name in install(TARGETS ...), etc.
        string(REGEX REPLACE "TARGETS\\s+${TARGET_NAME}" "TARGETS ${NEW_NAME}" _content "${_content}")
        file(WRITE "${_cm}" "${_content}")
        message(STATUS "Patched CMakeLists.txt for ${TARGET_NAME}")
    else()
        message(WARNING "CMakeLists.txt not found in ${SOURCE_DIR}!")
    endif()
endfunction()



