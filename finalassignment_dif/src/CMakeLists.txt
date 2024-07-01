add_subdirectory(fhe)
add_subdirectory(circuits)
# add_compile_options(-march=core-avx2 -fomit-frame-pointer -std=c++17 -fopenmp)

INSTALL(TARGETS ${PROJECT_NAME}
       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
       LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
       ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY fhe
        DESTINATION include 
        FILES_MATCHING PATTERN "*.h")
install(DIRECTORY circuits
        DESTINATION include 
        FILES_MATCHING PATTERN "*.h")
