include_directories(${PROJECT_SOURCE_DIR}/src/libgwaspp ${PROJECT_SOURCE_DIR}/src ${Boost_INCLUDE_DIR})
link_directories(${PROJECT_BINARY_DIR}/lib ${R_LIBRARIES_DIR})

ADD_EXECUTABLE(Main main.cpp)
ADD_EXECUTABLE(GWAS gwas_basic.cpp)

TARGET_LINK_LIBRARIES(GWAS libgwaspp util gzstream ${GSL_LIBRARIES} ${R_RMATH} ${Boost_LIBRARIES} -lrt)
