cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(gem)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 检查并找到 Garfield++
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# 查找 ROOT 库
find_package(ROOT REQUIRED)

# 包含 ROOT 的头文件
include_directories(${ROOT_INCLUDE_DIRS})

# ---Build executables----------------------------------------------------------
# 添加主程序
add_executable(gem gem.C)
add_executable(analysis analysis.c)
add_executable(output output.c)
add_executable(position position.c)
set_source_files_properties(analysis.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(output.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(position.c PROPERTIES LANGUAGE CXX)

# 链接 Garfield++ 和 ROOT 库
target_link_libraries(gem Garfield::Garfield ${ROOT_LIBRARIES})
target_link_libraries(analysis Garfield::Garfield ${ROOT_LIBRARIES})  # 新增analysis的链接
target_link_libraries(output Garfield::Garfield ${ROOT_LIBRARIES})
target_link_libraries(position Garfield::Garfield ${ROOT_LIBRARIES})

# ---处理数据文件---------------------------------------------------------------
# 创建构建目录中的 simdata 文件夹
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/simdata)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/result)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/comsol)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/fig)
# 定义源数据文件路径和目标路径
set(SIMDATA_FILES 
    ${CMAKE_SOURCE_DIR}/simdata/dielectrics.dat
    ${CMAKE_SOURCE_DIR}/simdata/THGEM.txt
    ${CMAKE_SOURCE_DIR}/simdata/THGEM.mphtxt
)

set(COMSOL_FILES 
    ${CMAKE_SOURCE_DIR}/comsol/THGEM.mph
)

# 复制数据文件到构建目录的 simdata 文件夹
add_custom_command(
    TARGET thgem POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SIMDATA_FILES}
        ${CMAKE_BINARY_DIR}/simdata/
    COMMENT "Copying simulation data files to build directory"
)
# 复制数据文件到构建目录的 comsol 文件夹
add_custom_command(
    TARGET thgem POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${COMSOL_FILES}
        ${CMAKE_BINARY_DIR}/comsol/
    COMMENT "Copying comsol files to build directory"
)

# 复制脚本文件到构建目录中
foreach(_file simulate.sh output.txt)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
