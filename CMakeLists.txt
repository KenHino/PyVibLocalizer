cmake_minimum_required(VERSION 3.14) # CMakeのVersion指定
project(pyvib NONE) # fibonacciというプロジェクト NONEはC++でないということ
find_package(Python) # 実行可能なPython packageを探してくる

execute_process(COMMAND ${Python_EXECUTABLE} -m venv
 ${CMAKE_BINARY_DIR}/venv)
set(ENV{VIRTUAL_ENV} ${CMAKE_BINARY_DIR}/venv)
set(Python_FIND_VIRTUALENV FIRST)
unset(Python_EXECUTABLE)
find_package(Python)
execute_process(COMMAND ${Python_EXECUTABLE} -m pip install
 --upgrade pip)
execute_process(COMMAND ${Python_EXECUTABLE} -m pip install -r
 ${CMAKE_SOURCE_DIR}/requirements.txt)
enable_testing()

foreach(name foo modules)
 add_test(test_${name} ${Python_EXECUTABLE} -m pytest
 ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_${name}.py)
endforeach(name)

