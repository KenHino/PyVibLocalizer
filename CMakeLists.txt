cmake_minimum_required(VERSION 3.14) # CMakeのVersion指定
project(BlenderPythonVenv)
##### ここにBlender 関連のPATHを指定してください 
##### Modify PATH related to Blender (below is a example of MacOS)
set(BPython "/Applications/Blender.app/Contents/Resources/3.5/python/bin/python3.10")
set(BPythonDir "/Applications/Blender.app/Contents/Resources/3.5/python/")
set(Blender "/Applications/Blender.app/Contents/MacOS/Blender")
#####
#####
execute_process(COMMAND env PATH=${BPythonDir}/bin:$ENV{PATH} ${BPython} -m venv ${CMAKE_BINARY_DIR}/venv)
set(ENV{VIRTUAL_ENV} ${CMAKE_BINARY_DIR}/venv)
set(Python_FIND_VIRTUALENV FIRST)
set(BPython ${CMAKE_BINARY_DIR}/venv/bin/python)
execute_process(COMMAND ${BPython} -m pip install --upgrade pip)
execute_process(COMMAND ${BPython} -m pip install -r ${CMAKE_SOURCE_DIR}/requirements.txt)
execute_process(COMMAND blender_notebook install --blender-exec=${Blender})

file(APPEND ${CMAKE_BINARY_DIR}/venv/bin/activate
    "alias blender='"
    ${Blender}
    " --python-use-system-env'\n"
    "PYTHONPATH="
    ${CMAKE_SOURCE_DIR}
    "\n"
)

file(GLOB PYTHON_SITE_PACKAGES "${CMAKE_BINARY_DIR}/venv/lib/python*/site-packages/")
foreach(PYTHON_SITE_PACKAGE IN LISTS PYTHON_SITE_PACKAGES)
    file(WRITE ${PYTHON_SITE_PACKAGE}/pyvib.pth "${CMAKE_SOURCE_DIR}")
endforeach()


unset(Python_EXECUTABLE)
find_package(Python)
enable_testing()
foreach(name foo modules)
 add_test(test_${name} ${Python_EXECUTABLE} -m pytest
 ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_${name}.py)
 set_property(TEST test_${name} PROPERTY LABELS ${name})
endforeach(name)

