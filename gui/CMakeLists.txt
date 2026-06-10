cmake_minimum_required(VERSION 3.16)
project(ChikaterasuGUI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Qt6 preferred; falls back to Qt5
find_package(Qt6 QUIET COMPONENTS Widgets)
if(Qt6_FOUND)
    set(QT_LIBS Qt6::Widgets)
else()
    find_package(Qt5 REQUIRED COMPONENTS Widgets)
    set(QT_LIBS Qt5::Widgets)
endif()

add_executable(ChikaterasuGUI
    main.cpp
    mainwindow.cpp
    mainwindow.h
)

target_link_libraries(ChikaterasuGUI PRIVATE ${QT_LIBS})

# macOS bundle
set_target_properties(ChikaterasuGUI PROPERTIES
    MACOSX_BUNDLE TRUE
    MACOSX_BUNDLE_GUI_IDENTIFIER "jp.kyoto-u.chikaterasu"
    MACOSX_BUNDLE_BUNDLE_NAME "ChikaterasuGUI"
)
