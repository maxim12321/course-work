QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS
# Tell the qcustomplot header that it will be used as library:
DEFINES += QCUSTOMPLOT_USE_LIBRARY

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Link with debug version of qcustomplot if compiling in debug mode, else with release library:
CONFIG(debug, release|debug) {
  win32:QCPLIB = qcustomplotd2
  else: QCPLIB = qcustomplotd
} else {
  win32:QCPLIB = qcustomplot2
  else: QCPLIB = qcustomplot
}
LIBS += -L../course-work/qcustomplot/ -l$$QCPLIB

SOURCES += \
    grid_data_loader.cpp \
    grid_data_processor.cpp \
    heatmap.cpp \
    main.cpp \
    main_window.cpp \
    menu.cpp \
    properties_item.cpp \
    properties_widget.cpp \
    solution_runner.cpp \
    solver/column_iteration_solver.cpp \
    solver/properties_manager.cpp \
    solver/row_iteration_solver.cpp \
    solver/solver.cpp \
    solver/tridiagonal.cpp \
    solver/utils/matrix.cpp \
    solver/utils/vector.cpp

HEADERS += \
    grid_data_loader.h \
    grid_data_processor.h \
    heatmap.h \
    main_window.h \
    menu.h \
    properties_item.h \
    properties_widget.h \
    solution_runner.h \
    solver/column_iteration_solver.h \
    solver/properties_manager.h \
    solver/row_iteration_solver.h \
    solver/solver.h \
    solver/tridiagonal.h \
    solver/utils/matrix.h \
    solver/utils/vector.h

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    resources.qrc
