#include <QApplication>
#include "mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    app.setApplicationName("ChikaterasuGUI");
    app.setOrganizationName("KyotoU-GradMed");

    MainWindow w;
    w.show();
    return app.exec();
}
