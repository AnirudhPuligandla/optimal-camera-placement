//#include "glwidget.h"
#include "vtkwidget.h"
#include "window.h"
#include "mainwindow.h"
#include <QKeyEvent>
#include <QDesktopWidget>
#include <QApplication>
#include <QMessageBox>
#include <QHBoxLayout>


Window::Window(MainWindow *mw) : mainWindow(mw)
{
    //glWidget = new GLWidget;
    vtkWidget = new VTKWidget();
    selectCamPos = new QCheckBox("select vertices");
    loadModel = new QPushButton("Load model");
    dilate = new QPushButton("Initialize data");
    optimize = new QPushButton("Optimize");

    //connect(selectCamPos, &QCheckBox::toggled, glWidget, &GLWidget::setSelection);
    connect(loadModel, &QPushButton::pressed, vtkWidget, &VTKWidget::loadModelPressed);
    connect(dilate, &QPushButton::pressed, vtkWidget, &VTKWidget::dilatePressed);
    //connect(optimize, &QPushButton::pressed, glWidget, &GLWidget::optimizePressed);

    QVBoxLayout *container = new QVBoxLayout;
    container->addWidget(vtkWidget);
    container->addWidget(selectCamPos);
    container->addWidget(loadModel);
    container->addWidget(dilate);
    container->addWidget(optimize);
    setLayout(container);
    show();

    setWindowTitle(tr("OCP Simulation"));
}

void Window::keyPressEvent(QKeyEvent *e)
{
    if (e->key() == Qt::Key_Escape)
        close();
    else
        QWidget::keyPressEvent(e);
}
