#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QCheckBox>
#include <QPushButton>

class VTKWidget;
//class GLWidget;
class MainWindow;

class Window : public QWidget
{
    Q_OBJECT

public:
    Window(MainWindow *mw);

protected:
    void keyPressEvent(QKeyEvent *event) override;

private:
    //GLWidget *glWidget;
    VTKWidget *vtkWidget;
    MainWindow *mainWindow;
    QCheckBox *selectCamPos;
    QPushButton *loadModel;
    QPushButton *dilate;
    QPushButton *optimize;
};

#endif // WINDOW_H
