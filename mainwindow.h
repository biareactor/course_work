#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "solver.h"
#include <tuple>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    void plot_eq();

    ~MainWindow();

private:
    Ui::MainWindow *ui;

    Solver::Params get_params();
//    bool save_csv(const std::string& dir_name, const Solver::Result& res, const std::string& file_name = "file");
};

#endif // MAINWINDOW_H
