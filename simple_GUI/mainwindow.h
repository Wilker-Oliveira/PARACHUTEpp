#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QLineEdit>
#include <QtCharts>
#include <QComboBox>

#include <cmath>
#include <complex>
#include <chrono>
#include <random>
#include "../MRC_1_IS/Models/JakesModel.hpp"
#include "../MRC_1_IS/Models/MCModel.hpp"
#include "../MRC_1_IS/Models/MEAModel.hpp"
#include "../MRC_1_IS/Models/MEDModel.hpp"
#include "../MRC_1_IS/Models/MEDSModel.hpp"
#include "../MRC_1_IS/Models/MSEModel.hpp"


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:

    void recalcSeries();

    Ui::MainWindow *ui;

    //Object and widgets
    QComboBox *funcType = nullptr;
    QPushButton *m_button = nullptr;
    QChart *chart = nullptr;
    QChartView *chartView = nullptr;
    QLineSeries *series = nullptr;

    //Axis control
    QValueAxis *axisX = nullptr;
    QValueAxis *axisY = nullptr;

    //Paramters
    QLineEdit *sig = nullptr;
    QLineEdit *fmax = nullptr;
    QLineEdit *finalTime = nullptr;
    QLineEdit *stepTime = nullptr;

    //Statistics
    QLineEdit *meanLE = nullptr;
    QLineEdit *varLE = nullptr;
    QLineEdit *meanPLE = nullptr;


};

#endif // MAINWINDOW_H
