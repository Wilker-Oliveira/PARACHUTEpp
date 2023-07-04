#ifndef MAINWINDOW_H
#define MAINWINDOW_H

/** Necessary Qt libraries. */
#include <QMainWindow>
#include <QPushButton>
#include <QLineEdit>
#include <QtCharts>
#include <QComboBox>
#include <QTabWidget>

/** Standard C++ libraries. */
#include <cmath>
#include <vector>
#include <complex>
#include <chrono>
#include <random>

/** Channel model classes. */
#include "../MRC_1_IS/Models/JakesModel.hpp"
#include "../MRC_1_IS/Models/MCModel.hpp"
#include "../MRC_1_IS/Models/MEAModel.hpp"
#include "../MRC_1_IS/Models/MEDModel.hpp"
#include "../MRC_1_IS/Models/MEDSModel.hpp"
#include "../MRC_1_IS/Models/MSEModel.hpp"


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

/**
 * @brief The MainWindow class
 * @brief Manages all the events in the GUI.
*/
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    /**
     * @brief tabManager
     * @brief Creates the tabs of the GUI.
     */
    QTabWidget *tabManager;
    Ui::MainWindow *ui;

};

/**
 * @brief The ChannelTimeResponse class
 * @brief Object class related to the first tab.
 * @brief This tab shows the channel response response trhough time and the channel statistics given the chosen model and respective parameters.
 */
class ChannelTimeResponse : public QWidget
{
    Q_OBJECT

public:
    explicit ChannelTimeResponse(QWidget *parent = nullptr);

private:

    void recalcSeries();

    /** Object and widgets. */
    QComboBox *funcType = nullptr;
    QPushButton *m_button = nullptr;
    QChart *chart = nullptr;
    QChartView *chartView = nullptr;
    QLineSeries *series = nullptr;

    /** Axis control. */
    QValueAxis *axisX = nullptr;
    QValueAxis *axisY = nullptr;

    /** Paramters. */
    QLineEdit *sig = nullptr;         /** Standard deviation. */
    QLineEdit *fmax = nullptr;        /** Maximum doppler frequency. */
    QLineEdit *finalTime = nullptr;   /** Number of coherence time. */
    QLineEdit *stepTime = nullptr;    /** Time step. */

    /** Statistics. */
    QLineEdit *meanLE = nullptr;      /** Mean of the components. */
    QLineEdit *varLE = nullptr;       /** Variance of the process. */
    QLineEdit *meanPLE = nullptr;     /** Mean of the process*/
};

//Second tab object Class
class ChannelAutoCorrelation : public QWidget
{
    Q_OBJECT

public:
    explicit ChannelAutoCorrelation(QWidget *parent = nullptr);

private:

    void calcCorr();

    //Object and widgets
    QPushButton *m_button = nullptr;
    QChart *chart = nullptr;
    QChartView *chartView = nullptr;
    QLineSeries *Parametric_series = nullptr;
    QLineSeries *Empiric_series = nullptr;

    //Axis control
    QValueAxis *axisX = nullptr;
    QValueAxis *axisY = nullptr;

    //Paramters
    QLineEdit *sig = nullptr;
    QLineEdit *fmax = nullptr;
    QLineEdit *finalTime = nullptr;
    QLineEdit *stepTime = nullptr;

};


//Useful Function on for calculate the Theoretical ACF
std::vector<float> envACF(float sig, std::vector<float> acf1, std::vector<float> acf2);

//Useful Function on for calculate the Empirical ACF
std::vector<std::complex<float>> autoCorr(std::vector<std::complex<float>> u, float lag);

//Third tab object Class
class ChannelFrequencyResponse : public QWidget
{
    Q_OBJECT

public:
    explicit ChannelFrequencyResponse(QWidget *parent = nullptr);

private:

    void calcFreqRes();

    //Object and widgets
    QPushButton *m_button = nullptr;
    QChart *chart = nullptr;
    QChartView *chartView = nullptr;
    QScatterSeries *Parametric_series = nullptr;

    //Axis control
    QValueAxis *axisX = nullptr;
    QValueAxis *axisY = nullptr;

    //Paramters
    QLineEdit *sig = nullptr;
    QLineEdit *fmax = nullptr;

};


#endif // MAINWINDOW_H
