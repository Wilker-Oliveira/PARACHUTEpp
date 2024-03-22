#ifndef MAINWINDOW_H
#define MAINWINDOW_H

/** Necessary Qt libraries. */
#include <fftw3.h>
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
 * @brief The ChannelTimeResponse class.
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
    QLineEdit *meanPLE = nullptr;     /** Mean of the process.*/
};

/**
 * @brief The ChannelAutoCorrelation class.
 * @brief Object related to the second tab.
 * @brief This tabs shows the channel's parametrical and empirical autocorrelation.
 */
class ChannelAutoCorrelation : public QWidget
{
    Q_OBJECT

public:
    explicit ChannelAutoCorrelation(QWidget *parent = nullptr);

private:

    void calcCorr();

    /** Object and widgets. */
    QPushButton *m_button = nullptr;
    QChart *chart = nullptr;
    QChartView *chartView = nullptr;
    QLineSeries *Parametric_series = nullptr;
    QLineSeries *Empiric_series = nullptr;

    /** Axis control. */
    QValueAxis *axisX = nullptr;
    QValueAxis *axisY = nullptr;

    /** Paramters. */
    QLineEdit *sig = nullptr;          /** Standard deviation. */
    QLineEdit *fmax = nullptr;         /** Maximum doppler frequency. */
    QLineEdit *finalTime = nullptr;    /** Number of coherence time. */
    QLineEdit *stepTime = nullptr;     /** Time step. */

};


/**
 * @brief envACF function that computes the theoretical autocorrelation.
 * @param sig float, standard deviation of the channel.
 * @param acf1 vector of floats, autocorrelation of the real component of the process.
 * @param acf2 vector of floats, autocorrelation of the imaginary component of the process.
 * @return vector of floats containing the autocorrelation of the process.
 */
std::vector<float> envACF(float sig, std::vector<float> acf1, std::vector<float> acf2);

/**
 * @brief autoCorr function that computes the empirical autocorrelation.
 * @param u vector of complex floats, the channel realization.
 * @param lag float, delay between time instants.
 * @return vector of complex floats containing the autocorrelation of the process.
 */
std::vector<std::complex<float>> autoCorr(std::vector<std::complex<float>> u, float lag);

void fft_shift(std::complex<double>* A, unsigned long size);

/**
 * @brief The ChannelFrequencyResponse class.
 * @brief Object class relative to the third tab.
 * @brief This tab shows the plot of the doppler frequency spectrum.
 */
class ChannelFrequencyResponse : public QWidget
{
    Q_OBJECT

public:
    explicit ChannelFrequencyResponse(QWidget *parent = nullptr);

private:

    void calcFreqRes();

    /** Object and widgets. */
    QComboBox *funcType = nullptr;
    QPushButton *m_button = nullptr;
    QChart *chart = nullptr;
    QChartView *chartView = nullptr;
    QScatterSeries *Parametric_series = nullptr;
    QScatterSeries *Empiric_series = nullptr;

    /** Axis control. */
    QValueAxis *axisX = nullptr;
    QValueAxis *axisY = nullptr;

    /** Paramters. */
    QLineEdit *sig = nullptr;          /** Standard deviation. */
    QLineEdit *fmax = nullptr;         /** Maximum doppler frequency. */

};


#endif // MAINWINDOW_H
