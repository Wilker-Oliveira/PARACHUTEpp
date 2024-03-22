#include "mainwindow.h"

/** Third tab Class Constructor. */
ChannelFrequencyResponse::ChannelFrequencyResponse(QWidget *parent)
    : QWidget(parent)
{
    /** Parameters of the process. */
    Parametric_series = new QScatterSeries();               /** @brief Parametric_Series stores the parametric doppler spectrum of the channel. */
    Empiric_series = new QScatterSeries();
    sig = new QLineEdit();                                  /** @brief sig Standard deviation of the process. */
    fmax = new QLineEdit();                                 /** @brief fmax Maximum doppler frequency. */

    /** Objects of the widget. */
    funcType = new QComboBox;                            /***/
    chart = new QChart();
    chartView = new QChartView(chart);
    m_button = new QPushButton("Refresh series", this);

    /** Setup of the Standard Deviation parameter.*/
    QLabel *sigLabel = new QLabel(tr("Standard deviation: "));
    sig->setText(QString::number(1));                              /** Iniciates the standard deviation in 0.7071 */
    sig->setValidator( new QDoubleValidator(0, 10, 4, this) );          /** Restricts the input value.*/

    /** Setup of the Maximum Doppler Frequency. */
    QLabel *fmaxLabel = new QLabel(tr("Maximum doppler frequency: "));
    fmax->setText(QString::number(91));                                 /** Iniciates the maximum doppler frequency in 21 */
    fmax->setValidator( new QDoubleValidator(10, 1000, 4, this) );      /** Restricts the input value. */

    /**
     * @brief Setup of the ComboBox, shows the methods of calculating the process available for choice.
     */
    QLabel *funcLabel = new QLabel(tr("Choose the model: "));
    funcType->addItem(tr("MED model"));
    funcType->addItem(tr("MEDS model"));
    funcType->addItem(tr("MSE model"));
    funcType->addItem(tr("MEA model"));
    funcType->addItem(tr("MCM model"));

    /** Setup of the chart view. */
    chart->legend()->hide();                                            /** Hides the legend of the figure. */
    chart->addSeries(Parametric_series);                                /** Adds the parametric doppler spectrum graph to the figure. */
    chart->addSeries(Empiric_series);
    chart->setTitle("Channel Doppler Spectrum.");                       /** Sets the title of the figure. */
    chartView->setRenderHint(QPainter::Antialiasing);

    /** Creating a connection to the refresh button to recalculate the doppler spectrum when clicked. */
    QObject::connect(m_button, &QPushButton::clicked, this, &ChannelFrequencyResponse::calcFreqRes);

    /**
     * @brief Declaring the function calcFreqRes.
     */
    ChannelFrequencyResponse::calcFreqRes();

    /**
     * @brief controlWidget, object that ...
     * @brief controlLayout, object that manages the position of the objects in the widget.
     */
    QWidget *controlWidget = new QWidget(this);
    QGridLayout *controlLayout = new QGridLayout(controlWidget);
    controlLayout->addWidget(funcLabel,0,0);
    controlLayout->addWidget(funcType,0,1);
    controlLayout->addWidget(sigLabel,1,0);
    controlLayout->addWidget(sig,1,1);
    controlLayout->addWidget(fmaxLabel,2,0);
    controlLayout->addWidget(fmax,2,1);


    controlLayout->addWidget(m_button,8,0,1,2);

    /**
     * @brief mainLayout object that disposes the widgets horizontally inside the window.
     */
    QHBoxLayout *mainLayout = new QHBoxLayout;

    mainLayout->addWidget(chartView);
    mainLayout->setStretch(0, 1);
    mainLayout->addWidget(controlWidget);

    setLayout(mainLayout);
}

/**
 * @brief calcFreqRes, function that refreshes the channel and the doppler spectrum.
 */
 void ChannelFrequencyResponse::calcFreqRes(){

    /** Clear the series. */
    Parametric_series->clear();
    Empiric_series->clear();

    /**
     * @brief Iniciating the variables necessary for the computation of the doppler spectrum.
     * @brief flim a double, the maximum doppler frequency.
     */
    double flim=fmax->text().toDouble()+1;

    std::vector<float> rp;           /** @brief rp a vector of floats, stores the real process component. */
    std::vector<float> ip;           /** @brief ip a vector of floats, stores the imaginary process component. */

    MEDModel<25,float> u1(sig->text().toDouble(),fmax->text().toDouble());   /** Declares real process component for the MED model. */
    MEDModel<25,float> u1i(sig->text().toDouble(),fmax->text().toDouble());  /** Declares imaginary process component for the MED model. */

    MEDSModel<25,float> u2(sig->text().toDouble(),fmax->text().toDouble());  /** Declares real process component for the MEDS model. */
    MEDSModel<25,float> u2i(sig->text().toDouble(),fmax->text().toDouble()); /** Declares imaginary process component for the MEDS model. */

    MSEModel<25,float> u3(sig->text().toDouble(),fmax->text().toDouble());  /** Declares real process component for the MEDS model. */
    MSEModel<25,float> u3i(sig->text().toDouble(),fmax->text().toDouble()); /** Declares imaginary process component for the MEDS model. */

    MEAModel<25,float> u4(sig->text().toDouble(),fmax->text().toDouble());  /** Declares real process component for the MEDS model. */
    MEAModel<25,float> u4i(sig->text().toDouble(),fmax->text().toDouble()); /** Declares imaginary process component for the MEDS model. */

    MCModel<25,float> u5(sig->text().toDouble(),fmax->text().toDouble());  /** Declares real process component for the MEDS model. */
    MCModel<25,float> u5i(sig->text().toDouble(),fmax->text().toDouble()); /** Declares imaginary process component for the MEDS model. */


    std::vector<double> x;                                                    /** @brief x vector of floats, stores the values of the x axis. */
    std::vector<double> f;                                                    /** @brief f vector of floats, stores the values of the y axis. */

    /** Set size of the vectors. */
    x.reserve(50);
    f.reserve(50);

    /** If the chosen model is the MED. */
    if (funcType->currentText()==("MED model")){
        rp = u1.CalcDiscreteDopplerSpectrum();   /** @brief rp vector of floats, stores the doppler spectrum values of the real process. */
        ip = u1i.CalcDiscreteDopplerSpectrum();  /** @brief rp vector of floats, stores the doppler spectrum values of the imaginary process. */
    }

    /** If the chosen model is the MEDS. */
    else if (funcType->currentText()==("MEDS model")){
        rp = u2.CalcDiscreteDopplerSpectrum();   /** @brief rp vector of floats, stores the doppler spectrum values of the real process. */
        ip = u2i.CalcDiscreteDopplerSpectrum();  /** @brief rp vector of floats, stores the doppler spectrum values of the imaginary process. */
    }

    /** If the chosen model is the MSE. */
    else if (funcType->currentText()==("MSE model")){
        rp = u3.CalcDiscreteDopplerSpectrum();   /** @brief rp vector of floats, stores the doppler spectrum values of the real process. */
        ip = u3i.CalcDiscreteDopplerSpectrum();  /** @brief rp vector of floats, stores the doppler spectrum values of the imaginary process. */
    }

    /** If the chosen model is the MEA. */
    else if (funcType->currentText()==("MEA model")){
        rp = u4.CalcDiscreteDopplerSpectrum();   /** @brief rp vector of floats, stores the doppler spectrum values of the real process. */
        ip = u4i.CalcDiscreteDopplerSpectrum();  /** @brief rp vector of floats, stores the doppler spectrum values of the imaginary process. */
    }

    /** If the chosen model is the MCM. */
    else {
        rp = u5.CalcDiscreteDopplerSpectrum();   /** @brief rp vector of floats, stores the doppler spectrum values of the real process. */
        ip = u5i.CalcDiscreteDopplerSpectrum();  /** @brief rp vector of floats, stores the doppler spectrum values of the imaginary process. */
    }

    /** Computes the Parametric PSD. */
    for(int i=0;i < rp.size()/2;i++){

        f[i]=rp[2*i];
        x[i]=2*rp[1+2*i];
        Parametric_series->append(f[i],(x[i]));               /** Adds the discrete doppler spectrum values in the series. */
    }


    Parametric_series->setMarkerSize(10);

    /** Setup of the horizontal axis of the plot. */
    axisX=new QValueAxis;
    axisX->setRange(-flim, flim);
    axisX->setTickCount(10);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleText("f(Hz)");
    chartView->chart()->setAxisX(axisX, Parametric_series);

    /** Setup of the vertical axis of the plot. */
    axisY=new QValueAxis;
    axisY->setRange(*min_element(x.begin(), x.end())*0.95, *max_element(x.begin(), x.end())*2.8);
    axisY->setTickCount(10);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleText("Channel frequency response");
    chartView->chart()->setAxisY(axisY, Parametric_series);

    /** Clearing the variables and objects of the widget.*/
    rp.clear();
    ip.clear();
    x.clear();
    f.clear();
}
