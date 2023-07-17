#include "mainwindow.h"

/** First tab Class Constructor. */
ChannelTimeResponse::ChannelTimeResponse(QWidget *parent)
    : QWidget(parent)
{
    /** Parameters of the process. */
    series = new QLineSeries();                          /** @brief series Stores the channel gains. */
    sig = new QLineEdit();                               /** @brief sig Standard deviation of the process. */
    fmax = new QLineEdit();                              /** @brief fmax Maximum doppler frequency. */
    finalTime = new QLineEdit();                         /** @brief finalTime Number of coherence time. */
    stepTime = new QLineEdit();                          /** @brief stepTime Time step. */

    /** Statistics of the process. */
    meanLE = new QLineEdit();                            /** @brief meanLE Mean of the channel. */
    varLE = new QLineEdit();                             /** @brief varLE Variance of the channel. */
    meanPLE = new QLineEdit();                           /** @brief meanPLE Mean of the process. */

    /** Objects of the widget. */
    funcType = new QComboBox;                            /***/
    chart = new QChart();                                /** Object to plot the series. */
    chartView = new QChartView(chart);                   /** Object to visualize the plot of the series. */
    m_button = new QPushButton("Refresh series", this);  /** Button to apply the parameters and compute the channel. */

    /** Setup of the Standard Deviation parameter.*/
    QLabel *sigLabel = new QLabel(tr("Standard deviation: "));
    sig->setText(QString::number(1));                                   /** Iniciates the standard deviation in 1 */
    sig->setValidator( new QDoubleValidator(0, 10, 4, this) );          /** Restricts the input value.*/

    /** Setup of the Maximum Doppler Frequency. */
    QLabel *fmaxLabel = new QLabel(tr("Maximum doppler frequency: "));
    fmax->setText(QString::number(91));                                 /** Iniciates the maximum doppler frequency in 91 */
    fmax->setValidator( new QDoubleValidator(10, 1000, 4, this) );      /** Restricts the input value. */

    /** Setup of the Coherence Time Number. */
    QLabel *finalLabel = new QLabel(tr("Number of coherence time: "));
    finalTime->setText(QString::number(5));                             /** Iniciates the coherence time number in 5. */
    finalTime->setValidator( new QDoubleValidator(1, 100, 3, this) );   /** Restricts the input value. */

    /** Setup of the Step Time. */
    QLabel *stepLabel = new QLabel(tr("Step time: "));
    stepTime->setText(QString::number(0.0001));                         /** Iniciates the step time in 0.0001 seconds. */
    stepTime->setValidator( new QDoubleValidator(0.1, 10, 6, this) );   /** Restricts the input value. */

    /**
     * @brief meanLabel, label of the mean of the channel.
     * @brief This variable is a "Read only", as shows the mean of the computed channel and it can't be edited.
     */
    QLabel *meanLabel = new QLabel(tr("Mean of the signal: "));
    meanLE->setReadOnly(true);
    meanLE->setText(QString::number(0));

    /**
     * @brief varLabel, label of the variance of the channel.
     * @brief This variable is a "Read only", as shows the variance of the computed channel and it can't be edited.
     */
    QLabel *varLabel = new QLabel(tr("Variance of the signal: "));
    varLE->setReadOnly(true);
    varLE->setText(QString::number(0));

    /**
     * @brief meanPLabel, label of the mean of the process.
     * @brief This variable is a "Read only", as shows the mean of the computed process and it can't be edited.
     */
    QLabel *meanPLabel = new QLabel(tr("Mean power of the signal: "));
    meanPLE->setReadOnly(true);
    meanPLE->setText(QString::number(0));

    /**
     * @brief Setup of the ComboBox, shows the methods of calculating the process available for choice.
     */
    QLabel *funcLabel = new QLabel(tr("Choose the model: "));
    funcType->addItem(tr("MED model"));
    funcType->addItem(tr("MEDS model"));

    /** Setup of the chart view. */
    chart->legend()->hide();                            /** Hides the legend of the figure. */
    chart->addSeries(series);                           /** Adds the channel gain graph to the figure. */
    chart->setTitle("Channel Time Response.");          /** Sets the title of the figure. */
    chartView->setRenderHint(QPainter::Antialiasing);

    /** Creating a connection to the refresh button to recalculate the channel when clicked. */
    QObject::connect(m_button, &QPushButton::clicked, this, &ChannelTimeResponse::recalcSeries);

    /**
     * @brief Declaring the function recalcSeries.
     */
    ChannelTimeResponse::recalcSeries();

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
    controlLayout->addWidget(finalLabel,3,0);
    controlLayout->addWidget(finalTime,3,1);
    controlLayout->addWidget(stepLabel,4,0);
    controlLayout->addWidget(stepTime,4,1);

    controlLayout->addWidget(meanLabel,5,0);
    controlLayout->addWidget(meanLE,5,1);
    controlLayout->addWidget(varLabel,6,0);
    controlLayout->addWidget(varLE,6,1);
    controlLayout->addWidget(meanPLabel,7,0);
    controlLayout->addWidget(meanPLE,7,1);
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
 * @brief recalcSeries, function that refreshes the plot of the channel.
 */
void ChannelTimeResponse::recalcSeries(){
    series->clear();                 /** Clear the series. */

    auto *x=new std::vector<float>;  /** @brief x a pointer, vector of floats that stores the values of the x axis. */
    std::vector<float> rp;           /** @brief rp a vector of floats, stores the real process component. */
    std::vector<float> ip;           /** @brief ip a vector of floats, stores the imaginary process component. */

    /**
     *  @brief Iniciating variables.
     *  @brief mean a double, the mean value of the channel.
     *  @brief var a double, the variance of the channel.
     *  @brief meanP a double, the mean value of the process.
     *  @brief ti a double, the time instant.
     *  @brief dt a double, the time step.
     *  @brief cohtime a double, the coherence time of the process.
     */
    double mean = 0,var = 0, meanP=0,ti=0,tf=finalTime->text().toFloat(),dt=stepTime->text().toFloat();
    double cohtime=1/fmax->text().toDouble();

    /** @brief t a vector of floats, the time vector. */
    std::vector<float> t ((tf*cohtime)/dt +1);

    for(int i=0;i<(tf*cohtime)/dt;i++) t[i]=ti + i*dt;

    MEDModel<20,float> u1(sig->text().toDouble(),fmax->text().toDouble());   /** Declares real process component for the MED model. */
    MEDModel<21,float> u1i(sig->text().toDouble(),fmax->text().toDouble());  /** Declares imaginary process component for the MED model. */

    MEDSModel<20,float> u2(sig->text().toDouble(),fmax->text().toDouble());  /** Declares real process component for the MEDS model. */
    MEDSModel<21,float> u2i(sig->text().toDouble(),fmax->text().toDouble()); /** Declares imaginary process component for the MEDS model. */

    /** If the chosen model is the MED. */
    if (funcType->currentText()==("MED model")){
        rp=u1.CalcProcess(t);                                            /** Computles the real process. */
        ip=u1i.CalcProcess(t);                                           /** Computes the imaginary process. */
        for(int i=0; i<=tf*cohtime/dt; i++){
            x->push_back(std::abs(std::complex<double>(rp[i],ip[i])));   /** Computes channel gain and adds it to the vector object. */
            series->append(t[i],10*log10(x->at(i)));                     /** Converts the gains to dB and includes includes it in the series. */
        }


    }
    /** If the chosen model is the MEDS. */
    else {
        rp=u2.CalcProcess(t);                                            /** Computles the real process. */
        ip=u2i.CalcProcess(t);                                           /** Computes the imaginary process. */
        for(int i=0; i<=tf*cohtime/dt; i++){
            x->push_back(std::abs(std::complex<double>(rp[i],ip[i])));   /** Computes channel gain and adds it to the vector object. */
            series->append(t[i],10*log10(x->at(i)));                     /** Converts the gains to dB and includes includes it in the series. */
        }


    }

    for (auto &v : *x)
        mean+=v;                /** Computes the mean of the channel. */
    mean/=x->size();

    for (auto &v : *x)
        var+=pow(v-mean,2);     /** Computes the variance of the channel. */
    var=(var/x->size());

    for (auto &v : *x)
        meanP+=pow(fabs(v),2);  /** Computes the mean of the process. */
    meanP/=x->size();


    meanLE->setText(QString::number(mean));     /** Sets the value of the mean of channel to show at the window. */
    varLE->setText(QString::number(var));       /** Sets the value of the variance of channel to show at the window. */
    meanPLE->setText(QString::number(meanP));   /** Sets the value of the mean of process to show at the window. */

    /** Setup of the horizontal axis of the plot. */
    axisX=new QValueAxis;
    axisX->setRange(0, tf*cohtime);
    axisX->setTickCount(10);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleText("t(ms)");
    chartView->chart()->setAxisX(axisX, series);

    /** Setup of the vertical axis of the plot. */
    axisY=new QValueAxis;
    axisY->setRange(10*log10(*min_element(x->begin(), x->end()))*1.05, 10*log10(*max_element(x->begin(), x->end()))*1.05);
    axisY->setTickCount(10);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleText("channel gain (dB)");
    chartView->chart()->setAxisY(axisY, series);

    /** Clear the vectors and objects. */
    rp.clear();
    ip.clear();
    x->clear();
    t.clear();
    delete x;
}
