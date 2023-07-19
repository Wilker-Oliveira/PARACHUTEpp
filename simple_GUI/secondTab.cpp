#include "mainwindow.h"

/** Second tab Class Constructor. */
ChannelAutoCorrelation::ChannelAutoCorrelation(QWidget *parent)
    : QWidget(parent)
{
    /** Parameters of the process. */
    Parametric_series = new QLineSeries();               /** @brief Parametric_Series stores the parametric autocorrelation of the channel. */
    Empiric_series = new QLineSeries();                  /** @brief Empiric_series stores the empirical autocorrelation of the channel. */
    sig = new QLineEdit();                               /** @brief sig Standard deviation of the process. */
    fmax = new QLineEdit();                              /** @brief fmax Maximum doppler frequency. */
    finalTime = new QLineEdit();                         /** @brief finalTime Number of coherence time. */
    stepTime = new QLineEdit();                          /** @brief stepTime Time step. */

    /** Objects of the widget. */
    chart = new QChart();
    chartView = new QChartView(chart);
    m_button = new QPushButton("Refresh series", this);

    /** Setup of the Standard Deviation parameter.*/
    QLabel *sigLabel = new QLabel(tr("Standard deviation: "));
    sig->setText(QString::number(1));                              /** Iniciates the standard deviation in 0.7071 */
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
    stepTime->setText(QString::number(0.0005));                         /** Iniciates the step time in 0.0005 seconds. */
    stepTime->setValidator( new QDoubleValidator(0.1, 10, 6, this) );   /** Restricts the input value. */

    /** Setup of the chart view. */
    chart->legend()->hide();                                            /** Hides the legend of the figure. */
    chart->addSeries(Parametric_series);                                /** Adds the parametric autocorrelation graph to the figure. */
    chart->addSeries(Empiric_series);                                   /** Adds the empirical autocorrelation graph to the figure. */
    chart->setTitle("Channel Autocorrelation.");                        /** Sets the title of the figure. */
    chartView->setRenderHint(QPainter::Antialiasing);

    /** Creating a connection to the refresh button to recalculate the autocorrelation when clicked. */
    QObject::connect(m_button, &QPushButton::clicked, this, &ChannelAutoCorrelation::calcCorr);

    /**
     * @brief Declaring the function calcCorr.
     */
    ChannelAutoCorrelation::calcCorr();

    /**
     * @brief controlWidget, object that ...
     * @brief controlLayout, object that manages the position of the objects in the widget.
     */
    QWidget *controlWidget = new QWidget(this);
    QGridLayout *controlLayout = new QGridLayout(controlWidget);
    controlLayout->addWidget(sigLabel,1,0);
    controlLayout->addWidget(sig,1,1);
    controlLayout->addWidget(fmaxLabel,2,0);
    controlLayout->addWidget(fmax,2,1);
    controlLayout->addWidget(finalLabel,3,0);
    controlLayout->addWidget(finalTime,3,1);
    controlLayout->addWidget(stepLabel,4,0);
    controlLayout->addWidget(stepTime,4,1);

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
 * @brief calcCorr, function that refreshes the channel and the autocorrelations.
 */
void ChannelAutoCorrelation::calcCorr(){

    /** Clear the series. */
    Parametric_series->clear();
    Empiric_series->clear();

    /**
     *  @brief Iniciating variables to store the time instant of the process.
     *  Receives the input values of the parameters.
     *  @brief ti a double, the time instant.
     *  @brief tf a double, the final time instant.
     *  @brief dt a double, the time step.
     *  @brief cohtime a double, the coherence time.
     */
    double ti=0,tf=finalTime->text().toFloat(),dt=stepTime->text().toFloat();
    double cohtime=1/fmax->text().toDouble();

    /** @brief s_length a double, the time length. */
    double s_length=ceil((tf*cohtime)/dt);

    auto *x=new std::vector<float>;           /** @brief x a pointer, vector of floats that stores the values of the x axis. */
    std::vector<float> t (s_length+1);        /** @brief t a vector of floats, the time vector. */
    std::vector<float> rp (s_length+1);       /** @brief rp a vector of floats, stores the real process component. */
    std::vector<float> ip (s_length+1);       /** @brief rp a vector of floats, stores the imaginary process component. */


    MEDModel<25,float> u1(sig->text().toDouble(),fmax->text().toDouble());   /** Declares real process component for the MEA model. */
    MEDModel<25,float> u1i(sig->text().toDouble(),fmax->text().toDouble());  /** Declares imaginary process component for the MEA model. */

    /** Computes the parametric autocorrelation. */
    for(int i=0;i<s_length;i++){
        t[i]=ti + i*dt;            /** Time instant. */
        rp[i]=u1.CalcACF(t[i]);    /** Autocorrelation of the real process component. */
        ip[i]=u1i.CalcACF(t[i]);   /** Autocorrelation of the imaginary process component. */
    }
    *x=envACF(sig->text().toDouble(),rp,ip);  /** Autocorrelation of the complex envelope. */

    double s_norm = x->at(0);
    for(int i=0; i<=s_length-1; i++){
        Parametric_series->append(t[i],(x->at(i))/s_norm);  /** Adds the autocorrelation values to the series. */
        x->at(i)/=s_norm;
    }
    Parametric_series->setName(QString("Parametric autocorrelation"));


    std::vector<float> u1value = u1.CalcProcess(t);    /** Computes the real process. */
    std::vector<float> u2value = u1i.CalcProcess(t);   /** Computes the imaginary process. */

    std::vector<std::complex<float>> S_ACFvec;         /** Vector to store the empirical autocorrelation values. */
    std::vector<std::complex<float>> uvalue,aux;       /** Vector of the channel realizations. */
    std::complex<float> ch;                            /** Vector to store the channel values. */

    for(int i = 0; i < u1value.size(); i++){
        ch = std::complex<float>(u1value[i], u2value[i]); /** Computes and stores the channel values. */
        uvalue.push_back(ch);
    }

    S_ACFvec = autoCorr(uvalue,s_length+1);    /** Takes the empirical autocorrelation of the channel. */
    uvalue.clear();                            /** Clear complex envelope vector. */

    /** Autocorrelation of 200 realizations. */
    for(int i =0;i<200;i++){
        /** Generate the random phases */
        u1.genPhases();
        u1i.genPhases();

        u1value = u1.CalcProcess(t);       /** Computes the real process. */
        u2value = u1i.CalcProcess(t);      /** Computes the imaginary process. */

        /** Generate the complex envelope */
        for(int i = 0; i < u1value.size(); i++){
            ch = std::complex<float>(u1value[i], u2value[i]); /** Computes and stores the channel values. */
            uvalue.push_back(ch);
        }

        aux = autoCorr(uvalue,s_length+1); /** Takes the empirical autocorrelation of the channel in auxiliar vector. */

        uvalue.clear(); /** Clear complex envelope vector. */

        for(float j=0;j<s_length+1;j++){
            S_ACFvec[j]= (S_ACFvec[j]+aux[j]); /** Stores the autocorrelation of the channel. */
        }

        aux.clear(); /** Clear auxiliar vector. */
    }

    s_norm = std::real(S_ACFvec[0]);                             /** Takes the norm as the real part of the first value of the autocorrelation. */
    for(int i=0; i<=s_length-1; i++){
        Empiric_series->append(t[i],std::real(S_ACFvec[i])/s_norm);     /** Adds the normalized autocorrelation values to the series. */
    }
    Empiric_series->setName(QString("Empirical autocorrelation with 200 realizations"));

    /** Setup of the horizontal axis of the plot. */
    axisX=new QValueAxis;
    axisX->setRange(0, tf*cohtime);
    axisX->setTickCount(10);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleText("t(ms)");
    chartView->chart()->setAxisX(axisX, Parametric_series);
    chartView->chart()->setAxisX(axisX, Empiric_series);

    /** Setup of the vertical axis of the plot. */
    axisY=new QValueAxis;
    axisY->setRange(*min_element(x->begin(), x->end())*1.05, *max_element(x->begin(), x->end())*1.05);
    axisY->setTickCount(10);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleText("channel response correlation");
    chartView->chart()->setAxisY(axisY, Parametric_series);
    chartView->chart()->setAxisY(axisY, Empiric_series);
    chartView->chart()->legend()->setVisible(true);

    /** Clearing the variables and objects of the widget.*/
    rp.clear();
    ip.clear();
    x->clear();
    t.clear();
    delete x;
}

/**
 * @brief envACF function that computes the theoretical autocorrelation.
 * @param sig float, standard deviation of the channel.
 * @param acf1 vector of floats, autocorrelation of the real component of the process.
 * @param acf2 vector of floats, autocorrelation of the imaginary component of the process.
 * @return vector of floats containing the autocorrelation of the process.
 */
std::vector<float> envACF(float sig, std::vector<float> acf1, std::vector<float> acf2){
    std::vector<float> en_acf;
    float aux;

    for(float i = 0; i< acf1.size(); i++){
        aux = acf1[i] + acf2[i];
        en_acf.push_back(aux);
    }
    return en_acf;
}

/**
 * @brief autoCorr function that computes the unormalized empirical autocorrelation.
 * @param u vector of complex floats, the channel realization.
 * @param lag float, delay between time instants.
 * @return vector of complex floats containing the autocorrelation of the process.
 */
std::vector<std::complex<float>> autoCorr(std::vector<std::complex<float>> u, float lag){

    std::complex<float> aux (0.0, 0.0);
    std::vector<std::complex<float>> ac;

    for(float i = 0; i < lag; i++){
        for(float j = 0; j < lag-i; j++){
            aux += u[j+i]*std::conj(u[j]);
        }
        ac.push_back(aux);
        aux = 0;
    }
    return ac;
}
