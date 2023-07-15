#include "mainwindow.h"

//Second tab Class Constructor
ChannelAutoCorrelation::ChannelAutoCorrelation(QWidget *parent)
    : QWidget(parent)
{
    //QT init
    Parametric_series = new QLineSeries();
    Empiric_series = new QLineSeries();
    sig = new QLineEdit();
    fmax = new QLineEdit();
    finalTime = new QLineEdit();
    stepTime = new QLineEdit();

    chart = new QChart();
    chartView = new QChartView(chart);
    m_button = new QPushButton("Refresh series", this);

    //Parameters setup
    QLabel *sigLabel = new QLabel(tr("Standard deviation: "));
    sig->setText(QString::number(0.7071));
    sig->setValidator( new QDoubleValidator(0, 10, 4, this) );

    QLabel *fmaxLabel = new QLabel(tr("Maximum doppler frequency: "));
    fmax->setText(QString::number(91));
    fmax->setValidator( new QDoubleValidator(10, 1000, 4, this) );

    QLabel *finalLabel = new QLabel(tr("Number of coherence time: "));
    finalTime->setText(QString::number(5));
    finalTime->setValidator( new QDoubleValidator(1, 100, 3, this) );

    QLabel *stepLabel = new QLabel(tr("Step time: "));
    stepTime->setText(QString::number(0.0005));
    stepTime->setValidator( new QDoubleValidator(0.1, 10, 6, this) );

    //Chart/Chartview setup
    chart->legend()->hide();
    chart->addSeries(Parametric_series);
    chart->addSeries(Empiric_series);
    chart->setTitle("Channel Autocorrelation.");
    chartView->setRenderHint(QPainter::Antialiasing);

    //Refresh button logic
    QObject::connect(m_button, &QPushButton::clicked, this, &ChannelAutoCorrelation::calcCorr);

    ChannelAutoCorrelation::calcCorr();

    //Positioning Setup
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

    QHBoxLayout *mainLayout = new QHBoxLayout;

    mainLayout->addWidget(chartView);
    mainLayout->setStretch(0, 1);
    mainLayout->addWidget(controlWidget);

    setLayout(mainLayout);

}

//Second tab Class refresh graph function
void ChannelAutoCorrelation::calcCorr(){

    Parametric_series->clear();
    Empiric_series->clear();

    //init variables
    double ti=0,tf=finalTime->text().toFloat(),dt=stepTime->text().toFloat();
    double cohtime=1/fmax->text().toDouble();

    double s_lenght=ceil((tf*cohtime)/dt);

    auto *x=new std::vector<float>;
    std::vector<float> t (s_lenght+1);
    std::vector<float> rp (s_lenght+1);
    std::vector<float> ip (s_lenght+1);


    MEAModel<20,float> u1(sig->text().toDouble(),fmax->text().toDouble());
    MEAModel<22,float> u1i(sig->text().toDouble(),fmax->text().toDouble());

    //calc Theoretical Corr
    for(int i=0;i<s_lenght;i++){
        t[i]=ti + i*dt;
        rp[i]=u1.CalcACF(t[i]);
        ip[i]=u1i.CalcACF(t[i]);
    }
    *x=envACF(sig->text().toDouble(),rp,ip);

    for(int i=0; i<=s_lenght-1; i++){
        Parametric_series->append(t[i],(x->at(i)));
    }
    Parametric_series->setName(QString("Parametric autocorrelation"));


    std::vector<float> u1value = u1.CalcProcess(t);
    std::vector<float> u2value = u1i.CalcProcess(t);

    std::vector<std::complex<float>> S_ACFvec;    // vector to store the empirical autocorrelation values
    std::vector<std::complex<float>> uvalue,aux;
    std::complex<float> ch;

    for(int i = 0; i < u1value.size(); i++){
        ch = std::complex<float>(u1value[i], u2value[i]);
        uvalue.push_back(ch);
    }

    S_ACFvec = autoCorr(uvalue,s_lenght+1);
    uvalue.clear();

    /**Taking the mean of the empirical autocorrelation function.*/

    for(int i =0;i<200;i++){
        // Generate the random phases
        u1.genPhases();
        u1i.genPhases();

        // Compute the processes
        u1value = u1.CalcProcess(t);
        u2value = u1i.CalcProcess(t);

        // Generate the complex envelope
        for(int i = 0; i < u1value.size(); i++){
            ch = std::complex<float>(u1value[i], u2value[i]);
            uvalue.push_back(ch);
        }

        // Compute the empirical autocorrelation in auxiliar vector
        aux = autoCorr(uvalue,s_lenght+1);

        // Clear complex envelope vector
        uvalue.clear();

        // Store the autocorrelation function result
        for(float j=0;j<s_lenght+1;j++){
            S_ACFvec[j]= (S_ACFvec[j]+aux[j]);
        }

        // Clear auxiliar vector
        aux.clear();
    }

    double s_norm = std::real(S_ACFvec[0]);
    for(int i=0; i<=s_lenght-1; i++){
        Empiric_series->append(t[i],std::real(S_ACFvec[i])/s_norm);
    }
    Empiric_series->setName(QString("Empiric autocorrelation with 200 realizations"));

    axisX=new QValueAxis;
    axisX->setRange(0, tf*cohtime);
    axisX->setTickCount(10);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleText("t(ms)");
    chartView->chart()->setAxisX(axisX, Parametric_series);
    chartView->chart()->setAxisX(axisX, Empiric_series);

    axisY=new QValueAxis;
    axisY->setRange(*min_element(x->begin(), x->end())*1.05, *max_element(x->begin(), x->end())*1.05);
    axisY->setTickCount(10);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleText("channel response correlation");
    chartView->chart()->setAxisY(axisY, Parametric_series);
    chartView->chart()->setAxisY(axisY, Empiric_series);
    chartView->chart()->legend()->setVisible(true);

    rp.clear();
    ip.clear();
    x->clear();
    t.clear();
    delete x;
}

//Useful utility function for calculate the Theoretical ACF
std::vector<float> envACF(float sig, std::vector<float> acf1, std::vector<float> acf2){
    std::vector<float> en_acf;
    float aux;

    for(float i = 0; i< acf1.size(); i++){
        aux = acf1[i] + acf2[i];
        en_acf.push_back(aux);
    }
    return en_acf;
}

//Useful utility function for calculate the Empirical ACF
//The output of this function is unormalized!!!!
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
