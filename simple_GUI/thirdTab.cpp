//Third tab Class Constructor
ChannelFrequencyResponse::ChannelFrequencyResponse(QWidget *parent)
    : QWidget(parent)
{
    //QT init
    series = new QLineSeries();
    sig = new QLineEdit();
    fmax = new QLineEdit();

    chart = new QChart();
    chartView = new QChartView(chart);
    m_button = new QPushButton("Refresh series", this);

    //Parameters setup
    QLabel *sigLabel = new QLabel(tr("Standard deviation: "));
    sig->setText(QString::number(0.7071));
    sig->setValidator( new QDoubleValidator(0, 10, 4, this) );

    QLabel *fmaxLabel = new QLabel(tr("Maximum doppler frequency: "));
    fmax->setText(QString::number(21));
    fmax->setValidator( new QDoubleValidator(10, 1000, 4, this) );

    //Chart/Chartview setup
    chart->legend()->hide();
    chart->addSeries(series);
    chart->setTitle("I hate to admit it, but this interface is really useful");
    chartView->setRenderHint(QPainter::Antialiasing);

    //Refresh button logic
    QObject::connect(m_button, &QPushButton::clicked, this, &ChannelFrequencyResponse::calcFreqRes);

    ChannelFrequencyResponse::calcFreqRes();

    //Positioning Setup
    QWidget *controlWidget = new QWidget(this);
    QGridLayout *controlLayout = new QGridLayout(controlWidget);
    controlLayout->addWidget(sigLabel,1,0);
    controlLayout->addWidget(sig,1,1);
    controlLayout->addWidget(fmaxLabel,2,0);
    controlLayout->addWidget(fmax,2,1);


    controlLayout->addWidget(m_button,8,0,1,2);

    QHBoxLayout *mainLayout = new QHBoxLayout;

    mainLayout->addWidget(chartView);
    mainLayout->setStretch(0, 1);
    mainLayout->addWidget(controlWidget);

    setLayout(mainLayout);
}

//Third tab Class refresh graph function
void ChannelFrequencyResponse::calcFreqRes(){

    series->clear();

    //init variables
    double dt=0.01;
    double flim=fmax->text().toDouble()+1;

    std::vector<float> x (2*(flim/dt+1));
    std::vector<float> f (2*(flim/dt+1));
    std::vector<float> rp (2*(flim/dt+1));
    std::vector<float> ip (2*(flim/dt+1));

    MEAModel<20,float> u1(sig->text().toDouble(),fmax->text().toDouble());
    MEAModel<21,float> u1i(sig->text().toDouble(),fmax->text().toDouble());

    //calc parametric PSD
    for(int i=0;i < (2*(flim/dt+1));i++){
        f[i]= -flim + i*dt;
        rp[i]=u1.CalcPSD(f[i]);
        ip[i]=u1i.CalcPSD(f[i]);
        x[i]=rp[i]+ip[i];
        series->append(f[i],(x[i]));
    }

    axisX=new QValueAxis;
    axisX->setRange(-flim, flim);
    axisX->setTickCount(10);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleText("t(ms)");
    chartView->chart()->setAxisX(axisX, series);

    axisY=new QValueAxis;
    axisY->setRange(*min_element(x.begin(), x.end())*1.05, *max_element(x.begin(), x.end())*1.05);
    axisY->setTickCount(10);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleText("channel Frequency response");
    chartView->chart()->setAxisY(axisY, series);

    rp.clear();
    ip.clear();
    x.clear();
    f.clear();
}

/*
 void ChannelFrequencyResponse::calcFreqRes(){

    series->clear();

    //init variables
    double flim=fmax->text().toDouble()+1;
\
    MEAModel<20,float> u1(sig->text().toDouble(),fmax->text().toDouble());
    MEAModel<21,float> u1i(sig->text().toDouble(),fmax->text().toDouble());

    std::vector<float> x;
    std::vector<float> f;
    x.reserve(40);
    f.reserve(40);
    std::vector<float> rp = u1.CalcDiscreteDopplerSpectrum();
    std::vector<float> ip = u1i.CalcDiscreteDopplerSpectrum();

    //calc parametric PSD
    for(int i=0;i < (2*20);i++){

        f[i]=rp[2*i];
        x[i]=rp[1+2*i]+ip[1+2*i];
        series->append(f[i],(x[i]));
    }

    axisX=new QValueAxis;
    axisX->setRange(-flim, flim);
    axisX->setTickCount(10);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleText("t(ms)");
    chartView->chart()->setAxisX(axisX, series);

    axisY=new QValueAxis;
    axisY->setRange(*min_element(x.begin(), x.end())*1.05, *max_element(x.begin(), x.end())*1.05);
    axisY->setTickCount(10);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleText("channel Frequency response");
    chartView->chart()->setAxisY(axisY, series);

    rp.clear();
    ip.clear();
    x.clear();
    f.clear();
}
*/
