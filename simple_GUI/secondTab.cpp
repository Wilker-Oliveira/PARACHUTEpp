//Second tab Class Constructor
ChannelAutoCorrelation::ChannelAutoCorrelation(QWidget *parent)
    : QWidget(parent)
{
    //QT init
    series = new QLineSeries();
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
    stepTime->setText(QString::number(0.0001));
    stepTime->setValidator( new QDoubleValidator(0.1, 10, 6, this) );

    //Chart/Chartview setup
    chart->legend()->hide();
    chart->addSeries(series);
    chart->setTitle("I hate to admit it, but this interface is really useful");
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

    series->clear();

    //init variables
    double ti=0,tf=finalTime->text().toFloat(),dt=stepTime->text().toFloat();
    double cohtime=1/fmax->text().toDouble();

    auto *x=new std::vector<float>;
    std::vector<float> t ((tf*cohtime)/dt +1);
    std::vector<float> rp ((tf*cohtime)/dt +1);
    std::vector<float> ip ((tf*cohtime)/dt +1);

    MEAModel<20,float> u1(sig->text().toDouble(),fmax->text().toDouble());
    MEAModel<21,float> u1i(sig->text().toDouble(),fmax->text().toDouble());

    //calc Theoretical Corr
    for(int i=0;i<(tf*cohtime)/dt;i++){
        t[i]=ti + i*dt;
        rp[i]=u1.CalcACF(t[i]);
        ip[i]=u1i.CalcACF(t[i]);
    }
    *x=envACF(sig->text().toDouble(),rp,ip);



    for(int i=0; i<=tf*cohtime/dt; i++){
        series->append(t[i],(x->at(i)));
    }

    axisX=new QValueAxis;
    axisX->setRange(0, tf*cohtime);
    axisX->setTickCount(10);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleText("t(ms)");
    chartView->chart()->setAxisX(axisX, series);

    axisY=new QValueAxis;
    axisY->setRange(*min_element(x->begin(), x->end())*1.05, *max_element(x->begin(), x->end())*1.05);
    axisY->setTickCount(10);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleText("channel response correlation");
    chartView->chart()->setAxisY(axisY, series);

    rp.clear();
    ip.clear();
    x->clear();
    t.clear();
    delete x;
}


