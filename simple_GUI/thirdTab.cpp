#include "mainwindow.h"

//Third tab Class Constructor
ChannelFrequencyResponse::ChannelFrequencyResponse(QWidget *parent)
    : QWidget(parent)
{
    //QT init
    Parametric_series = new QScatterSeries();
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
    chart->addSeries(Parametric_series);
    chart->setTitle("Channel Doppler Spectrum.");
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


 void ChannelFrequencyResponse::calcFreqRes(){

    Parametric_series->clear();

    //init variables
    double flim=fmax->text().toDouble()+1;

    MEDModel<20,float> u1(sig->text().toDouble(),fmax->text().toDouble());
    MEDModel<21,float> u1i(sig->text().toDouble(),fmax->text().toDouble());

    std::vector<float> x;
    std::vector<float> f;
    x.reserve(40);
    f.reserve(40);
    std::vector<float> rp = u1.CalcDiscreteDopplerSpectrum();
    std::vector<float> ip = u1i.CalcDiscreteDopplerSpectrum();

    //calc parametric PSD
    for(int i=0;i < 40;i++){

        f[i]=rp[2*i];
        x[i]=rp[1+2*i];
        Parametric_series->append(f[i],(x[i]));
    }

    Parametric_series->setMarkerSize(10);

    axisX=new QValueAxis;
    axisX->setRange(-flim, flim);
    axisX->setTickCount(10);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleText("t(ms)");
    chartView->chart()->setAxisX(axisX, Parametric_series);

    axisY=new QValueAxis;
    //axisY->setRange(0, *max_element(x.begin(), x.end())*1.5);
    axisY->setRange(0,0.06);
    axisY->setTickCount(10);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleText("channel Frequency response");
    chartView->chart()->setAxisY(axisY, Parametric_series);

    rp.clear();
    ip.clear();
    x.clear();
    f.clear();
}
