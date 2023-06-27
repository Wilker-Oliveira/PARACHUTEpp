#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include <QLabel>
#include <numeric>
#include <algorithm>

//Main Class Constructor
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    //Variables
    ui->setupUi(this);

    tabManager = new QTabWidget;
    tabManager->addTab(new ChannelImpulseResponse, tr("Channel response in time"));
    tabManager->addTab(new ChannelAutoCorrelation, tr("Channel autocorrelation"));
    tabManager->addTab(new ChannelFrequencyResponse, tr("Channel doppler spectrum"));


    QWidget *mainWidget = new QWidget(this);
    QVBoxLayout *mainLayout = new QVBoxLayout(mainWidget);
    mainLayout->addWidget(tabManager);

    setCentralWidget(mainWidget);
    setWindowTitle(tr("Simple GUI"));

}
//Main Class Destructor
MainWindow::~MainWindow()
{
    delete ui;
}

//First tab Class Constructor
ChannelImpulseResponse::ChannelImpulseResponse(QWidget *parent)
    : QWidget(parent)
{
    //QT init
    series = new QLineSeries();
    sig = new QLineEdit();
    fmax = new QLineEdit();
    finalTime = new QLineEdit();
    stepTime = new QLineEdit();

    meanLE = new QLineEdit();
    varLE = new QLineEdit();
    meanPLE = new QLineEdit();
    funcType = new QComboBox;
    chart = new QChart();
    chartView = new QChartView(chart);
    m_button = new QPushButton("Refresh series", this);

    //Parameters setup
    QLabel *sigLabel = new QLabel(tr("Standard deviation: "));
    sig->setText(QString::number(1));
    sig->setValidator( new QDoubleValidator(0, 10, 4, this) );

    QLabel *fmaxLabel = new QLabel(tr("Maximun doppler frequency: "));
    fmax->setText(QString::number(91));
    fmax->setValidator( new QDoubleValidator(10, 1000, 4, this) );

    QLabel *finalLabel = new QLabel(tr("Number of coherence time: "));
    finalTime->setText(QString::number(5));
    finalTime->setValidator( new QDoubleValidator(1, 100, 3, this) );

    QLabel *stepLabel = new QLabel(tr("Step time: "));
    stepTime->setText(QString::number(0.0001));
    stepTime->setValidator( new QDoubleValidator(0.1, 10, 6, this) );

    //Statistics setup
    QLabel *meanLabel = new QLabel(tr("Mean of the signal: "));
    meanLE->setReadOnly(true);
    meanLE->setText(QString::number(0));
    QLabel *varLabel = new QLabel(tr("Variance of the signal: "));
    varLE->setReadOnly(true);
    varLE->setText(QString::number(0));
    QLabel *meanPLabel = new QLabel(tr("Mean power of the signal: "));
    meanPLE->setReadOnly(true);
    meanPLE->setText(QString::number(0));

    //ComboBox setup
    QLabel *funcLabel = new QLabel(tr("Choose the model: "));
    funcType->addItem(tr("MED model"));
    funcType->addItem(tr("MEDS model"));

    //Chart/Chartview setup
    chart->legend()->hide();
    chart->addSeries(series);
    chart->setTitle("I'll not kill myself doing this, but i really thought about it");
    chartView->setRenderHint(QPainter::Antialiasing);

    //Refresh button logic
    QObject::connect(m_button, &QPushButton::clicked, this, &ChannelImpulseResponse::recalcSeries);

    ChannelImpulseResponse::recalcSeries();

    //Positioning Setup
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

    QHBoxLayout *mainLayout = new QHBoxLayout;

    mainLayout->addWidget(chartView);
    mainLayout->setStretch(0, 1);
    mainLayout->addWidget(controlWidget);

    setLayout(mainLayout);

}

//First tab Class refresh graph function
void ChannelImpulseResponse::recalcSeries(){
    series->clear();

    auto *x=new std::vector<float>;
    std::vector<float> rp;
    std::vector<float> ip;

    double mean = 0,var = 0, meanP=0,ti=0,tf=finalTime->text().toFloat(),dt=stepTime->text().toFloat();
    double cohtime=1/fmax->text().toDouble();

    std::vector<float> t ((tf*cohtime)/dt +1);

    for(int i=0;i<(tf*cohtime)/dt;i++) t[i]=ti + i*dt;

    MEDModel<20,float> u1(sig->text().toDouble(),fmax->text().toDouble());
    MEDModel<21,float> u1i(sig->text().toDouble(),fmax->text().toDouble());

    MEDSModel<20,float> u2(sig->text().toDouble(),fmax->text().toDouble());
    MEDSModel<21,float> u2i(sig->text().toDouble(),fmax->text().toDouble());

    if (funcType->currentText()==("MED model")){
        rp=u1.CalcProcess(t);
        ip=u1i.CalcProcess(t);
        for(int i=0; i<=tf*cohtime/dt; i++){
            x->push_back(std::abs(std::complex<double>(rp[i],ip[i])));
            series->append(t[i],10*log10(x->at(i)));
        }


    }
    else {
        rp=u2.CalcProcess(t);
        ip=u2i.CalcProcess(t);
        for(int i=0; i<=tf*cohtime/dt; i++){
            x->push_back(std::abs(std::complex<double>(rp[i],ip[i])));
            series->append(t[i],10*log10(x->at(i)));
        }


    }

    for (auto &v : *x)
        mean+=v;
    mean/=x->size();
    for (auto &v : *x)
        var+=pow(v-mean,2);
    var=(var/x->size());

    for (auto &v : *x)
        meanP+=pow(fabs(v),2);
    meanP/=x->size();


    meanLE->setText(QString::number(mean));
    varLE->setText(QString::number(var));
    meanPLE->setText(QString::number(meanP));

    axisX=new QValueAxis;
    axisX->setRange(0, tf*cohtime);
    axisX->setTickCount(10);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleText("t(ms)");
    chartView->chart()->setAxisX(axisX, series);

    axisY=new QValueAxis;
    axisY->setRange(10*log10(*min_element(x->begin(), x->end()))*1.05, 10*log10(*max_element(x->begin(), x->end()))*1.05);
    axisY->setTickCount(10);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleText("channel gain (dB)");
    chartView->chart()->setAxisY(axisY, series);

    rp.clear();
    ip.clear();
    x->clear();
    t.clear();
    delete x;
}

//Useful utility function for calculate the ACF
std::vector<float> envACF(float sig, std::vector<float> acf1, std::vector<float> acf2){
    std::vector<float> en_acf;
    float aux;

    for(float i = 0; i< acf1.size(); i++){
        aux = acf1[i] + acf2[i];
        en_acf.push_back(aux);
    }
    return en_acf;
}

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

    QLabel *fmaxLabel = new QLabel(tr("Maximun doppler frequency: "));
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

    QLabel *fmaxLabel = new QLabel(tr("Maximun doppler frequency: "));
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
