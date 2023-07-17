#include "mainwindow.h"

/** First tab Class Constructor. */
ChannelTimeResponse::ChannelTimeResponse(QWidget *parent)
    : QWidget(parent)
{
    /** Parameters of the process. */
    series = new QLineSeries();                          /** Stores the channel gains. */
    sig = new QLineEdit();                               /** Standard deviation of the process. */
    fmax = new QLineEdit();                              /** Maximum doppler frequency. */
    finalTime = new QLineEdit();                         /** Number of coherence time. */
    stepTime = new QLineEdit();                          /** Time step. */

    /** Statistics of the process. */
    meanLE = new QLineEdit();                            /** Mean of the channel. */
    varLE = new QLineEdit();                             /** Variance of the channel. */
    meanPLE = new QLineEdit();                           /** Mean of the process. */

    /** Objects of the widget. */
    funcType = new QComboBox;                            /***/
    chart = new QChart();                                /** Object to plot the series. */
    chartView = new QChartView(chart);                   /** Object to visualize the plot of the series. */
    m_button = new QPushButton("Refresh series", this);  /** Button to apply the parameters and compute the channel. */

    //Parameters setup
    QLabel *sigLabel = new QLabel(tr("Standard deviation: "));
    sig->setText(QString::number(1));
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
    chart->setTitle("Channel Time Response.");
    chartView->setRenderHint(QPainter::Antialiasing);

    //Refresh button logic
    QObject::connect(m_button, &QPushButton::clicked, this, &ChannelTimeResponse::recalcSeries);

    ChannelTimeResponse::recalcSeries();

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
void ChannelTimeResponse::recalcSeries(){
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
