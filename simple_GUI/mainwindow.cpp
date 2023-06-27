#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "firstTab.cpp"
#include "secondTab.cpp"
#include "thirdTab.cpp"

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

