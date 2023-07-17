/** Window and tab classes. */
#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "./firstTab.cpp"
#include "./secondTab.cpp"
#include "./thirdTab.cpp"

#include <QLabel>
#include <numeric>
#include <algorithm>

/** Main Class Constructor. */
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{

    ui->setupUi(this);

    tabManager = new QTabWidget;                                                        /** Tab object. */
    tabManager->addTab(new ChannelTimeResponse, tr("Channel response in time"));        /** First tab: Channel response in time. */
    tabManager->addTab(new ChannelAutoCorrelation, tr("Channel autocorrelation"));      /** Second tab: Channel autocorrelation. */
    tabManager->addTab(new ChannelFrequencyResponse, tr("Channel doppler spectrum"));   /** Third tab: Channel doppler spectrum. */


    QWidget *mainWidget = new QWidget(this);                 /** Main widget object. */
    QVBoxLayout *mainLayout = new QVBoxLayout(mainWidget);   /** Layout of the widget object. */
    mainLayout->addWidget(tabManager);

    setCentralWidget(mainWidget);      /** Making central widget at the window. */
    setWindowTitle(tr("Simple GUI"));  /** Title of the window. */

}
/** Main Class Destructor. */
MainWindow::~MainWindow()
{
    delete ui;
}

