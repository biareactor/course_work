#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include "solver.h"

#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <variant>

#include <QtMath>
#include <QString>
#include <QTableWidgetItem>
#include <QtCharts/QLineSeries>
#include <QValueAxis>
#include <QAbstractAxis>
#include <QFont>
#include <QFileDialog>
#include <QFile>
#include <QDirIterator>
#include <QStringList>

#include <limits>
typedef std::numeric_limits< double > dbl;

static QString approx(double num)
{
    std::ostringstream streamObj;
    streamObj << std::defaultfloat << num;

    return QString::fromStdString(streamObj.str());
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->solve_btn, &QPushButton::clicked, this, &MainWindow::plot_eq);
}

MainWindow::~MainWindow()
{
    delete ui;
}

Solver::Params MainWindow::get_params()
{
    Solver::Params p;
    p.x0 = ui->x0->text().toDouble();
    p.u1 = ui->u1->text().toDouble();
    p.u2 = ui->u2->text().toDouble();
    p.h0 = ui->h0->text().toDouble();
    p.border = ui->border->text().toDouble();
    p.eps_border = ui->eps_border->text().toDouble();
    p.n = ui->n->text().toUInt();
    p.a = ui->a->text().toDouble();
    p.alpha = ui->alpha->text().toDouble();
    p.A = ui->A->text().toDouble();
    p.omega = ui->omega->text().toDouble();
    p.D = ui->D->text().toDouble();

    return p;
}

static void set_axes(QAbstractAxis* axisX, QAbstractAxis* axisY, QString label1, QString label2, size_t fontsize)
{
    auto font = QFont();
    font.setPixelSize(24);

    axisX->setTitleText(label1);
    axisX->setTitleFont(font);

    axisY->setTitleText(label2);
    axisY->setTitleFont(font);
}

void MainWindow::plot_eq()
{
    // Clear table
    while (ui->table->rowCount() > 0)
    {
        ui->table->removeRow(0);
    }

    /*---------------solve by numerical method---------------------*/

    auto p = get_params();
    auto solver = new Solver(p);

    auto res = solver->eq_solve(
    [&p] (double x, double u1, double u2)
    {
        return u2;
    },
    [&p] (double x, double u1, double u2)
    {
        return -p.alpha * u2 + p.a - std::sin(u1) + p.A * std::sin(p.omega * x);
    });

    /*---------------save in file---------------------*/

    if (ui->save->isChecked())
    {
        save_csv("./csv_saves", res, ui->file_name->text().toStdString());
    }

    /*---------------table---------------------*/

    for (size_t i = 0; i < res.x.size(); i++)
    {
        ui->table->insertRow(i);

        size_t j = 0;
        for (const auto& value : res.get_values(i))
        {
            ui->table->setItem(i, j, new QTableWidgetItem(approx(value)));
            ++j;
        }
    }

    /*---------------phase portrait---------------------*/

    auto* series1 = new QLineSeries();
    auto* chart1 = new QChart();

    for (size_t i = res.x.size()/16*15; i < res.x.size(); i++)
    {
        series1->append(res.v[0][i], res.v[1][i]);
    }

    chart1->addSeries(series1);
    chart1->createDefaultAxes();

    auto axisX = chart1->axes(Qt::Horizontal).back();
    auto axisY = chart1->axes(Qt::Vertical).back();
    set_axes(axisX, axisY, "x", QChar(0x1E8B), 24);

    chart1->legend()->hide();

    ui->phase_portrait->setRubberBand(QChartView::RectangleRubberBand);
    ui->phase_portrait->setRenderHint(QPainter::Antialiasing);
    ui->phase_portrait->setChart(chart1);

    /*---------------plot of coordinate versus time---------------------*/

    auto* series2 = new QLineSeries();
    auto* chart2 = new QChart();

    for (size_t i = res.x.size()/16*15; i < res.x.size(); i++)
    {
        series2->append(res.x[i], res.v[0][i]);
    }

    chart2->addSeries(series2);
    chart2->createDefaultAxes();

    axisX = chart2->axes(Qt::Horizontal).back();
    axisY = chart2->axes(Qt::Vertical).back();
    set_axes(axisX, axisY, "t", "x", 24);

    chart2->legend()->hide();

    ui->coord_chart->setRubberBand(QChartView::RectangleRubberBand);
    ui->coord_chart->setRenderHint(QPainter::Antialiasing);
    ui->coord_chart->setChart(chart2);

    /*---------------plot of speed versus time---------------------*/

    auto* series3 = new QLineSeries();
    auto* chart3 = new QChart();

    for (size_t i = res.x.size()/16*15; i < res.x.size(); i++)
    {
        series3->append(res.x[i], res.v[1][i]);
    }

    chart3->addSeries(series3);
    chart3->createDefaultAxes();

    axisX = chart3->axes(Qt::Horizontal).back();
    axisY = chart3->axes(Qt::Vertical).back();
    set_axes(axisX, axisY, "t", QChar(0x1E8B), 24);

    chart3->legend()->hide();

    ui->speed_chart->setRubberBand(QChartView::RectangleRubberBand);
    ui->speed_chart->setRenderHint(QPainter::Antialiasing);
    ui->speed_chart->setChart(chart3);

    /*---------------- potential ------------------------*/

    auto* series4 = new QLineSeries();
    auto* chart4 = new QChart();

    for (size_t i = 0; i < res.x.size(); i++)
    {
        series4->append(res.x[i], -(p.a * res.x[i] + std::cos(res.x[i])));
    }

    chart4->addSeries(series4);
    chart4->createDefaultAxes();

    axisX = chart4->axes(Qt::Horizontal).back();
    axisY = chart4->axes(Qt::Vertical).back();
    set_axes(axisX, axisY, "t", "x", 24);

    chart4->legend()->hide();

    ui->potential->setRubberBand(QChartView::RectangleRubberBand);
    ui->potential->setRenderHint(QPainter::Antialiasing);
    ui->potential->setChart(chart4);
}

static QString make_fullname(const std::string& dir_name, const std::string& file_name)
{
    QDirIterator it(dir_name.c_str(), QStringList() << (file_name + "*.csv").c_str(), QDir::Files);

    size_t total_files = 0;
    while (it.next() != "")
        ++total_files;

    return QString((dir_name + "/" + file_name + (total_files == 0 ? "" : std::to_string(total_files + 1)) + ".csv").c_str());
}

static std::string put_comma(bool condition)
{
    return condition ? "," : "";
}

bool MainWindow::save_csv(const std::string& dir_name, const Solver::Result& res, const std::string& file_name)
{
    bool success = false;

    QFile data(make_fullname(dir_name, file_name));
    if (data.open(QIODevice::ReadWrite))
    {
        QTextStream output(&data);
        output << res.get_fields_csv() + '\n';
        for (size_t i = 0; i < res.x.size(); ++i)
            output << res.get_values_csv(i) + '\n';

        success = true;
    }

    data.close();
    return success;
}
