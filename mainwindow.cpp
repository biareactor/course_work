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

static size_t get_method_order(const QString& method)
{
    size_t order = 0;

    if (method == "Хюна")
        order = 2;
    else if (method == "РК 4")
        order = 4;

    return order;
}

Solver::Params MainWindow::get_params()
{
    Solver::Params params;
    params.x0 = ui->x0->text().toDouble();
    params.u1 = ui->u1->text().toDouble();
    params.u2 = ui->u2->text().toDouble();
    params.h0 = ui->h0->text().toDouble();
    params.border = ui->border->text().toDouble();
    params.eps_border = ui->eps_border->text().toDouble();
    params.n = ui->n->text().toUInt();
    params.a = ui->a->text().toDouble();
    params.alpha = ui->alpha->text().toDouble();
    params.A = ui->A->text().toDouble();
    params.omega = ui->omega->text().toDouble();
    params.step_control = ui->step_control->isChecked();
    params.method = ui->method->currentText();
    //    p.D = ui->D->text().toDouble();

    params.delta = std::pow(10, -10);
    params.eps = std::pow(10, -3);
    params.p = get_method_order(params.method); //method order

    return params;
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

static double modulo(double a)
{
    double r = 2 * M_PI;
    return a - floor (a / r) * r;
}

static std::pair<double, double> find_idx_in_matrix(const std::vector<std::vector<double>>& v2, size_t start)
{
    size_t count = 0;
    size_t i = 0;
    size_t j = 0;

    while (i < v2.size())
    {
        for (j = 0; j < v2[i].size() && count != start; ++j)
            ++count;

        if (count != start)
            ++i;
        else
            break;
    }

    return {i, j};
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

//    if (ui->save->isChecked())
//    {
//        save_csv("./csv_saves", res, ui->file_name->text().toStdString());
//    }

    /*---------------table---------------------*/

//    for (size_t i = 0; i < res.x.size(); i++)
//    {
//        ui->table->insertRow(i);

//        size_t j = 0;
//        for (const auto& value : res.get_values(i))
//        {
//            ui->table->setItem(i, j, new QTableWidgetItem(approx(value)));
//            ++j;
//        }
//    }

    size_t part = ui->part->text().toUInt();
    QString select_border = ui->select_border->currentText();
    size_t start = select_border == "Начало" ? 0 : res.x.size()/part*(part-1);
    size_t end = select_border == "Начало" ? res.x.size()/part : res.x.size();

    /*---------------phase portrait---------------------*/

    auto* chart1 = new QChart();
    auto* series1 = new QLineSeries();
    series1->setColor(QColorConstants::Blue);

    auto[idx1, idx2] = find_idx_in_matrix(res.v2, start);

    for (size_t i = start; i < end; i++)
    {
        series1->append(modulo(res.v1[i]), res.v2[idx1][idx2]);

        if (idx2 + 1 == res.v2[idx1].size())
        {
            chart1->addSeries(series1);
            series1 = new QLineSeries();
            series1->setColor(QColorConstants::Blue);

            ++idx1;
            idx2 = 0;
        }
        else
            ++idx2;
    }

    chart1->createDefaultAxes();

    auto axisX = chart1->axes(Qt::Horizontal).back();
    auto axisY = chart1->axes(Qt::Vertical).back();
    axisX->setRange(0,6.29);
    axisY->setRange(0,8);
    set_axes(axisX, axisY, "x", QChar(0x1E8B), 24);

    chart1->legend()->hide();

    ui->phase_portrait->setRubberBand(QChartView::RectangleRubberBand);
    ui->phase_portrait->setRenderHint(QPainter::Antialiasing);
    ui->phase_portrait->setChart(chart1);

    /*---------------plot of coordinate versus time---------------------*/

//    auto* series2 = new QLineSeries();
//    auto* chart2 = new QChart();

//    for (size_t i = start; i < end; i++)
//    {
//        series2->append(res.x[i], res.v1[i]);
//    }

//    chart2->addSeries(series2);
//    chart2->createDefaultAxes();

//    axisX = chart2->axes(Qt::Horizontal).back();
//    axisY = chart2->axes(Qt::Vertical).back();
//    set_axes(axisX, axisY, "t", "x", 24);

//    chart2->legend()->hide();

//    ui->coord_chart->setRubberBand(QChartView::RectangleRubberBand);
//    ui->coord_chart->setRenderHint(QPainter::Antialiasing);
//    ui->coord_chart->setChart(chart2);

//    /*---------------plot of speed versus time---------------------*/

//    auto* series3 = new QLineSeries();
//    auto* chart3 = new QChart();

//    for (size_t i = start; i < end; i++)
//    {
//        idx1 = 0;
//        idx2 = 0;
//        series3->append(res.x[i], res.v2[idx1][idx2]);
//    }

//    chart3->addSeries(series3);
//    chart3->createDefaultAxes();

//    axisX = chart3->axes(Qt::Horizontal).back();
//    axisY = chart3->axes(Qt::Vertical).back();
//    set_axes(axisX, axisY, "t", QChar(0x1E8B), 24);

//    chart3->legend()->hide();

//    ui->speed_chart->setRubberBand(QChartView::RectangleRubberBand);
//    ui->speed_chart->setRenderHint(QPainter::Antialiasing);
//    ui->speed_chart->setChart(chart3);

//    /*---------------- potential ------------------------*/

//    auto* series4 = new QLineSeries();
//    auto* chart4 = new QChart();

//    for (size_t i = 0; i < res.x.size(); i++)
//    {
//        series4->append(res.x[i], -(p.a * res.x[i] + std::cos(res.x[i])));
//    }

//    chart4->addSeries(series4);
//    chart4->createDefaultAxes();

//    axisX = chart4->axes(Qt::Horizontal).back();
//    axisY = chart4->axes(Qt::Vertical).back();
//    set_axes(axisX, axisY, "t", "x", 24);

//    chart4->legend()->hide();

//    ui->potential->setRubberBand(QChartView::RectangleRubberBand);
//    ui->potential->setRenderHint(QPainter::Antialiasing);
//    ui->potential->setChart(chart4);
}

//static QString make_fullname(const std::string& dir_name, const std::string& file_name)
//{
//    QDirIterator it(dir_name.c_str(), QStringList() << (file_name + "*.csv").c_str(), QDir::Files);

//    size_t total_files = 0;
//    while (it.next() != "")
//        ++total_files;

//    return QString((dir_name + "/" + file_name + (total_files == 0 ? "" : std::to_string(total_files + 1)) + ".csv").c_str());
//}

//static std::string put_comma(bool condition)
//{
//    return condition ? "," : "";
//}

//bool MainWindow::save_csv(const std::string& dir_name, const Solver::Result& res, const std::string& file_name)
//{
//    bool success = false;

//    QFile data(make_fullname(dir_name, file_name));
//    if (data.open(QIODevice::ReadWrite))
//    {
//        QTextStream output(&data);
//        output << res.get_fields_csv() + '\n';
//        for (size_t i = 0; i < res.x.size(); ++i)
//            output << res.get_values_csv(i) + '\n';

//        success = true;
//    }

//    data.close();
//    return success;
//}
