#include <Mahi/Gui.hpp>
#include <Mahi/Util.hpp>
#include <Eigen/Core> // для Vector2d
#include <cmath>

using namespace mahi::gui;
using namespace mahi::util;
using namespace Eigen;


enum CurveType {
    Normal, LimitingCycle
};

// кривая - в том виде, в котором её примет ImPlot
struct Curve {
    std::vector<double> xs;
    std::vector<double> ys;
    CurveType type = Normal;
};



// само приложение
class CauchyVis : public Application {
public:
    std::vector<Curve> curves;

    // правая часть нашей задачи
    double alpha = 1;

    Vector2d f(Vector2d p) {
        double x = p.x();
        double y = p.y();

        return Vector2d( x + y - alpha*x*x*x,
                        -x + y - y*y*y);
    }



    // метод Рунге-Кутты  5го порядка (формулы из Бахвалова)
    Vector2d rungekutta(Vector2d p, double& h, double eps) {
        while(1) {
            Vector2d k1 = h*f(p);
            Vector2d k2 = h*f(p + k1/2);
            Vector2d k3 = h*f(p + k1/4 + k2/4);
            Vector2d k4 = h*f(p - k2 + 2*k3);
            Vector2d k5 = h*f(p + (7*k1 + 10*k2 + k4)/27);
            Vector2d k6 = h*f(p + (28*k1 - 125*k2 + 546*k3 + 54*k4 - 378*k5));
            Vector2d dp = (k1 + 4*k3 + k4)/6;
            Vector2d delta_ = -(42*k1 + 224*k3 + 21*k4 - 162*k5 - 125*k6)/336;
            double delta = delta_.norm();

            if (delta < eps) {
                //printf("%f\n", h);
                return dp;
            }

            double chi = pow(delta/eps, 1.0/6);
            chi = std::fmin(chi, 10);
            chi = std::fmax(chi, 0.1);
            h = 0.95*h/chi;
        }
    }



    double eps = 0.01;

    void new_curve(double x, double y, double dist) {
        Vector2d p = {x, y};
        Curve curve = {{p.x()}, {p.y()}};

        double h = 1;
        if (dist < 0) {
            h = -1;
            dist = -dist;
        }

        while (dist > 0) {
            Vector2d dp = rungekutta(p, h, eps);
            dist -= std::fmax(std::fabs(h), dp.norm());

            p += dp;
            curve.xs.push_back(p.x());
            curve.ys.push_back(p.y());

            //printf("step %f\t%f\t%f\n", dist, p.x(), p.y());
        }

        curves.push_back(curve);
    }



    double cycle_y  = 0;
    double cycle_x1 = 1;
    double cycle_x2 = 2;

    // обход по "циклу" один раз - т.е. пока горизонтальная прямая не пересечется ещё раз
    // возвращает сдвиг по оси x
    // кривую создает по указателю, если передан не NULL
    double search_cycle_loop(double x, double y, double h, Curve* result) {
        Vector2d p = {x, y};

        if (result) {
            result->xs = {x};
            result->ys = {y};
            result->type = LimitingCycle;
        }

        // первый шаг: чтобы понять, идем вверх или вниз
        Vector2d first_dp = rungekutta(p, h, eps);
        bool up = first_dp.y() > 0;

        p = p + first_dp;
        while(1) {
            Vector2d p_orig = p;
            Vector2d dp = rungekutta(p, h, eps);
            p += dp;

            bool crossing = (up && p_orig.y() < y && p.y() > y) || (!up && p_orig.y() > y && p.y() < y);
            if (crossing) {
                if (std::abs(dp.x()) < eps) {
                    // считаем, что пересечение здесь
                    printf("dx: %f\n", p.x() - x);
                    return p.x() - x;
                } else {
                    // откатываемся на шаг назад
                    // и уменьшаем шаг вдвое
                    p = p_orig;
                    h = h/2;
                    printf("back: %f\n", h);
                }
            }

            if (result) {
                result->xs.push_back(p.x());
                result->ys.push_back(p.y());
            }
        }
    }

    // поиск предельного цикла с помощью предыдущей функции делением пополам
    // (предполагается, что цикл устойчивый или неустойчивый)
    void search_cycle(double h) {
        double x1 = cycle_x1;
        double x2 = cycle_x2;

        double dx1 = search_cycle_loop(x1, cycle_y, h, NULL);
        double dx2 = search_cycle_loop(x2, cycle_y, h, NULL);

        double final_x;

        while(1) {
            double x3 = (x1 + x2) / 2;
            double dx3 = search_cycle_loop(x3, cycle_y, h, NULL);
            printf("cycle x: %f\n", x3);

            if (std::abs(dx3) < eps) {
                printf("final cycle x: %f\n", x3);

                final_x = x3;
                break;
            } else {
                if (std::signbit(dx3) == std::signbit(dx1)) {
                    x1 = x3;
                    dx1 = dx3;
                } else {
                    x2 = x3;
                    dx2 = dx3;
                }
            }
        }

        Curve curve;
        search_cycle_loop(final_x, cycle_y, h, &curve);
        curves.push_back(curve);
    }



    CauchyVis() : Application(1024, 720,"CauchyVis") {}

    // тут интерфейс
    // рисуются кривые
    void update() override {
        ImGuiViewport* main_viewport = ImGui::GetMainViewport();
        ImVec2 workPos = main_viewport->GetWorkPos();

        ImGui::SetNextWindowPos(ImVec2(workPos.x + 750, workPos.y + 10), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(250,500), ImGuiCond_Once);
        ImGui::Begin("Options");

        ImGui::TextWrapped("Ctrl+Click creates an integral curve. Left click - positive direction, right click - negative direction");
        ImGui::TextWrapped("Ctrl+Shift+Click for both directions");

        if (ImGui::Button("Clear")) {
            curves.clear();
        }

        static double curveLength = 2;
        ImGui::InputDouble("length", &curveLength, 0.5, 2);

        static int qfm = 10;
        static int qfn = 10;
        bool quickfill = ImGui::Button("Quick fill");
        ImGui::InputInt("m", &qfm);
        ImGui::InputInt("n", &qfn);

        static bool showStartPoints = true;
        ImGui::Checkbox("show start points", &showStartPoints);

        ImGui::Separator();
        ImGui::Text("Approximation options:");
        ImGui::InputDouble("eps", &eps, 0.0001, 0.01, "%.6f");

        ImGui::Separator();
        ImGui::Text("Problem parameters:");
        ImGui::InputDouble("alpha", &alpha, 0.01, 1, "%.3f");

        ImGui::Separator();
        static bool showCycleSearchControls = false;
        if (showCycleSearchControls = ImGui::CollapsingHeader("Limiting cycle search")) {
            if (ImGui::Button("Search forward")) {
                search_cycle(1);
            }
            if (ImGui::Button("Search backward")) {
                search_cycle(-1);
            }
        }

        ImGui::End();


        ImGui::SetNextWindowPos(ImVec2(workPos.x + 10, workPos.y + 10), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(700, 700), ImGuiCond_Once);
        ImGui::Begin("Plots");

        bool ctrl = ImGui::GetIO().KeyCtrl;
        bool shift = ImGui::GetIO().KeyShift;
        int flags = ctrl ? ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect : 0;
        int axisflags = ctrl ? ImPlotAxisFlags_Lock : 0;

        if (ImPlot::BeginPlot("Plot", "x", "y", ImVec2(-1, -1), flags, axisflags, axisflags)) {
            // Ctrl+Left/Right Click - выпустить кривую

            if (ctrl && ImPlot::IsPlotHovered()) {
                ImPlotPoint pt = ImPlot::GetPlotMousePos();

                if (ImGui::IsMouseClicked(0) || (ImGui::IsMouseClicked(1) && shift)) {
                    new_curve(pt.x, pt.y, curveLength);
                }

                if (ImGui::IsMouseClicked(1) || (ImGui::IsMouseClicked(0) && shift)) {
                    new_curve(pt.x, pt.y, -curveLength);
                }
            }

            // quick fill

            if (quickfill) {
                ImPlotLimits limits = ImPlot::GetPlotLimits();
                double w = limits.X.Size();
                double h = limits.Y.Size();

                for (int i = 0; i <= qfm; i++) {
                    for (int j = 0; j <= qfn; j++) {
                        double x = limits.X.Min + i*w/qfm;
                        double y = limits.Y.Min + j*h/qfn;

                        new_curve(x, y, curveLength);
                        new_curve(x, y, -curveLength);
                    }
                }
            }

            // рисуем кривые

            for (Curve &curve : curves) {
                const char* desc = curve.type == LimitingCycle ? "limiting cycle" : "integral curve";

                ImPlot::PlotLine(desc, &curve.xs[0], &curve.ys[0], curve.xs.size());
                if (showStartPoints) {
                    ImPlot::PlotScatter(desc, &curve.xs[0], &curve.ys[0], 1);
                }
            }

            // управление для поиска предельного цикла

            if (showCycleSearchControls) {
                ImPlot::DragLineY("cycle search y", &cycle_y);
                ImPlot::DragLineX("cycle search x1", &cycle_x1);
                ImPlot::DragLineX("cycle search x2", &cycle_x2);
            }

            ImPlot::EndPlot();
        }

        ImGui::End();
    }
};



int main() {
    CauchyVis app;
    app.run();
    return 0;
}
