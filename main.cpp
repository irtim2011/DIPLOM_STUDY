#include <iostream>
#include <fstream>
#include <math.h>
#define _CRT_SECURE_NO_WARNINGS
constexpr double FR = 0.1; // 2 5
constexpr double KAPPA = (1 / FR); //0.2  0.125 0.1
constexpr double RE = 50.0;//1 10 100
constexpr double EPS = 0.0000001;


using namespace std;

int main()
{
    constexpr int N = 70; /* size of mesh */
    constexpr double delta_x = 1.e-7;
    constexpr double x_dump_step = 0.001;
    constexpr int k_step = int(x_dump_step / delta_x);
    double x = 0.0;
    double L = 10.01;
    int n_steps = L / delta_x;

    double* u_old = new double[N];
    double* u_new = new double[N];

    double* h_old = new double[N];
    double* h_new = new double[N];

    double eq_x = 0.0; //moment of equalization of the speed profile
    int check = 0; //for finding this moment once
    double c;//constant in asimptotic

    /* set starting conditions */
    double q = 1.0;
    for (int i = 0; i < N; i++)
    {
        h_old[i] = i * 1.0 / (N - 1);
        u_old[i] = 1.0 - 1.0 * h_old[i] * h_old[i];
        //u_old[i] = 1.0;
    }

    ofstream out1("start.txt");
    char filename_out_h[1000], filename_out_u[1000];
    sprintf_s(filename_out_h, "h_FR%g_RE%g.csv", FR, RE);
    sprintf_s(filename_out_u, "u_FR%g_RE%g.csv", FR, RE);
    ofstream out_all_h(filename_out_h);
    ofstream out_all_u(filename_out_u);
    /*save-start*/
    for (int i = 0; i < N; i++)
    {

        out1 << u_old[i] << " " << h_old[i] << endl;
    }

    ofstream out2("h_n.txt");
    ofstream out3("h_asimtotic.txt");

    /* calculating part */
    int m = 0;
    for (int k = 0; k < n_steps; k++)
    {
        ///u_n
        for (int i = 1; i < N - 1; i++)
        {
            double u_y = (u_old[i + 1] - u_old[i]) / (h_old[i + 1] - h_old[i]);
            double u_yy = ((u_old[i + 1] - u_old[i]) / (h_old[i + 1] - h_old[i])
                - (u_old[i] - u_old[i - 1]) / (h_old[i] - h_old[i - 1])) / (h_old[i + 1] - h_old[i]);
            u_new[i] = u_old[i] + delta_x * 1.0 / u_old[i] * (1.0 + 1.0 / (KAPPA * RE) * (u_yy + u_y / h_old[i]));
        }
        u_new[N - 1] = u_new[N - 2];
        u_new[0] = u_new[1];

        ///h_n
        h_new[0] = 0.0;
        for (int i = 1; i < N; i++)
        {
            double a = u_new[i - 1] * h_new[i - 1] - u_new[i] * h_new[i - 1];
            double b = u_new[i - 1] * h_new[i - 1] * h_new[i - 1];
            double c = (u_old[i - 1] * h_old[i - 1] + u_old[i] * h_old[i]) * (h_old[i] - h_old[i - 1]);
            h_new[i] = (-a + sqrt(a * a + 4.0 * u_new[i] * (b + c))) / (2.0 * u_new[i]);
        }


        /* rewriting*/
        //upd 17_02_23: + calculate average value
        double average = 0.0; //for calculate an average value
        for (int i = 0; i < N; i++)
        {
            u_old[i] = u_new[i];
            h_old[i] = h_new[i];
            average = average + u_new[i];
        }



        if (k % 1000000 == 0.0)
        {
            //            for(int i = 0; i < N; i++)
            //            {
            //                out2<<h_new[i]<<" ";
            //            }
            //            out2<<x<<endl;
            out2 << h_new[N - 1] << " " << x << endl;
            cout << x << endl;
        }

        if (k % k_step == 0)// dump to the file - all
        {
            if (k == 0) {
                out_all_u << "x;"; out_all_h << "x;";

                for (int j = 0; j < N; j++) {
                    out_all_u << "u_" << j << ";";
                    out_all_h << "h_" << j << ";";
                }
                out_all_u << "\n";
                out_all_h << "\n";
            }
            out_all_u << (delta_x * k) << ";";
            out_all_h << (delta_x * k) << ";";
            for (int j = 0; j < N; j++) {
                out_all_u << u_new[j] << ";";
                out_all_h << h_new[j] << ";";
            }
            out_all_u << "\n";
            out_all_h << "\n";
        }

        ///*we find the moment (double eq_x) of equalization of the speed profile*/
        //for this we define parameter of accuracy (EPS)
        average = average / N;

        double maximum = 0.0; //maximum of accuracy from average

        for (int i = 0; i < N; i++)
        {
            if (abs(average - u_old[i]) > maximum)
            {
                maximum = abs(average - u_old[i]);
            }
        }

        if (maximum < EPS && check == 0)
        {
            eq_x = x;
            check = 1;
            cout << "The profile has leveled off! Moment: " << eq_x << "  RE = " << RE << "  KAPPA = " << KAPPA << endl;

            //now we can calculate a constant in asimptotic
            c = pow(h_old[N - 1], -4.0) * pow(q, -2.0) - 2.0 * eq_x;
            cout << "c = " << c << endl;


            //this end to leveled so write result
            ofstream out4("end.txt");
            for (int i = 0; i < N; i++)
            {
                out4 << u_old[i] << " " << h_old[i] << endl;
            }
            out4.close();
        }

        double h_asimptotic;
        if (check == 1)
        {
            h_asimptotic = pow(q, -0.5) * pow((2.0 * x + c), -0.25);
            if (k % 1000000 == 0.0)
            {
                out3 << h_asimptotic << " " << x << endl;
            }
        }

        x = x + delta_x;

    }
    // cout << N << " " << c << endl;

    out1.close();
    out2.close();
    out3.close();


    return 0;
}
