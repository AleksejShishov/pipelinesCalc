#include <iostream>
#include <vector>
#include <chrono>


using namespace std;

const int P0Count = 4; //число входных трубопроводов
const double k = 0.5;  //dummy моск расчета функции гидр. сопротивления 
const int N = 12; // расчетных параметров


void Gauss (vector<vector<double>>&, vector<double>&);

int main()
{
    setlocale(LC_ALL, "");

    vector<vector<double>> a =
    {
        {1, 0, 0, 0, k, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, k, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, k, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, -k, 0, 0, 0},
        {0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, -1, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1},
        {-1, 1, 0, 0, 0, k, 0, 0, 0, 0, 0, 0},
        {0, -1, 0, 1, 0, 0, 0, k, 0, 0, 0, 0},
        {-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, k},
        {0, 0, 0, 1, -1, 0, 0, k, 0, k, 0, 0},
    };

    double p0Const[4]{};
    int i{};
    do
    {
        cout << "Введите входное давление на трубопроводе номер " << i<< ": ";
        cin >> p0Const[i];
        if (p0Const[i] >= 0) 
            i++;
    } while (i<P0Count);

    vector<double> b = { p0Const[0], p0Const[1], p0Const[2], p0Const[3], 0, 0, 0, 0, 0, 0, 0, 0};

    auto start = chrono::high_resolution_clock::now(); // запоминаем время начала выполнения функции
    Gauss(a, b);
    auto end = chrono::high_resolution_clock::now(); // запоминаем время окончания выполнения функции
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start); // вычисляем время выполнения функции в миллисекундах

    cout << "\nВремя расчетов: " << duration.count() << " миллисек." << endl; // выводим время выполнения функции

    return 0;
}

void Gauss(vector<vector<double>>& a, vector<double>& b)
{
    for (int i = 0; i < N; i++)
    {
        int maxRow = i;
        for (int j = i + 1; j < N; j++) ////находим максимальный элемент в столбце i
        {
            if (abs(a[j][i]) > abs(a[maxRow][i]))
            {
                maxRow = j;
            }
        }

        swap(a[i], a[maxRow]); //меняем местами текущую строку и строку с максимальным элементом
        swap(b[i], b[maxRow]);

        for (int j = i + 1; j < N; j++) //приводим матрицу к треугольному виду
        {
            double coef = a[j][i] / a[i][i];
            for (int k = i + 1; k < N; k++)
            {
                a[j][k] -= coef * a[i][k];
            }
            b[j] -= coef * b[i];
            a[j][i] = 0;
        }
    }

    vector<double> x(N);
    char calcValues[][20] = { "P1", "P2", "P3", "P4", "Q01", "Q12", "Q02", "Q24", "Q04", "Q34", "Q03", "Q13" };


    for (int i = N - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < N; j++) {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
    
    cout << endl;
    for (int i = 0; i < N; i++)
    {
        cout << i + 1 << ".\t" << calcValues[i] << " = " << "\t" << x[i] << endl;
    }
} 