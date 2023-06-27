#include <iostream>
#include <vector>
#include <chrono>


using namespace std;

const double Pi{ 3.14159 };
const double ReCr{ 2400 };

//физическая модель жидкости
struct Liquid
{
    double nu{ 1.004e-6 };  //кинематическая вязкость воды при 293 К
    double rho{ 1000 };     // Плотность жидкости в трубе
 };
//гидравлическая модель трубопровода
struct Pipe
{
    Liquid water;           //Физические параметры воды

    double L{ 10 };         //длина лебра трубы
    double D{ 0.45 };       //диаметр трубы

    double ksi = 64 * L / (ReCr * D);
    double reLam = 16 * Pi * D / 0.03;
    double lam = 8 / (Pi * Pi * D * D * D * D);

public:
    double kTurb = ksi * lam / water.rho;       //функция проводимости - гидр. сопротивления для турбулентного режима
    double k= ksi * lam * water.nu * reLam;    //функция проводимости - гидр. сопротивления для ламинарного режима
};

//структура для хранения сконвертированной (сжатой) разреженной матрицы
struct ConvertedMatrix
{
    vector<double> values;  // ненулевые элементы матрицы
    vector<int> rows;       // индексы строк ненулевых элементов
    vector<int> columns;    // индексы столбцов ненулевых элементов
};


//класс для описания гидравлической сети
class HydraulicNet
{
    const int pInputCount{ 4 };         //число входных трубопроводов
    const int N{ 12 };                  //число расчетных параметров

    vector<vector<double>> a;
    vector<double> b;
    vector<double> pInput{0,0,0,0};           //входные значения давлений
    vector<string> calcValues { "P1", "P2", "P3", "P4", "Q01", "Q12", "Q02", "Q24", "Q04", "Q34", "Q03", "Q13" }; //подписи для результата

    //возвращает число не нулевых элементов в матрице а
    int CountNonZeroElements()
    {
        int numNonZeroElements{};
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (a[i][j] != 0) {
                    numNonZeroElements++;
                }
            }
        }
        return numNonZeroElements;
    }

public:
    ConvertedMatrix sparseRowA;
    ConvertedMatrix coordinateVectorsA;

    //ввод начальных давлений
    void InputValues()
    {
        int i{};
        do
        {
            cout << "Введите входное давление на трубопроводе номер " << i << " в Паскалях: ";
            cin >> pInput[i];
            if (pInput[i] >= 0)
                i++;
        } while (i < pInputCount);
    }

     //конвертирует разреженную матрицу в формат Sparse Row
    ConvertedMatrix ConvertToSparseRow()
    {
        int numNonZeroElements = CountNonZeroElements();
        vector<double>values (numNonZeroElements);
        vector<int> columns (numNonZeroElements);
        vector<int> rowsCumulative;
        ConvertedMatrix sparseRow;

        int c = 0;                              //кумулятивное число не нулевых элементов в строке
        for (int i = 0; i < N; i++)
        {
            rowsCumulative.push_back(c);
            for (int j = 0; j < N; j++) 
            {
                if (a[i][j] != 0)
                {
                    values[c] = a[i][j];
                    columns[c] = j;
                    c++;
                }
            }
        }
        rowsCumulative.push_back(c);

        sparseRow.columns = columns;
        sparseRow.rows = rowsCumulative;
        sparseRow.values = values;
        return sparseRow;
    }

    //конвертирует разреженную матрицу в формат 3-х векторов: значения, координаты i, координаты j
    ConvertedMatrix ConvertToCoordinateVectors()
    {
        int k{}, numNonZeroElements = CountNonZeroElements();
        vector<double>values(numNonZeroElements);
        vector<int> columns(numNonZeroElements);
        vector<int> rows (numNonZeroElements);
        ConvertedMatrix coordinateVectors;

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (a[i][j] != 0)
                {
                    values[k] = a[i][j];
                    rows[k] = i;
                    columns[k] = j;
                    k++;
                }
            }
        }
 
        coordinateVectors.columns = columns;
        coordinateVectors.rows = rows;
        coordinateVectors.values = values;
        return coordinateVectors;
    }


    //выводит на консоль матрицу в формате Sparse Row
    void ShowSparseRowMatrix(const ConvertedMatrix matrix)
    {
        cout << "\nЗначения: ";
        for (double value : matrix.values) {
            cout << value << "\t";
        }
        cout << endl;

        cout << "Ряды: ";
        for (int column : matrix.columns) {
            cout << column << " ";
        }
        cout << endl;

        cout << "Строки: ";
        for (int rowsCum : matrix.rows) {
            cout << rowsCum << " ";
        }
        cout << endl;
    }
     
    //считает матрицу a и вектор b по методу Гаусса
    void Gauss()
    {
        for (int i = 0; i < N; i++)
        {
            int maxRow = i;
            for (int j = i + 1; j < N; j++) //находим максимальный элемент в столбце i
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

        vector<double> x(N); //результирующий вектор
 
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
            cout << i + 1 << ".\t" << calcValues[i] << " = " << "\t" << x[i] <<endl;
        }
    }
    
    //считает матрицу Sparse Row и вектор b по методу Гауса 
    void GaussWithSparseRowMatrix()
    {}

    //коструктор
    HydraulicNet(vector<vector<double>>& inputAMatrix, int numberOfVertex) : a{inputAMatrix}, pInputCount { numberOfVertex }
    {
        InputValues();
        b = { pInput[0], pInput[1], pInput[2], pInput[3], 0, 0, 0, 0, 0, 0, 0, 0 };
    }
};

//прототип; вычисляет время выполнения ф-ции, в качестве параметров передатся объект и его метод, время выполнения которого измеряется
double GetCalculationRuntime(HydraulicNet, void (HydraulicNet::* calculation)());

int main()
{
    setlocale(LC_ALL, "");

    Pipe pipe;
    double k{ pipe.k };                         //для системы линейных уравнений сразу вычисляем коэфф. проводимости - гидравлическое сопротивление

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

    HydraulicNet net(a, 4);                                                  //гидравлическая сеть с 4-я узлами

    net.sparseRowA = net.ConvertToSparseRow();                               //создаем и выводим матрицу со сжатыми рядами        
    net.ShowSparseRowMatrix(net.sparseRowA);

    net.coordinateVectorsA = net.ConvertToCoordinateVectors();               //создаем и выводим векторы значений и координат 
    net.ShowSparseRowMatrix(net.coordinateVectorsA);

    double funcRuntime = GetCalculationRuntime(net, &HydraulicNet::Gauss);   // выводим время выполнения функции
    cout << "\nВремя расчетов простым Гаусом, мс: " << funcRuntime << endl;

    funcRuntime = GetCalculationRuntime(net, &HydraulicNet::GaussWithSparseRowMatrix);    // выводим время выполнения функции
    cout << "\nВремя расчетов  Гаусом c матрицей SparseRow, мс: " << funcRuntime << endl;

    return 0;
} 

double GetCalculationRuntime(HydraulicNet net, void (HydraulicNet::* calculation)())
{
    auto start = chrono::high_resolution_clock::now();                                    // запоминаем время начала выполнения функции

    (net.*calculation)();

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);             // вычисляем время выполнения функции в миллисекундах
    return duration.count();
}
