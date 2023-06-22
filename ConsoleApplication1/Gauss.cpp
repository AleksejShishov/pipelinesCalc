#include <iostream>
#include <vector>
#include <chrono>


using namespace std;

const double Pi{ 3.14159 };
const double ReCr{ 2400 };

class Liquid
{
public:
    double nu{ 1.004e-6 }; //�������������� �������� ���� ��� 293 �
    double rho{ 1000 };           // ��������� �������� � �����
 };

struct Pipe
{
    Liquid water;

    double L{ 10 };
    double D{ 0.45 };

    double ksi = 64 * L / (ReCr * D);
    double reLam = 16 * Pi * D / 0.03;
    double lam = 8 / (Pi * Pi * D * D * D * D);

public:
    double kTurb = ksi * lam / water.rho;  //������� ������������ - ����. ������������� ��� ������������� ������
    double k = ksi * lam * water.nu * reLam;  //������� ������������ - ����. ������������� ��� ����������� ������
};

struct SparseMatrix {
    vector<double> values; // ��������� �������� �������
    vector<int> rowsCum; // ������� ����� ��������� ���������
    vector<int> columns; // ������� �������� ��������� ���������
};

class HydraulicNet
{
    const int inputsCount; //����� ������� �������������
    const int N = 12; // ����� ��������� ����������
    vector<double> pInput{0, 0, 0, 0}; //������� �������� ��������

    vector<vector<double>> a;
    vector<double> b;
    vector<string> calcValues { "P1", "P2", "P3", "P4", "Q01", "Q12", "Q02", "Q24", "Q04", "Q34", "Q03", "Q13" }; //������� ��� ����������

public:
 
    void InputValues()
    {
        int i{};
        do
        {
            cout << "������� ������� �������� �� ������������ ����� " << i << ": ";
            cin >> pInput[i];
            if (pInput[i] >= 0)
                i++;
        } while (i < inputsCount);
    }

    int FindNonZero()
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

    SparseMatrix ConvertToSparseRow()
    {
        int numNonZeroElements = FindNonZero();
        vector<double> values (numNonZeroElements);
        vector<int> columns (numNonZeroElements);
        vector<int> rowsCum;

        int c = 0;                              //������������ ����� �� 0 ��������� � ������
        for (int i = 0; i < N; i++) {
            rowsCum.push_back(c);
            for (int j = 0; j < N; j++) {
                if (a[i][j] != 0) {
                    values[c] = a[i][j];
                    columns[c] = j;
                    c++;
                }
            }
        }
        rowsCum.push_back(c);

        SparseMatrix matrix;
        matrix.columns = columns;
        matrix.rowsCum = rowsCum;
        matrix.values = values;

        return matrix;
    }

    void ShowSparseRowMatrix(SparseMatrix matrix)
    {
        cout << "\n��������: ";
        for (double value : matrix.values) {
            cout << value << "\t";
        }
        cout << endl;

        cout << "����: \t\t";
        for (int column : matrix.columns) {
            cout << column << " ";
        }
        cout << endl;

        cout << "������ ������: \t";
        for (int rowsCum : matrix.rowsCum) {
            cout << rowsCum << " ";
        }
        cout << endl;
    }
     
    void Gauss()
    {
        for (int i = 0; i < N; i++)
        {
            int maxRow = i;
            for (int j = i + 1; j < N; j++) //������� ������������ ������� � ������� i
            {
                if (abs(a[j][i]) > abs(a[maxRow][i]))
                {
                    maxRow = j;
                }
            }

            swap(a[i], a[maxRow]); //������ ������� ������� ������ � ������ � ������������ ���������
            swap(b[i], b[maxRow]);

            for (int j = i + 1; j < N; j++) //�������� ������� � ������������ ����
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

        vector<double> x(N); //�������������� ������
 
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

    void GaussWithSparseRowMatrix(SparseMatrix a)
    {}

    HydraulicNet(vector<vector<double>>& inputLeftMatrix, int numberOfVertex) : a{inputLeftMatrix}, inputsCount { numberOfVertex }
    {
        InputValues();
        b = { pInput[0], pInput[1], pInput[2], pInput[3], 0, 0, 0, 0, 0, 0, 0, 0 };
    }
};



int main()
{
    setlocale(LC_ALL, "");

    Pipe pipe;
    double k{ pipe.k };

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

    HydraulicNet net4Vertex(a, 4);

    SparseMatrix aSparseRowMatrix = net4Vertex.ConvertToSparseRow();

    net4Vertex.ShowSparseRowMatrix(aSparseRowMatrix);

    auto start = chrono::high_resolution_clock::now(); // ���������� ����� ������ ���������� �������
    net4Vertex.Gauss();
    auto end = chrono::high_resolution_clock::now(); // ���������� ����� ��������� ���������� �������
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start); // ��������� ����� ���������� ������� � �������������

    cout << "����� �������� ������� ������: " << duration.count() << " ��������." << endl; // ������� ����� ���������� �������


    start = chrono::high_resolution_clock::now(); 
    net4Vertex.GaussWithSparseRowMatrix(aSparseRowMatrix);
    end = chrono::high_resolution_clock::now(); 
    duration = chrono::duration_cast<chrono::milliseconds>(end - start); 

    cout << "����� �������� ������ �� ������ ��������: " << duration.count() << " ��������." << endl; // ������� ����� ���������� �������

    return 0;
}

 