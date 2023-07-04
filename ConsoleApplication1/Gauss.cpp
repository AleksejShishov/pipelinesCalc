#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>

#include "gauss.h"

using namespace std;

const string line(100, '=');                    //строка разделитель для вывода на консоль
const string filenameM = "matrix.txt";          //матрицу можно загрузить из файла
const string filenameResultJSON = "result.json";//расчетные параметры сохраняются в формате JSON 
const vector<vector<double>> ERROR = { {1.7e308}, {0} };


vector<string> signatureX { "P1", "P2", "P3", "P4", "Q01", "Q12", "Q02", "Q24", "Q04", "Q34", "Q03", "Q13" }; //расчетные параметры, подписи
vector<string> signatureXN { "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34", "Q01", "Q02", "Q03", "Q04" }; //расчетные параметры, подписи
//double Q01{ 0 }, Q12{ 0 }, Q02{ 0 }, Q24{ 0 }, Q04{ 0 }, Q34{ 0 }, Q03{ 0 }, Q13{ 0 }, p1, p2, p3, p4;

const int PInputCount{ 4 };                     //число входных трубопроводов, 4 по умолчанию
const int N{ 12 };                              //число расчетных параметров, 12 по умолчанию, передаётся в конструкторе в модель

//режимы для отладки
enum class Mode
{
    Debug = false,
    Release = true
};

//
// ФИЗИЧЕСКАЯ МОДЕЛЬ
// 
const double Pi{ 3.14159 };
const double ReCr{ 2400 };

//физическая модель жидкости
struct Liquid
{
    double nu{ 1.004e-6 };                      //кинематическая вязкость воды при 293 К
    double rho{ 1000 };                         // Плотность жидкости в трубе
 };

//гидравлическая модель трубопровода
struct Pipe
{
    Liquid water;                               //Физические параметры воды

    double L{ 10 };                             //длина лебра трубы
    double D{ 0.45 };                           //диаметр трубы

    double ksi = 64 * L / (ReCr * D);
    double reLam = 16 * Pi * D / 0.03;
    double lam = 8 / (Pi * Pi * D * D * D * D);

public:
    double kTurb = ksi * lam / water.rho;      //функция проводимости - гидр. сопротивления для турбулентного режима
    double k= ksi * lam * water.nu * reLam;    //функция проводимости - гидр. сопротивления для ламинарного режима
};

//
// СТРУКТУРЫ ДАННЫХ ДЛЯ ХРАНЕНИЯ РАЗРЕЖЕННОЙ МАТРИЦЫ
// 
enum class DataType
{
    CLA,
    CSRA,
    CSRAA,
    KNUTHRowsList,
    KNUTHColumnsList
};

//матрица в виде 3-ёх векторов: значения, координаты (2 варианта использования: просто вектор i - CLA / вектор i сжат - CSRA)
struct CompressedMatrix
{
    vector<double> values;                      //ненулевые элементы матрицы
    vector<int> rows;                           //индексы строк ненулевых элементов
    vector<int> columns;                        //индексы столбцов ненулевых элементов
};

//матрица в виде 4 векторов: значений, координат, суммы элементов в строке (CSR matrix + vector j)
struct CompressedMatrix4V : CompressedMatrix
{
   vector<int> nextRow;                         //указатели на номер элемента - начала строки
};

//матрица в виде связных списков по схеме Кнута
struct Node {
    double value;
    int column;
    int row;
    Node* nextInRow;
    Node* nextInColumn;
};

//
// МАТЕМАТИЧЕСКАЯ МОДЕЛЬ
//
//класс для описания гидравлической сети. 
class HydraulicNet
{
    const int N{ 12 };                          //число расчетных параметров, 12 по умолчанию, передаётся в конструкторе
    const double k{ 0.5 };                      //коэфициент проводимости, гидравлическое сопротивление, передаётся в кострукторе

    vector<vector<double>> a;                   //левая часть матрицы, представляющей СЛАУ
 
    int numNonZero;
    int CountNonZeroElements(vector<vector<double>>& matrix)
    {
        int numNonZeroElements{ 0 };
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (matrix[i][j] != 0) {
                    numNonZeroElements++;
                }
            }
        }
        return numNonZeroElements;
    }

public:
    CompressedMatrix CLA;                       //Coordinate List matrix A, 3 вектора не 0 элементов: значенния, координаты ряда, столбца
    CompressedMatrix СSRA;                      //Compressed Sparse Row matrix A, , вектор координат рядов сжат
    CompressedMatrix4V CSRAA;                   //формат 4-х векторов: значения, координаты i, координаты j, суммы элементов в строке, комбинация CLA и CRSA
    vector<Node*> KNUTHrows;                    //Связные списки по модели Кнута, добавляются векторы входа в новый ряд и новый столбец
    vector<Node*> KNUTHcolumns;

    vector<string> signatureX { "P1", "P2", "P3", "P4", "Q01", "Q12", "Q02", "Q24", "Q04", "Q34", "Q03", "Q13" }; //расчетные параметры, подписи для Гаусса
    vector<string> signatureXN { "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34", "Q01", "Q02", "Q03", "Q04" }; //расчетные параметры, подписи для Ньютона
 //   double Q01{ 0 }, Q12{ 0 }, Q02{ 0 }, Q24{ 0 }, Q04{ 0 }, Q34{ 0 }, Q03{ 0 }, Q13{ 0 };

    //
    //********************************* КОНВЕРТАЦИЯ РАЗРЕЖЕННОЙ МАТРИЦЫ *************
    //
    //конвертирует разреженную матрицу в формат 3-х векторов: значения, координаты i, координаты j
    CompressedMatrix ConvertToCLA(const vector<vector<double>>& matrix)
    {
        vector<double>values;
        vector<int> columns;
        vector<int> rows;
        CompressedMatrix coordinateVectors;

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (matrix[i][j] != 0)
                {
                    values.push_back(matrix[i][j]);
                    rows.push_back(i);
                    columns.push_back(j);
                }
            }
        }

        coordinateVectors.columns = columns;
        coordinateVectors.rows = rows;
        coordinateVectors.values = values;
        return coordinateVectors;
    }

    //конвертирует разреженную матрицу в формат Sparse Row
    CompressedMatrix ConvertToCSRA(const vector<vector<double>>& matrix)
    {
        vector<double>values;
        vector<int> columns;
        vector<int> rowsCumulative;        
        CompressedMatrix CSRA;

        int c = 0;                              //кумулятивное число не нулевых элементов в строке
        rowsCumulative.push_back(0);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++) 
            {
                if (matrix[i][j] != 0)
                {
                    values.push_back(matrix[i][j]);
                    columns.push_back(j);
                    c++;
                }
            }
        rowsCumulative.push_back(c);
        }

        CSRA.columns = columns;
        CSRA.rows = rowsCumulative;
        CSRA.values = values;
        return CSRA;
    }

    //конвертирует разреженную матрицу в формат 4-х векторов: значения, координаты i, координаты j, суммы элементов в строке
    CompressedMatrix4V ConvertToCSRAAdapted(const vector<vector<double>>& matrix)
    {
        vector<double> values;
        vector<int> columns;
        vector<int> rows;
        vector<int> nextRow;
        CompressedMatrix4V CSRAA;

        int c = 0;
        nextRow.push_back(0);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (matrix[i][j] != 0)
                {
                    values.push_back(matrix[i][j]);
                    columns.push_back(j);
                    rows.push_back(i);
                    c++;
                }
            }
            nextRow.push_back(c);
        }

        CSRAA.columns = columns;
        CSRAA.rows = rows;
        CSRAA.nextRow = nextRow;
        CSRAA.values = values;
        return CSRAA;
    }

    //конвертирует разреженную матрицу в формат Кнута - векторы начала столбцов, начала строк, СИЛЬНО СВЯЗНЫЙ СПИСОК
    void ConvertToKnuthLists(const vector<vector<double>>& matrix) 
    { 
        KNUTHrows.resize(N, nullptr);
        KNUTHcolumns.resize(N, nullptr);

        for (int i = 0; i < N; i++) 
        {
            for (int j = 0; j < N; j++) 
            {
                double value = matrix[i][j];
                if (value != 0)
                {
                    Node* element = new Node();
                    element->row = i;
                    element->column = j;
                    element->value = value;

                    if (KNUTHrows[i] == nullptr)                             //вход в новую строку
                    {
                        KNUTHrows[i] = element;
                    }
                    else
                    {
                        Node* currentRow = KNUTHrows[i];
                        while (currentRow->nextInRow != nullptr)             //перемещаемся в конец строки, чтобы добавит новый не начаьльный элемент
                        {
                            currentRow = currentRow->nextInRow;
                        }
                        currentRow->nextInRow = element;
                    }

                    if (KNUTHcolumns[j] == nullptr)
                    {
                        KNUTHcolumns[j] = element;
                    }
                    else 
                    {
                        Node* currentColumn = KNUTHcolumns[j];
                        while (currentColumn->nextInColumn != nullptr)
                        {
                            currentColumn = currentColumn->nextInColumn;
                        }
                        currentColumn->nextInColumn = element;
                    }
                }
            }
        }
    }


    //
    //********************************* ВЫЧИСЛИТЕЛЬНЫЙ БЛОК **************************
    // 
    //считает матрицу a и вектор b по методу Гаусса    
    vector<double> Gauss(vector<double>& b)
    {
        vector<double> x(N, 0.0);               //собираем результат в вектор x

        for (int i = 0; i < N; i++)
        {
            int maxRow{ i }; 
            for (int j = i + 1; j < N; j++)     //находим максимальный элемент в столбце i, алгоритм с главным элементом
            {
                if (abs(a[j][i]) > abs(a[maxRow][i]))
                {
                    maxRow = j;
                }
            }

            if (a[maxRow][i] == 0.0)            //нулевой столбец
            {
                continue;
            }

            if (maxRow != i)
            {
                swap(a[i], a[maxRow]);          //меняем местами текущую строку и строку с максимальным элементом
                swap(b[i], b[maxRow]);
            }

            for (int j = i + 1; j < N; j++)     //прямой ход 
            {
                if (a[j][i] != 0)
                {
                    double coef = a[j][i] / a[i][i];
                    for (int k = i + 1; k < N; k++)
                    {
                        a[j][k] -= coef * a[i][k];
                    }
                    b[j] -= coef * b[i];
                    a[j][i] = 0;               //обнуляем параметры вручную, это: -1 итерация и точный 0,  вместо Aji - Aji*Aii/Aii
                }
            }
        }
 
        for (int i = N - 1; i >= 0; i--)       //обратный ход 
        {     
            double sum{ 0.0 };
            for (int j = i + 1; j < N; j++)
            {
                sum += a[i][j] * x[j];          //подставляем найденные х в строку выше
            }
            x[i] = (b[i] - sum) / a[i][i];      //результирующий вектор 
        }
        return x;
    }

    //считает матрицу в виде Coordinate List и вектор b по методу Гауса 
    vector<double> GaussWithCoordinateVectors(vector<double>& b)
     {
         vector<double> x(N, 0.0);                              //собираем результат в вектор x

        for (int i = 0; i < N; i++)
        {
            double mainElemAbs{ 0.0 }, mainElem{ 0.0 };
            int mainRow{ 0 };                                   //ряд главного элемента
            int mainIndex{ 0 };                                 //индекс главного элемента
            for (int j = 0; j < numNonZero; j++)                //находим максимальный элемент в столбце i, алгоритм с главным элементом
            {
                if (CLA.columns[j] == i && abs(CLA.values[j]) > mainElemAbs && CLA.rows[j] >= i)
                {
                    mainElemAbs = abs(CLA.values[j]);
                    mainElem = CLA.values[j];
                    mainRow = CLA.rows[j];
                    mainIndex = j;
                }
            }

            if (mainElem == 0.0)                                //нулевой столбец
            {
                continue;
            }
 
            if (mainRow != i)
            {
                for (int j = i; j < numNonZero;j++)              //меняем местами текущую строку и строку с максимальным элементом
                {
                    if (CLA.rows[j] == i)
                    {
                        CLA.rows[j] = mainRow;
                        continue;
                    }
                    if (CLA.rows[j] == mainRow)
                    {
                        CLA.rows[j] = i;
                    }
                }
                double temp{ b[i] };                            //меняем местами вектор свободных членов                                
                b[i] = b[mainRow];
                b[mainRow] = temp;
            }  

            bool zeroAtCurrentRow;                              //признак, что в строке есть 0 элемент, который после вычитания строе !=0
            vector<double> newElements;                         //элементы, которые  перестают быть 0 после вычитания строк
            vector<int> newElColumns;                           //номер столбца таких элементов

            for (int j = 0; j < numNonZero; j++)                //прямой ход 
            {
                if ((CLA.columns[j] == i) && (CLA.rows[j] > i))
                {
                    double coef = CLA.values[j] / mainElem;
                    int currentRow = CLA.rows[j];
                    for (int l = mainIndex; CLA.rows[l] == mainRow; l++)             //строка mainRow с текущим главным элементом
                    {
                        zeroAtCurrentRow = true;
                        if (CLA.rows[l] == i  && CLA.columns[l] > i)
                        {
                            for (int k = j; CLA.rows[k] == currentRow; k++)          //строка workRow для обнуления и вычитания
                            {
                                if (CLA.columns[k] == CLA.columns[l])                //не 0 элементы в одном столбце mainRow и workRow
                                {
                                    CLA.values[k] -= coef * CLA.values[l];
                                    zeroAtCurrentRow = false;
                                }
                                if (k + 1 == numNonZero)
                                {
                                    break;
                                }
                            }
                            if (zeroAtCurrentRow)                                    //не 0 элемент только в mainRow
                            {
                                newElements.push_back(-coef * CLA.values[l]);
                                newElColumns.push_back(CLA.columns[l]);
                            }                            
                        }
                        if (l + 1 == numNonZero)
                        {
                            break;
                        }
                    }

                    b[currentRow] -= coef * b[i];

                    CLA.values.erase(CLA.values.begin() + j );                       //не обнуляем нижние элементы, а удаляем из списка
                    CLA.rows.erase(CLA.rows.begin() + j );
                    CLA.columns.erase(CLA.columns.begin() + j );
                    --numNonZero;

                    vector<int>::iterator it;                                        //добавляем новые не 0 элементы, равые 0 - coef*a[i][l]
                    vector<double>::iterator itValue;

                    itValue = CLA.values.begin();
                    CLA.values.insert(itValue + j, newElements.begin(), newElements.end());

                    it = CLA.columns.begin();
                    CLA.columns.insert(it + j, newElColumns.begin(), newElColumns.end());

                    it = CLA.rows.begin();
                    CLA.rows.insert(it + j, newElements.size() , currentRow);
                    numNonZero += newElements.size();
                    newElColumns.clear();
                    newElements.clear();
                }
            }
        }

        int diagonal{ 0 };
        for (int i = N - 1; i >= 0; i--)                                             //обратный ход 
        {
            double sum{ 0 };
            for (int j = 0; j < numNonZero; j++)
            {
                if (CLA.rows[j] == i)
                {
                    sum += CLA.values[j] * x[CLA.columns[j]];
                    if(CLA.rows[j] == CLA.columns[j])
                    {
                        diagonal = j;
                    }
                }
            }
            x[i] = (b[i] - sum) / CLA.values[diagonal];                             //результирующий вектор 
        }
        return x;
     }

     //считает матрицу a и вектор b по методу Гаусса    
    vector<double> GaussWithCRSAAdapted(vector<double>& b)
    {
        vector<double> x(N, 0.0);                                   //собираем результат в вектор x

        for (int i = 0; i < N; i++)
        {
            int mainRow{ i };
            double mainElem{ 0.0 };
            int diagonal{ 0 };
            for (int j = CSRAA.nextRow[i]; j < numNonZero; j++)     //находим максимальный элемент в столбце i, алгоритм с главным элементом
            {
                if (CSRAA.columns[j] == i)
                    if (abs(CSRAA.values[j]) > abs(mainElem))
                    {
                        mainRow = j;
                        mainElem = CSRAA.values[j];
                        if (CLA.rows[j] == i)
                        {
                            diagonal = j;
                        }
                    }
            }
            if (mainElem == 0.0)                                    //нулевой столбец
            {
                continue;
            }

            if (mainRow != i)
            {
                for (int j = diagonal; j < numNonZero;j++)          //меняем местами текущую строку и строку с максимальным элементом
                {
                    if (CSRAA.rows[j] == i)
                    {
                        CSRAA.rows[j] = mainRow;
                        CSRAA.nextRow[i] = CSRAA.nextRow[i];
                        continue;
                    }
                    if (CSRAA.rows[j] == mainRow)
                    {
                        CSRAA.rows[j] = i;
                    }
                }
                double temp{ b[i] };                                //меняем местами вектор свободных членов                                
                b[i] = b[mainRow];
                b[mainRow] = temp;

                //************** TODO
            } 
        }
    }

     //считает матрицу в виде представления Кнута списков и вектор b по методу Гауса
     vector<double> GaussWithKnuthLists(vector<double>& b)
     {
         vector<double> x(N, 0.0);

         for (int i = 0; i < N; i++)
         {
             double mainElementValue{ 0.0 };
             Node* mainElem = nullptr;

             for (Node* element = KNUTHcolumns[i]; element->nextInColumn != nullptr; element = element->nextInColumn)
             {
                 if (abs(element->value) > mainElementValue && element->row >= i)
                 {
                     mainElementValue = abs(element->value);
                     mainElem = element;
                 }
             }
         

/*             if (mainElem->row != KNUTHrows[i]->row)                   //****** TODO
             {
                 swap(KNUTHrows[i], KNUTHrows[mainElem->row]);

                 for (Node* element = KNUTHrows[i]; element != nullptr; element = element->nextInRow)
                 {
                     swap(element->column, mainElem->column);
                 }
             }*/
/*
             for (Node* element = KNUTHrows[i]->nextInRow; element != nullptr; element = element->nextInRow)
             {
                 double coef = element->value / KNUTHrows[i]->value;

                 for (Node* pivot = KNUTHrows[i]->nextInColumn; pivot != nullptr; pivot = pivot->nextInColumn) 
                 {
                     bool found = false;

                     for (Node* colElement = pivot; colElement->nextInRow != nullptr; colElement = colElement->nextInRow) 
                     {
                         if (colElement->row == element->row) 
                         {
                             colElement->value -= coef * pivot->value;
                             found = true;
                             break;
                         }
                     }

                     if (!found) 
                     {
                         Node* newElement = new Node{ -coef * pivot->value, element->row, pivot->column, nullptr, nullptr };
                         Node* currentElement = KNUTHrows[element->row];

                         while (currentElement->nextInRow != nullptr && currentElement->nextInRow->column < newElement->column) 
                         {
                             currentElement = currentElement->nextInRow;
                         }

                         newElement->nextInRow = currentElement->nextInRow;
                         currentElement->nextInRow = newElement;

                         currentElement = KNUTHcolumns[pivot->column];

                         while (currentElement->nextInColumn != nullptr && currentElement->nextInColumn->row < newElement->row)

                         {
                             currentElement = currentElement->nextInColumn;
                         }

                         newElement->nextInColumn = currentElement->nextInColumn;
                         currentElement->nextInColumn = newElement;
                     }
                 }
             }*/ 
         }

         for (int i = N - 1; i >= 0; i--)                 //обратный ход, готов
         {
             double sum{ 0.0 };
             Node* currentRow = KNUTHrows[i];
             while (currentRow != nullptr) 
             {
                 int j = currentRow->column;
                 sum += currentRow->value * x[j];
                 currentRow = currentRow->nextInRow;
             }
             x[i] = (b[i] - sum) / KNUTHrows[i]->value;
         }
         return x;
     }

     // Функция для нахождения обратной матрицы
     vector<vector<double>> InverseMatrix(vector<vector<double>> a)
     {
         int N{ (int)a.size() };
         vector<vector<double>> extended( N, vector<double>(N * 2,0.0));
         vector<vector<double>>inverse(N, vector<double>(N,0));
         for (int i = 0; i < N; i++)
         {                                                  //добавляем единичную матрицу справа
             for (int j = 0; j < N; j++)
             {
                 extended[i][j] = a[i][j];
             }
             extended[i][i + N] = 1;
         }
         
         for (int i = 0; i < N; i++)                        //приводим расширенную матрицу к диагональному виду, проверим
         {      
             if (extended[i][i] == 0)
             {
                 int j = i + 1;
                 while (j < N && extended[j][i] == 0) 
                 {
                     j++;
                 }
                 if (j == N) 
                 {
                     cout << "Матрица вырожденная, обратной матрицы не существует." << std::endl;
                     return ERROR;
                 }
                 for (int k = 0; k < N * 2; k++)
                 {
                     swap(extended[i][k], extended[j][k]);
                 }
             }
             
             for (int j = 0; j < N; j++)                    //прямой ход Гаусса
             {
                 if (j != i)
                 {
                     double coef = extended[j][i] / extended[i][i];
                     for (int k = 0; k < N * 2; k++) 
                     {
                         extended[j][k] -= coef * extended[i][k];
                     }
                 }
             }

             for (int i = 0; i < N; i++)                    //нормализуем строки расширенной матрицы, диагональ левой матрицы = 1
             {
                 double coef = extended[i][i];
                 for (int j = 0; j < N * 2; j++) 
                 {
                     extended[i][j] /= coef;
                 }
             }

             for (int i = 0; i < N; i++)                    //извлекаем обратную матрицу из расширенной матрицы
             {
                 for (int j = 0; j < N; j++) 
                 {
                     inverse[i][j] = extended[i][j + N];
                 }
             }
         }
         return inverse;
     }

     //считает СНАУ методом Ньютона; ограничено число итераций: 100, но можно перегрузить в параметрах эпсилон - приближение и число итераций
     vector<double> Newton(const vector<double>& p0, vector<double>& x, double eps = 1e-6, int itMax = 100, int N = 8)
     {
         int iterations{ 0 };                            //ограничим число итераций: 100, но можно перегрузить в параметрах эпсилон - приближение и iter

         for (int it = 0; it < itMax; it++)
         {        
             x[8] =  x[4] + x[5];
             x[9] =  x[6] - x[4];
             x[10] = x[7] - x[5];
             x[11] = x[6] + x[7];

             double f1 = p0[0] - x[0] - x[4]  * x[4]  * k;  //вычисление значений функций и их производных
             double f2 = x[0]  - x[1] - x[5]  * x[5]  * k;
             double f3 = x[1]  - x[3] - x[7]  * x[7]  * k;
             double f4 = x[0]  - x[2] - x[11] * x[11] * k;
             double f5 = x[2]  - x[3] - x[9]  * x[9]  * k;
             double f6 = x[3] - p0[3] - x[8]  * x[8]  * k;
             double f7 = p0[1] - x[1] - x[6]  * x[6]  * k;
             double f8 = p0[2] - x[2] - x[10] * x[10] * k;
/*
             double f9  = p0[0];                            //добавим уравнения, чтобы не менять размер матрица Якоби
             double f10 = p0[1];
             double f11 = p0[2];
             double f12 = p0[3];
*/
             double df1dp1  = -1;                           //производные
             double df1dQ01 = -2 * x[4] * k;
             double df2dp2  = -1;
             double df2dQ12 = -2 * x[5] * k;
             double df3dp2  =  1;
             double df3dp4  = -1;
             double df3dQ24 = -2 * x[7] * k;
             double df4dp1  =  1;
             double df4dp3  = -1;
             double df4dQ13 = -2 * x[11] * k;
             double df5dp3  =  1;
             double df5dp4  = -1;
             double df5dQ34 = -2 * x[9] * k;
             double df6dp4  =  1;
             double df6dQ04 = -2 * x[8] * k;
             double df7dp2  = -1;
             double df7dQ02 = -2 * x[6] * k;
             double df8dp3  = -1;
             double df8dQ03 = -2 * x[10] * k;

             //                                           ***     Матрица Якоби    ***
             //signatureX' :: i :{ "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34", "Q01", "Q02", "Q03", "Q04" }
             //                       0    1     2     3     4      5      6       7      8     9      10    11             // 
             //            :: j  :{f1 ..f12}
             // 
             
             vector<vector<double>> jm= 
             { 
                                          {df1dp1, 0, 0, 0,       0, 0, 0, 0},
                                          {0, df2dp2, 0, 0,       df2dQ12, 0, 0, 0},
                                          {0, df3dp2, 0, df3dp4,  0, 0, df3dQ24, 0},
                                          {df4dp1, 0, df4dp3, 0,  0, df4dQ13, 0, 0},
                                          {0, 0, df5dp3, df5dp4,  0, 0, 0, df5dQ34},
                                          {0, 0, 0, df6dp4,       0, 0, 0, 0},
                                          {0, df7dp2, 0, 0,       0, 0, 0, 0},
                                          {0, 0, df8dp3, 0,       0, 0, 0, 0}

 /*                    {df1dp1,  0,       0,     df4dp1,  0,      0,      0,      0,  0,0,0,0},  //p1
                     { 0,    df2dp2, df3dp2,    0,       0,      0, df7dp2,      0,  0,0,0,0},
                     { 0,       0,       0,     df4dp3, df5dp3,  0      , 0, df8dp3, 0,0,0,0},
                     { 0,       0,   df3dp4,    0,      df5dp4, df6dp4,   0,     0 , 0,0,0,0},
                     { 0,    df2dQ12,    0,     0,       0,      0,       0,     0,  0,0,0,0},  //Q12
                     { 0,       0,       0,     df4dQ13, 0,      0,       0,     0 , 0 ,0 ,0 ,0  },  //Q13
                     { 0,       0,   df3dQ24,   0,       0,      0,       0,     0  ,0,0,0,0},  //Q24
                     { df1dQ01, 0,       0,     0,      df5dQ34, 0,       0,     0,0,0,0,0},  //Q34
                     { 0,       0,       0,     0,       0,      0,       0,     0,0,0,0,0  },  //Q01
                     { 0,        0 ,      0,     0,       0,      0, df7dQ02,     0 ,0,0,0,0 },
                     { 0,       0,       0,     0,       0,       0,    0,    df8dQ03,0,0,0,0},
                     { 0,       0 ,      0 ,    0 ,      0 ,   df6dp4,    0,     0,0,0,0,0} 
                     */
             };
             vector<vector<double>>  inverseJ = InverseMatrix(jm);

             vector<double> F = { -f1, -f2, -f3, -f4, -f5, -f6, -f7, -f8 };        //вектор функций

             vector<double> dp(N, 0);                                               //решение системы уравнений J * dp = F

             for (int i = 0; i < N; i++)
             {
                 for (int j = 0; j < N; j++)
                 {
                     dp[i] +=  inverseJ[i][j] * F[j];
                 }
             }

             for (int i = 0; i < N; i++)                                            //обновление значений переменных
             {
                 x[i] += dp[i]; cout << "\ndp " << i << " " << dp[i];
             }


             if (abs(dp[0]) < eps && abs(dp[1]) < eps && abs(dp[2]) < eps)          //проверка условия окончания итераций по эпсилон
             {
                break;
             }
         }         
         cout << p0[0] << "p0 " << x[0] << " p1 " << x[3] << " p4 " << " q12: " << x[4] << " q01 " << x[8];
              return x;
     }

    //
    //********************************* СЕРВИСНЫЕ ФУНКЦИИ ****************************
    // 
    //ввод начальных значений давления 
     vector<double> InputValues(int pInputCount)
     {
         vector<double> pInput(pInputCount, 0.0);
         int i{};
         do
         {
             cout << "Введите входное давление на трубопроводе номер " << i << " в Паскалях: ";
             cin >> pInput[i];
             if (pInput[i] >= 0)
                 i++;
         } while (i < pInputCount);
         return pInput;
     }

    //выводит на консоль сжатую матрицу в формате векторов не 0-ых элементов
    void ShowCompressedMatrix(const void *matrix, DataType dataType)
    {  
        auto m = (CompressedMatrix4V*) matrix;
        for (double value : m->values) 
        {
            cout << value << " ";
        }
        cout << endl << endl;

        cout << "Столбцы: ";
        for (int column : m->columns)
        {
            cout << column << " ";
        }
        cout << endl;

        cout << "Строки: ";
        for (int rows : m->rows)
        {
            cout << rows << " ";
        }
        cout << endl;

        if (dataType == DataType::CSRAA) 
        {
            cout << "Cжатые сторки: ";
            for (int rows : m->nextRow)
            {
                cout << rows << " ";
            }
            cout << endl;
        }
        cout << endl << line << endl;
    }

    //выводит на консоль сжатую матрицу в Кнута
    void ShowKnuthList(const vector<Node*> list, DataType dataType)
    {
        for (int i = 0; i < N; i++)
        {
            cout << endl << i;
            dataType == DataType::KNUTHRowsList ? cout << " cтрока : " : cout << "cтолбец : ";
            Node* element = list[i];
            while (element)
            {
                cout.width(12);
                cout << left << element->value;
                cout << left << "i=" << element->column;
                cout << left << " j=" << element->row <<",  ";
                dataType == DataType::KNUTHRowsList ? element = element->nextInRow : element = element->nextInColumn;
            }
        }
    }

    //выводит на консоль результирующий вектор x
    void ShowResultX(vector<double> x, bool newton = false)                                              //собираем результат в вектор x    
    {
        cout << endl << endl;
        for (int i = 0; i < N; i++)
        {
            auto title = (newton) ? signatureXN[i] : signatureX[i];
            cout << i + 1 << ".\t" << title << " = " << "\t" << x[i] << endl;
        }
    }

    //функция для вывода матрицы на экран
    void ShowSqMatrix(const vector<vector<double>> matrix)
    {
        int N = (int)matrix.size();
        cout << endl;
        for (int i = 0; i < N; i++) 
        {
            for (int j = 0; j < N; j++)
            {
                cout << "\t" << matrix[i][j] << " ";
            }
            cout << endl;
        }
    }

    //переводит char в double для считываемых из файла элементов матрицы
    double ConvertSymbolToValue(char symbol)
    {
        switch (symbol) 
        {
        case '0':
            return 0.0;
        case '1':
            return 1.0;
        case '-':
            return -1.0;
        case 'k':
            return 0.5;
        default:
            cout << line << "\n\tНе правильный формат данных!" << endl << line;
            return (int) symbol;
        }
    }

    //ввод матрицы из файла в текстовом виде с переводом в double, по флагу printFlag выводится на экран
    vector<vector<double>> ReadMatrix(const string& filename, bool printFlag = true)
    {
        vector<vector<double>> matrix;
        ifstream file(filename);

        if (file.is_open())
        {
            string line;
            while (getline(file, line))
            {
                vector<double> row;
                for (char symbol : line)
                {
                    double value = ConvertSymbolToValue(symbol);
                    row.push_back(value);
                    printFlag&& cout << value << " ";
                }
                matrix.push_back(row);
                printFlag&& cout << "\n";
            }
            file.close();
        }
        else
        {
            cout << line << "\n\tНевозможно открыть файл с матрицей : " << filename << endl << line;
        }
        return matrix;
    }

    //сохраняет результирующий вектор  в формате JSON в файл, ключ - название параметра
    void WriteJSON(const string& filename, const vector<string>& title, const vector<double>& x)
    {
        string jsonData = ConvertVectorToJson(title, x);
        ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            outputFile << jsonData;
            outputFile.close();
        }
        else
        {
            cout << line << "\n\tНевозможно открыть файл для записи : " << filename << endl << line;
        }
    }    

    //преобразования векторов в формат JSON
    string ConvertVectorToJson(const vector<string>& title, const vector<double>& x) 
    {
        string json = "{\n";

        for (size_t i = 0; i < title.size(); i++) 
        {
            json += "    \"" + title[i] + "\": " + to_string(x[i]);                  //собираем название параметра и его значение в строку json 

            if (i != title.size() - 1) 
            {
                json += ",";                                                         //строки JSON пишутся через запятую, кроме последней
            }
            json += "\n";
        }
        json += "}";
        return json;
    }

    //
    //коструктор
    //
    HydraulicNet(vector<vector<double>>& aMatrix, int numOfParameters, double resistK): a{ aMatrix }, N{numOfParameters}, k{resistK}
    {
        numNonZero = CountNonZeroElements(aMatrix);
    }
};

//прототип; вычисляет время выполнения ф-ции, в качестве параметров передатся объект и его метод, время выполнения которого измеряется
double GetCalculationRuntime(HydraulicNet, vector<double> (HydraulicNet::* calculation)(vector<double>&), vector<double>, vector<double>*);

//
//************************************* MAIN *******************************************
//
int main()
{
    setlocale(LC_ALL, "");

    Pipe pipe;
    double k = pipe.k;                          //для системы линейных уравнений сразу вычисляем коэфф. проводимости - гидравлическое сопротивление
 
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
        {0, 0, 0, 1, -1, 0, 0, 0, 0, k, 0, 0},
    };
    HydraulicNet net(a, N, k);
    vector<double> b;
    vector<double> pInput(PInputCount, 0);                                                      //вводим входные давления, 4 * Pвхю

    Mode mode = Mode::Debug;                                                                    

//************************************ ДЛЯ ТЕСТОВ !!!!! ******************************          //для тестов харкодим Po 1..4

    (mode == Mode::Release) ? pInput = net.InputValues(PInputCount) : pInput = {150, 130, 100, 50};
    
    
    b = { pInput[0], pInput[1], pInput[2], pInput[3], 0, 0, 0, 0, 0, 0, 0, 0 }; 
    vector<double> x(N, 0.0);                                                                   //собираем результат в вектор x

    cout << endl << "\tПреобразование исходной матрицы к виду Coordinate List:\n" << endl;
    net.CLA = net.ConvertToCLA(a);                                                 
    net.ShowCompressedMatrix(&net.CLA, DataType::CLA);

    cout << endl << "\tПреобразование исходной матрицы к виду CompressedSparseRow:\n" << endl;
    net.СSRA = net.ConvertToCSRA(a);      
    net.ShowCompressedMatrix(&net.СSRA, DataType::CSRA);

    cout << endl << "\tПреобразование исходной матрицы к виду CompressedSparseRowAdapted:\n" << endl;
    net.CSRAA = net.ConvertToCSRAAdapted(a);
    net.ShowCompressedMatrix(&net.CSRAA, DataType::CSRAA);

    cout << endl << "\tПреобразование исходной матрицы к виду Кнута:\n" << endl;
    net.ConvertToKnuthLists(a);
    net.ShowKnuthList (net.KNUTHrows, DataType::KNUTHRowsList);
    net.ShowKnuthList (net.KNUTHcolumns, DataType::KNUTHColumnsList);


    cout << "\n\n" << line << "\n\tРасчет давлений в узлах и объёмных расходов ы трубопроводах.\n" << line << endl;
   
    double funcRuntime = GetCalculationRuntime(net, &HydraulicNet::GaussWithKnuthLists, b, &x);
    net.ShowResultX(x);
    cout << "\nВремя расчетов  Гаусом c матрицей в варианте связных списков по модели Кнута, мкс: " << funcRuntime << endl;
    cout << line;

     funcRuntime = GetCalculationRuntime(net, &HydraulicNet::GaussWithCoordinateVectors, b, &x); 
    net.ShowResultX(x);
    cout << "\nВремя расчетов  Гаусом c матрицей в варианте Coordinate List, мкс: " << funcRuntime << endl;
    cout << line;

    funcRuntime = GetCalculationRuntime(net, &HydraulicNet::Gauss, b, &x);               // выводим время выполнения метода Гауса
    net.ShowResultX(x);
    cout << "\nВремя расчетов простым Гаусом, мкс: " << funcRuntime << endl;
    cout << line;

    net.WriteJSON(filenameResultJSON, net.signatureX, x);
    a = net.ReadMatrix(filenameM, (bool)mode);

    x = {100, 100 , 100 , 100 , 100, 100, 100, 100, 100, 100, 100, 100};
    x = net.Newton(pInput,x);
    net.ShowResultX(x, true);

    return 0;
} 

//
//вычисляет время выполнения ф-ции, в качестве параметров передатся объект и его метод, время выполнения метода измеряется
//
double GetCalculationRuntime(HydraulicNet net, vector<double> (HydraulicNet::* calculation)(vector<double>&), vector<double> matrixB, vector<double> *x)
{
    auto start = chrono::high_resolution_clock::now();                                          // запоминаем время начала выполнения функции

    (*x) = (net.*calculation)(matrixB);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);                    // вычисляем время выполнения функции в миллисекундах
    return duration.count();
}