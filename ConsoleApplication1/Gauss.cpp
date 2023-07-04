#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>

#include "gauss.h"

using namespace std;

const string line(100, '=');                    //������ ����������� ��� ������ �� �������
const string filenameM = "matrix.txt";          //������� ����� ��������� �� �����
const string filenameResultJSON = "result.json";//��������� ��������� ����������� � ������� JSON 
const vector<vector<double>> ERROR = { {1.7e308}, {0} };


vector<string> signatureX { "P1", "P2", "P3", "P4", "Q01", "Q12", "Q02", "Q24", "Q04", "Q34", "Q03", "Q13" }; //��������� ���������, �������
vector<string> signatureXN { "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34", "Q01", "Q02", "Q03", "Q04" }; //��������� ���������, �������
//double Q01{ 0 }, Q12{ 0 }, Q02{ 0 }, Q24{ 0 }, Q04{ 0 }, Q34{ 0 }, Q03{ 0 }, Q13{ 0 }, p1, p2, p3, p4;

const int PInputCount{ 4 };                     //����� ������� �������������, 4 �� ���������
const int N{ 12 };                              //����� ��������� ����������, 12 �� ���������, ��������� � ������������ � ������

//������ ��� �������
enum class Mode
{
    Debug = false,
    Release = true
};

//
// ���������� ������
// 
const double Pi{ 3.14159 };
const double ReCr{ 2400 };

//���������� ������ ��������
struct Liquid
{
    double nu{ 1.004e-6 };                      //�������������� �������� ���� ��� 293 �
    double rho{ 1000 };                         // ��������� �������� � �����
 };

//�������������� ������ ������������
struct Pipe
{
    Liquid water;                               //���������� ��������� ����

    double L{ 10 };                             //����� ����� �����
    double D{ 0.45 };                           //������� �����

    double ksi = 64 * L / (ReCr * D);
    double reLam = 16 * Pi * D / 0.03;
    double lam = 8 / (Pi * Pi * D * D * D * D);

public:
    double kTurb = ksi * lam / water.rho;      //������� ������������ - ����. ������������� ��� ������������� ������
    double k= ksi * lam * water.nu * reLam;    //������� ������������ - ����. ������������� ��� ����������� ������
};

//
// ��������� ������ ��� �������� ����������� �������
// 
enum class DataType
{
    CLA,
    CSRA,
    CSRAA,
    KNUTHRowsList,
    KNUTHColumnsList
};

//������� � ���� 3-�� ��������: ��������, ���������� (2 �������� �������������: ������ ������ i - CLA / ������ i ���� - CSRA)
struct CompressedMatrix
{
    vector<double> values;                      //��������� �������� �������
    vector<int> rows;                           //������� ����� ��������� ���������
    vector<int> columns;                        //������� �������� ��������� ���������
};

//������� � ���� 4 ��������: ��������, ���������, ����� ��������� � ������ (CSR matrix + vector j)
struct CompressedMatrix4V : CompressedMatrix
{
   vector<int> nextRow;                         //��������� �� ����� �������� - ������ ������
};

//������� � ���� ������� ������� �� ����� �����
struct Node {
    double value;
    int column;
    int row;
    Node* nextInRow;
    Node* nextInColumn;
};

//
// �������������� ������
//
//����� ��� �������� �������������� ����. 
class HydraulicNet
{
    const int N{ 12 };                          //����� ��������� ����������, 12 �� ���������, ��������� � ������������
    const double k{ 0.5 };                      //���������� ������������, �������������� �������������, ��������� � �����������

    vector<vector<double>> a;                   //����� ����� �������, �������������� ����
 
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
    CompressedMatrix CLA;                       //Coordinate List matrix A, 3 ������� �� 0 ���������: ���������, ���������� ����, �������
    CompressedMatrix �SRA;                      //Compressed Sparse Row matrix A, , ������ ��������� ����� ����
    CompressedMatrix4V CSRAA;                   //������ 4-� ��������: ��������, ���������� i, ���������� j, ����� ��������� � ������, ���������� CLA � CRSA
    vector<Node*> KNUTHrows;                    //������� ������ �� ������ �����, ����������� ������� ����� � ����� ��� � ����� �������
    vector<Node*> KNUTHcolumns;

    vector<string> signatureX { "P1", "P2", "P3", "P4", "Q01", "Q12", "Q02", "Q24", "Q04", "Q34", "Q03", "Q13" }; //��������� ���������, ������� ��� ������
    vector<string> signatureXN { "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34", "Q01", "Q02", "Q03", "Q04" }; //��������� ���������, ������� ��� �������
 //   double Q01{ 0 }, Q12{ 0 }, Q02{ 0 }, Q24{ 0 }, Q04{ 0 }, Q34{ 0 }, Q03{ 0 }, Q13{ 0 };

    //
    //********************************* ����������� ����������� ������� *************
    //
    //������������ ����������� ������� � ������ 3-� ��������: ��������, ���������� i, ���������� j
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

    //������������ ����������� ������� � ������ Sparse Row
    CompressedMatrix ConvertToCSRA(const vector<vector<double>>& matrix)
    {
        vector<double>values;
        vector<int> columns;
        vector<int> rowsCumulative;        
        CompressedMatrix CSRA;

        int c = 0;                              //������������ ����� �� ������� ��������� � ������
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

    //������������ ����������� ������� � ������ 4-� ��������: ��������, ���������� i, ���������� j, ����� ��������� � ������
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

    //������������ ����������� ������� � ������ ����� - ������� ������ ��������, ������ �����, ������ ������� ������
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

                    if (KNUTHrows[i] == nullptr)                             //���� � ����� ������
                    {
                        KNUTHrows[i] = element;
                    }
                    else
                    {
                        Node* currentRow = KNUTHrows[i];
                        while (currentRow->nextInRow != nullptr)             //������������ � ����� ������, ����� ������� ����� �� ���������� �������
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
    //********************************* �������������� ���� **************************
    // 
    //������� ������� a � ������ b �� ������ ������    
    vector<double> Gauss(vector<double>& b)
    {
        vector<double> x(N, 0.0);               //�������� ��������� � ������ x

        for (int i = 0; i < N; i++)
        {
            int maxRow{ i }; 
            for (int j = i + 1; j < N; j++)     //������� ������������ ������� � ������� i, �������� � ������� ���������
            {
                if (abs(a[j][i]) > abs(a[maxRow][i]))
                {
                    maxRow = j;
                }
            }

            if (a[maxRow][i] == 0.0)            //������� �������
            {
                continue;
            }

            if (maxRow != i)
            {
                swap(a[i], a[maxRow]);          //������ ������� ������� ������ � ������ � ������������ ���������
                swap(b[i], b[maxRow]);
            }

            for (int j = i + 1; j < N; j++)     //������ ��� 
            {
                if (a[j][i] != 0)
                {
                    double coef = a[j][i] / a[i][i];
                    for (int k = i + 1; k < N; k++)
                    {
                        a[j][k] -= coef * a[i][k];
                    }
                    b[j] -= coef * b[i];
                    a[j][i] = 0;               //�������� ��������� �������, ���: -1 �������� � ������ 0,  ������ Aji - Aji*Aii/Aii
                }
            }
        }
 
        for (int i = N - 1; i >= 0; i--)       //�������� ��� 
        {     
            double sum{ 0.0 };
            for (int j = i + 1; j < N; j++)
            {
                sum += a[i][j] * x[j];          //����������� ��������� � � ������ ����
            }
            x[i] = (b[i] - sum) / a[i][i];      //�������������� ������ 
        }
        return x;
    }

    //������� ������� � ���� Coordinate List � ������ b �� ������ ����� 
    vector<double> GaussWithCoordinateVectors(vector<double>& b)
     {
         vector<double> x(N, 0.0);                              //�������� ��������� � ������ x

        for (int i = 0; i < N; i++)
        {
            double mainElemAbs{ 0.0 }, mainElem{ 0.0 };
            int mainRow{ 0 };                                   //��� �������� ��������
            int mainIndex{ 0 };                                 //������ �������� ��������
            for (int j = 0; j < numNonZero; j++)                //������� ������������ ������� � ������� i, �������� � ������� ���������
            {
                if (CLA.columns[j] == i && abs(CLA.values[j]) > mainElemAbs && CLA.rows[j] >= i)
                {
                    mainElemAbs = abs(CLA.values[j]);
                    mainElem = CLA.values[j];
                    mainRow = CLA.rows[j];
                    mainIndex = j;
                }
            }

            if (mainElem == 0.0)                                //������� �������
            {
                continue;
            }
 
            if (mainRow != i)
            {
                for (int j = i; j < numNonZero;j++)              //������ ������� ������� ������ � ������ � ������������ ���������
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
                double temp{ b[i] };                            //������ ������� ������ ��������� ������                                
                b[i] = b[mainRow];
                b[mainRow] = temp;
            }  

            bool zeroAtCurrentRow;                              //�������, ��� � ������ ���� 0 �������, ������� ����� ��������� ����� !=0
            vector<double> newElements;                         //��������, �������  ��������� ���� 0 ����� ��������� �����
            vector<int> newElColumns;                           //����� ������� ����� ���������

            for (int j = 0; j < numNonZero; j++)                //������ ��� 
            {
                if ((CLA.columns[j] == i) && (CLA.rows[j] > i))
                {
                    double coef = CLA.values[j] / mainElem;
                    int currentRow = CLA.rows[j];
                    for (int l = mainIndex; CLA.rows[l] == mainRow; l++)             //������ mainRow � ������� ������� ���������
                    {
                        zeroAtCurrentRow = true;
                        if (CLA.rows[l] == i  && CLA.columns[l] > i)
                        {
                            for (int k = j; CLA.rows[k] == currentRow; k++)          //������ workRow ��� ��������� � ���������
                            {
                                if (CLA.columns[k] == CLA.columns[l])                //�� 0 �������� � ����� ������� mainRow � workRow
                                {
                                    CLA.values[k] -= coef * CLA.values[l];
                                    zeroAtCurrentRow = false;
                                }
                                if (k + 1 == numNonZero)
                                {
                                    break;
                                }
                            }
                            if (zeroAtCurrentRow)                                    //�� 0 ������� ������ � mainRow
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

                    CLA.values.erase(CLA.values.begin() + j );                       //�� �������� ������ ��������, � ������� �� ������
                    CLA.rows.erase(CLA.rows.begin() + j );
                    CLA.columns.erase(CLA.columns.begin() + j );
                    --numNonZero;

                    vector<int>::iterator it;                                        //��������� ����� �� 0 ��������, ����� 0 - coef*a[i][l]
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
        for (int i = N - 1; i >= 0; i--)                                             //�������� ��� 
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
            x[i] = (b[i] - sum) / CLA.values[diagonal];                             //�������������� ������ 
        }
        return x;
     }

     //������� ������� a � ������ b �� ������ ������    
    vector<double> GaussWithCRSAAdapted(vector<double>& b)
    {
        vector<double> x(N, 0.0);                                   //�������� ��������� � ������ x

        for (int i = 0; i < N; i++)
        {
            int mainRow{ i };
            double mainElem{ 0.0 };
            int diagonal{ 0 };
            for (int j = CSRAA.nextRow[i]; j < numNonZero; j++)     //������� ������������ ������� � ������� i, �������� � ������� ���������
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
            if (mainElem == 0.0)                                    //������� �������
            {
                continue;
            }

            if (mainRow != i)
            {
                for (int j = diagonal; j < numNonZero;j++)          //������ ������� ������� ������ � ������ � ������������ ���������
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
                double temp{ b[i] };                                //������ ������� ������ ��������� ������                                
                b[i] = b[mainRow];
                b[mainRow] = temp;

                //************** TODO
            } 
        }
    }

     //������� ������� � ���� ������������� ����� ������� � ������ b �� ������ �����
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

         for (int i = N - 1; i >= 0; i--)                 //�������� ���, �����
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

     // ������� ��� ���������� �������� �������
     vector<vector<double>> InverseMatrix(vector<vector<double>> a)
     {
         int N{ (int)a.size() };
         vector<vector<double>> extended( N, vector<double>(N * 2,0.0));
         vector<vector<double>>inverse(N, vector<double>(N,0));
         for (int i = 0; i < N; i++)
         {                                                  //��������� ��������� ������� ������
             for (int j = 0; j < N; j++)
             {
                 extended[i][j] = a[i][j];
             }
             extended[i][i + N] = 1;
         }
         
         for (int i = 0; i < N; i++)                        //�������� ����������� ������� � ������������� ����, ��������
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
                     cout << "������� �����������, �������� ������� �� ����������." << std::endl;
                     return ERROR;
                 }
                 for (int k = 0; k < N * 2; k++)
                 {
                     swap(extended[i][k], extended[j][k]);
                 }
             }
             
             for (int j = 0; j < N; j++)                    //������ ��� ������
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

             for (int i = 0; i < N; i++)                    //����������� ������ ����������� �������, ��������� ����� ������� = 1
             {
                 double coef = extended[i][i];
                 for (int j = 0; j < N * 2; j++) 
                 {
                     extended[i][j] /= coef;
                 }
             }

             for (int i = 0; i < N; i++)                    //��������� �������� ������� �� ����������� �������
             {
                 for (int j = 0; j < N; j++) 
                 {
                     inverse[i][j] = extended[i][j + N];
                 }
             }
         }
         return inverse;
     }

     //������� ���� ������� �������; ���������� ����� ��������: 100, �� ����� ����������� � ���������� ������� - ����������� � ����� ��������
     vector<double> Newton(const vector<double>& p0, vector<double>& x, double eps = 1e-6, int itMax = 100, int N = 8)
     {
         int iterations{ 0 };                            //��������� ����� ��������: 100, �� ����� ����������� � ���������� ������� - ����������� � iter

         for (int it = 0; it < itMax; it++)
         {        
             x[8] =  x[4] + x[5];
             x[9] =  x[6] - x[4];
             x[10] = x[7] - x[5];
             x[11] = x[6] + x[7];

             double f1 = p0[0] - x[0] - x[4]  * x[4]  * k;  //���������� �������� ������� � �� �����������
             double f2 = x[0]  - x[1] - x[5]  * x[5]  * k;
             double f3 = x[1]  - x[3] - x[7]  * x[7]  * k;
             double f4 = x[0]  - x[2] - x[11] * x[11] * k;
             double f5 = x[2]  - x[3] - x[9]  * x[9]  * k;
             double f6 = x[3] - p0[3] - x[8]  * x[8]  * k;
             double f7 = p0[1] - x[1] - x[6]  * x[6]  * k;
             double f8 = p0[2] - x[2] - x[10] * x[10] * k;
/*
             double f9  = p0[0];                            //������� ���������, ����� �� ������ ������ ������� �����
             double f10 = p0[1];
             double f11 = p0[2];
             double f12 = p0[3];
*/
             double df1dp1  = -1;                           //�����������
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

             //                                           ***     ������� �����    ***
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

             vector<double> F = { -f1, -f2, -f3, -f4, -f5, -f6, -f7, -f8 };        //������ �������

             vector<double> dp(N, 0);                                               //������� ������� ��������� J * dp = F

             for (int i = 0; i < N; i++)
             {
                 for (int j = 0; j < N; j++)
                 {
                     dp[i] +=  inverseJ[i][j] * F[j];
                 }
             }

             for (int i = 0; i < N; i++)                                            //���������� �������� ����������
             {
                 x[i] += dp[i]; cout << "\ndp " << i << " " << dp[i];
             }


             if (abs(dp[0]) < eps && abs(dp[1]) < eps && abs(dp[2]) < eps)          //�������� ������� ��������� �������� �� �������
             {
                break;
             }
         }         
         cout << p0[0] << "p0 " << x[0] << " p1 " << x[3] << " p4 " << " q12: " << x[4] << " q01 " << x[8];
              return x;
     }

    //
    //********************************* ��������� ������� ****************************
    // 
    //���� ��������� �������� �������� 
     vector<double> InputValues(int pInputCount)
     {
         vector<double> pInput(pInputCount, 0.0);
         int i{};
         do
         {
             cout << "������� ������� �������� �� ������������ ����� " << i << " � ��������: ";
             cin >> pInput[i];
             if (pInput[i] >= 0)
                 i++;
         } while (i < pInputCount);
         return pInput;
     }

    //������� �� ������� ������ ������� � ������� �������� �� 0-�� ���������
    void ShowCompressedMatrix(const void *matrix, DataType dataType)
    {  
        auto m = (CompressedMatrix4V*) matrix;
        for (double value : m->values) 
        {
            cout << value << " ";
        }
        cout << endl << endl;

        cout << "�������: ";
        for (int column : m->columns)
        {
            cout << column << " ";
        }
        cout << endl;

        cout << "������: ";
        for (int rows : m->rows)
        {
            cout << rows << " ";
        }
        cout << endl;

        if (dataType == DataType::CSRAA) 
        {
            cout << "C����� ������: ";
            for (int rows : m->nextRow)
            {
                cout << rows << " ";
            }
            cout << endl;
        }
        cout << endl << line << endl;
    }

    //������� �� ������� ������ ������� � �����
    void ShowKnuthList(const vector<Node*> list, DataType dataType)
    {
        for (int i = 0; i < N; i++)
        {
            cout << endl << i;
            dataType == DataType::KNUTHRowsList ? cout << " c����� : " : cout << "c������ : ";
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

    //������� �� ������� �������������� ������ x
    void ShowResultX(vector<double> x, bool newton = false)                                              //�������� ��������� � ������ x    
    {
        cout << endl << endl;
        for (int i = 0; i < N; i++)
        {
            auto title = (newton) ? signatureXN[i] : signatureX[i];
            cout << i + 1 << ".\t" << title << " = " << "\t" << x[i] << endl;
        }
    }

    //������� ��� ������ ������� �� �����
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

    //��������� char � double ��� ����������� �� ����� ��������� �������
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
            cout << line << "\n\t�� ���������� ������ ������!" << endl << line;
            return (int) symbol;
        }
    }

    //���� ������� �� ����� � ��������� ���� � ��������� � double, �� ����� printFlag ��������� �� �����
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
            cout << line << "\n\t���������� ������� ���� � �������� : " << filename << endl << line;
        }
        return matrix;
    }

    //��������� �������������� ������  � ������� JSON � ����, ���� - �������� ���������
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
            cout << line << "\n\t���������� ������� ���� ��� ������ : " << filename << endl << line;
        }
    }    

    //�������������� �������� � ������ JSON
    string ConvertVectorToJson(const vector<string>& title, const vector<double>& x) 
    {
        string json = "{\n";

        for (size_t i = 0; i < title.size(); i++) 
        {
            json += "    \"" + title[i] + "\": " + to_string(x[i]);                  //�������� �������� ��������� � ��� �������� � ������ json 

            if (i != title.size() - 1) 
            {
                json += ",";                                                         //������ JSON ������� ����� �������, ����� ���������
            }
            json += "\n";
        }
        json += "}";
        return json;
    }

    //
    //����������
    //
    HydraulicNet(vector<vector<double>>& aMatrix, int numOfParameters, double resistK): a{ aMatrix }, N{numOfParameters}, k{resistK}
    {
        numNonZero = CountNonZeroElements(aMatrix);
    }
};

//��������; ��������� ����� ���������� �-���, � �������� ���������� ��������� ������ � ��� �����, ����� ���������� �������� ����������
double GetCalculationRuntime(HydraulicNet, vector<double> (HydraulicNet::* calculation)(vector<double>&), vector<double>, vector<double>*);

//
//************************************* MAIN *******************************************
//
int main()
{
    setlocale(LC_ALL, "");

    Pipe pipe;
    double k = pipe.k;                          //��� ������� �������� ��������� ����� ��������� �����. ������������ - �������������� �������������
 
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
    vector<double> pInput(PInputCount, 0);                                                      //������ ������� ��������, 4 * P���

    Mode mode = Mode::Debug;                                                                    

//************************************ ��� ������ !!!!! ******************************          //��� ������ �������� Po 1..4

    (mode == Mode::Release) ? pInput = net.InputValues(PInputCount) : pInput = {150, 130, 100, 50};
    
    
    b = { pInput[0], pInput[1], pInput[2], pInput[3], 0, 0, 0, 0, 0, 0, 0, 0 }; 
    vector<double> x(N, 0.0);                                                                   //�������� ��������� � ������ x

    cout << endl << "\t�������������� �������� ������� � ���� Coordinate List:\n" << endl;
    net.CLA = net.ConvertToCLA(a);                                                 
    net.ShowCompressedMatrix(&net.CLA, DataType::CLA);

    cout << endl << "\t�������������� �������� ������� � ���� CompressedSparseRow:\n" << endl;
    net.�SRA = net.ConvertToCSRA(a);      
    net.ShowCompressedMatrix(&net.�SRA, DataType::CSRA);

    cout << endl << "\t�������������� �������� ������� � ���� CompressedSparseRowAdapted:\n" << endl;
    net.CSRAA = net.ConvertToCSRAAdapted(a);
    net.ShowCompressedMatrix(&net.CSRAA, DataType::CSRAA);

    cout << endl << "\t�������������� �������� ������� � ���� �����:\n" << endl;
    net.ConvertToKnuthLists(a);
    net.ShowKnuthList (net.KNUTHrows, DataType::KNUTHRowsList);
    net.ShowKnuthList (net.KNUTHcolumns, DataType::KNUTHColumnsList);


    cout << "\n\n" << line << "\n\t������ �������� � ����� � �������� �������� � �������������.\n" << line << endl;
   
    double funcRuntime = GetCalculationRuntime(net, &HydraulicNet::GaussWithKnuthLists, b, &x);
    net.ShowResultX(x);
    cout << "\n����� ��������  ������ c �������� � �������� ������� ������� �� ������ �����, ���: " << funcRuntime << endl;
    cout << line;

     funcRuntime = GetCalculationRuntime(net, &HydraulicNet::GaussWithCoordinateVectors, b, &x); 
    net.ShowResultX(x);
    cout << "\n����� ��������  ������ c �������� � �������� Coordinate List, ���: " << funcRuntime << endl;
    cout << line;

    funcRuntime = GetCalculationRuntime(net, &HydraulicNet::Gauss, b, &x);               // ������� ����� ���������� ������ �����
    net.ShowResultX(x);
    cout << "\n����� �������� ������� ������, ���: " << funcRuntime << endl;
    cout << line;

    net.WriteJSON(filenameResultJSON, net.signatureX, x);
    a = net.ReadMatrix(filenameM, (bool)mode);

    x = {100, 100 , 100 , 100 , 100, 100, 100, 100, 100, 100, 100, 100};
    x = net.Newton(pInput,x);
    net.ShowResultX(x, true);

    return 0;
} 

//
//��������� ����� ���������� �-���, � �������� ���������� ��������� ������ � ��� �����, ����� ���������� ������ ����������
//
double GetCalculationRuntime(HydraulicNet net, vector<double> (HydraulicNet::* calculation)(vector<double>&), vector<double> matrixB, vector<double> *x)
{
    auto start = chrono::high_resolution_clock::now();                                          // ���������� ����� ������ ���������� �������

    (*x) = (net.*calculation)(matrixB);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);                    // ��������� ����� ���������� ������� � �������������
    return duration.count();
}