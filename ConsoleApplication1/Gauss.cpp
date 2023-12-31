#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdarg.h>

#include "gauss.h"

using namespace std;

string line(100, '=');                    //������ ����������� ��� ������ �� �������
const string filenameM = "matrix.txt";          //������� ����� ��������� �� �����
const string filenameResultJSON = "result.json";//��������� ��������� ����������� � ������� JSON 

const double error = 1.7e308;					//����, ��� ��������� ������ ERROR_TYPE � ������� ERRROR ��� ERROR 2
vector<double> ERROR = { 0 , 0 };
vector<vector<double>> ERROR2 = { {ERROR}, {0, 0} };

enum class ERROR_TYPE							
{
	GAUSS = 1,
	GAUSSCLA,
	GAUSSCSRAA,
	GAUSSKNUTH,
	NEWTON,
	NEWTONHARDCODE,
	INVERSE,
	FILE,
	Max
};

//������ ��� �������
enum class Mode
{
	DEBUG = false,
	RELEASE = true
};

const string LOGfile = "param.log";

//����
void LOG(string comment, const int num, ...)
{
	ofstream outputFile;
	outputFile.open(LOGfile, ios_base::app);

	if (outputFile.is_open())
	{
		string result{ comment };
		result += " : ";
		va_list factor;
		va_start(factor, num);
		for (int i = 0;i < num; i++)
		{
			result += to_string(va_arg(factor, double));						// �������� �������� �������� ��������� ���� int
			result += " # ";
		}
		va_end(factor);

		outputFile << result << endl;
		outputFile.close();
	}
	else
	{
		cout << "\n\t���������� ������� ���� ��� ������ Log'��: " << LOGfile << endl << line;
		ERROR[0] = error;
		ERROR[1] = (double)ERROR_TYPE::FILE;
	}
}

//enum �����
enum class Task
{
	Gauss = 0,
	GaussCLA,
	GaussCSRAA,
	GaussKnuth,
	Newton1,
	NewtonHardCode4,
	All,
	MaxTask = 7
};

const int PInputCount{ 4 };                     //����� ������� �������������, 4 �� ���������
const int N{ 12 };                              //����� ��������� ����������, 12 �� ���������, ��������� � ������������ � ������

//
// ��������� ������ ��� �������� ����������� �������
// 
//enum ������������ �������� ������
enum class DataType
{
	CLA,
	CSRA,
	CSRAA,
	MATRIX,
	KNUTHRowsList,
	KNUTHColumnsList,
	ALL
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

//������� � ���� ������� ������� �� ����� �����, Node - ������� ������
struct Node {
	double value;
	int column;
	int row;
	Node* nextInRow;
	Node* nextInColumn;
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
	double k = ksi * lam * water.nu * reLam;    //������� ������������ - ����. ������������� ��� ����������� ������
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

	vector<string> signatureXN { "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34", "Q01", "Q02", "Q03", "Q04" }; //��������� ���������, ������� ��� �������

	//
	//********************************* ����������� ����������� ������� *************
	//
	//
	//������������ ����������� ������� � ������ 3-� ��������: ��������, ���������� i, ���������� j
	CompressedMatrix ConvertToCLA(const vector<vector<double>>& matrix)
	{
 		vector<double> values;
		vector<int> columns;
		vector<int> rows;
		CompressedMatrix CLAVectors;

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

		CLAVectors.columns = columns;
		CLAVectors.rows = rows;
		CLAVectors.values = values;
		return CLAVectors;
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
	//������� ������� � � ������ � �� ������ ������    
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

	//������� �������, �������������� � ���� Coordinate List � ������ � �� ������ ����� 
	vector<double> GaussCLA(vector<double>& b)
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
				cout << "\n!\n";
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
				cout << endl << " change i to moain: " << i << " " << mainRow;
				mainRow = i;
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
						if (CLA.rows[l] == i && CLA.columns[l] >= i)
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

					CLA.values.erase(CLA.values.begin() + j);                        //�� �������� ������ ��������, � ������� �� ������
					CLA.rows.erase(CLA.rows.begin() + j);
					CLA.columns.erase(CLA.columns.begin() + j);
					--numNonZero;

					vector<int>::iterator it;                                        //��������� ����� �� 0 ��������, ����� 0 - coef*a[i][l]
					vector<double>::iterator itValue;

					itValue = CLA.values.begin();
					CLA.values.insert(itValue + j, newElements.begin(), newElements.end());

					it = CLA.columns.begin();
					CLA.columns.insert(it + j, newElColumns.begin(), newElColumns.end());

					it = CLA.rows.begin();
					CLA.rows.insert(it + j, newElements.size(), currentRow);
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
					if (CLA.rows[j] == CLA.columns[j])
					{
						diagonal = j;
					}
				}
			}
			x[i] = (b[i] - sum) / CLA.values[diagonal];                             //�������������� ������ 
		}
		return x;
	}

	//������� ������� � � ������ � �� ������ ������, ������� a � ���� ������ �� 0 ��������� � ���. ��������: ������������ ���-�� �� 0 ��������� � ������   
	vector<double> GaussCRSAA(vector<double>& b)
	{
		ERROR[0] = error;
		ERROR[1] = (double)ERROR_TYPE::GAUSSCSRAA;
		return ERROR;
	}

	//������� �������, �������������� � �������� ����� - ������ ������ ����� � ��������
	vector<double> GaussKnuth(vector<double>& b)
	{/*
		vector<double> x(N, 0.0);

		for (int i = 0; i < N; i++)
		{
			for (Node* pivot = KNUTHrows[i]->nextInColumn; pivot != nullptr; pivot = pivot->nextInColumn)
			{
				double coef = pivot->value / KNUTHrows[i]->value;
				for (Node* element = KNUTHrows[i]; element != nullptr; element = element->nextInRow)
				{
					bool zeroCurrentRow = true;
					for (Node* colElement = pivot; colElement != nullptr; colElement = colElement->nextInRow)
					{
						if (colElement->row == element->row)
						{
							colElement->value -= coef * element->value;
							zeroCurrentRow = false;
							break;
						}
					}
					
					if (zeroCurrentRow)
					{
						Node* newElement = new Node{ -coef * element->value, element->column, pivot->row, nullptr, nullptr };
						Node* currentElement = KNUTHrows[pivot->row];

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
				b[pivot->row] -= b[i] * coef;
				pivot->value = 0;
			}
			
			for (int j = i + 1; j < N; i++)
			{
				KNUTHrows[j] = KNUTHrows[j]->nextInRow;
			}
			KNUTHrows[i]->nextInColumn = nullptr;
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
		for (int i = 0; i < N; i++) 
		{
			for (int j = 0; j <N; j++)
			for (Node* e = KNUTHrows[i]; e != nullptr; e = e->nextInRow)
			{	
				(j == e->column) ? cout << "\t" << e->value : cout << 0;
			}
			cout << endl;
		}
		*/
		ERROR[0] = error;
		ERROR[1] = (double)ERROR_TYPE::GAUSSKNUTH;
		return ERROR;
	}

	//��� ������������ ������� ����� ������������ ��������� ����� deltaY / deltaX
	double GetDiff(double (*equation)(double), double x, double dp = 1.0e-6)
	{
		double derivative = (equation(x + dp) - equation(x - dp)) / (2 * dp);
		return derivative;
	}
	

	//������� �������� ���� ���������
	double Newton(double (*equation)(double), double x, double epsilon = 1e-06, int it = 100)
	{
		double dx{}, df{};
		for (int i = 0; i < it; i++) 
		{
			double f = equation(x);
			df = GetDiff(equation, x, epsilon);
			if (df != 0) {
				dx = f / df;
			}
			else 
			{
				ERROR[0] = error;
				ERROR[1] = (double)ERROR_TYPE::NEWTON;
				return 0;
			}
			x -= dx;
//			LOG("��� �������� �������1 f # df # dx # x ", 4, f, df, dx, x);
			if (abs(f) < epsilon)
			{
				return x;
			}
		}
		return x;
	}

	//������� ��� ���������� �������� �������
	vector<vector<double>> InverseMatrix(vector<vector<double>> a)
	{
		int N{ (int)a.size() };
		vector<vector<double>> extended(N, vector<double>(N * 2, 0.0));
		vector<vector<double>>inverse(N, vector<double>(N, 0));
		for (int i = 0; i < N; i++)
		{													  //��������� ��������� ������� ������
			for (int j = 0; j < N; j++)
			{
				extended[i][j] = a[i][j];
			}
			extended[i][i + N] = 1;
		}

		for (int i = 0; i < N; i++)						      //�������� ����������� ������� � ������������� ����, �������� ����������
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
					ERROR2[0][0] = error;
					ERROR2[0][1] = (double)ERROR_TYPE::INVERSE;
					return ERROR2;
				}
				swap(extended[i], extended[j]);
			}

			for (int j = i + 1; j < N; j++)                    //������ ��� ������
			{
				if (extended[j][i] != 0)
				{
					double coef = extended[j][i] / extended[i][i];
					for (int k = i + 1; k < N * 2; k++)
					{
						extended[j][k] -= coef * extended[i][k];
					}
					extended[j][i] = 0.0;
				}
			}
		}

		for (int i = 0; i < N; i++)								 //����������� ������ ����������� �������, ��������� ����� ������� = 1
		{
			double coef = extended[i][i];
			for (int j = i; j < N * 2; j++)
			{
				extended[i][j] /= coef;
			}
		}

		for (int i = 0; i < N; i++)								 //��������� �������� ������� �� ����������� �������
		{
			for (int j = 0; j < N; j++)
			{
				inverse[i][j] = extended[i][j + N];
			}
		}

		return inverse;
	}

	//TODO            typedef double (*FunctionArr[])(double);
	//������� ���� ������� �������; ���������� ����� ��������: 100, �� ����� ����������� � ���������� ������� - ����������� � ����� ��������
	vector<double> NewtonHardCode(vector<double>& x)
	{
		const double eps = 1e-6;			
		int N =8, itMax = 100;
		double stepK = 1;

		vector<double> p0 = { x[8], x[9], x[10], x[11] };
		int iterations{ 0 };						             //��������� ����� ��������: 100, �� ����� ����������� � ���������� ������� - ����������� � iter
		int sumF{ 0 };

		for (int it = 0; it < itMax; it++)
		{			
										     //��������� � ���� ������� �������� � Gauss.h
		//
		//      p01---->p1---->p2<----p02
		//              |      |
		//              v      v
		//      p03---->p3---->p4----->p04
		// 
			 //���������� �������� ������� � �� �����������
			double f1 = p0[0] - x[0]  - x[4] * x[4] * k - 2 * k * x[4] * x[5] - k * x[5] * x[5];			
			double f2 = x[0]  - x[1]  - x[4] * x[4] * k;
			double f3 = x[1]  - x[3]  - x[6] * x[6] * k;
			double f4 = x[0]  - x[2]  - x[5] * x[5] * k;
			double f5 = x[2]  - x[3]  - x[7] * x[7] * k;
			double f6 = x[3]  - p0[3] - x[6] * x[6] * k - 2 * k * x[6] * x[7] - k * x[7] * x[7];
			double f7 = p0[1] - x[1]  - x[6] * x[6] * k + 2 * k * x[6] * x[4] - k * x[4] * x[4];
			double f8 = p0[2] - x[2]  - x[7] * x[7] * k + 2 * k * x[7] * x[5] - k * x[5] * x[5];
			/*
						 double f9  = p0[0];                     //������� ���������, ����� �� ������ ������ ������� �����
						 double f10 = p0[1];
						 double f11 = p0[2];
						 double f12 = p0[3];
			*/
			double df1dp1 = -1;									 //�����������
			double df1dQ12 = -2 * x[4] * k - 2 * k * x[5];
			double df1dQ13 = -2 * x[5] * k - 2 * k * x[4];

			double df2dp1  =  1;
			double df2dp2  = -1;
			double df2dQ12 = -2 * x[4] * k;
			
			double df3dp2  = 1;
			double df3dp4 = -1;
			double df3dQ24= -2 * x[6] * k;
			
			double df4dp1  =  1;
			double df4dp3  = -1;
			double df4dQ13 = -2 * x[5] * k;
			
			double df5dp3  = 1;
			double df5dp4 = -1;
			double df5dQ34= -2 * x[7] * k;
			
			double df6dp4   = 1;
			double df6dQ24 = -2 * x[6] * k - 2 * k * x[7];
			double df6dQ34 = -2 * x[7] * k - 2 * k * x[6];

			double df7dp2  = -1;
			double df7dQ24 = -2 * x[6] * k + 2 * k * x[4];
			double df7dQ12 =  2 * x[6] * k - 2 * k * x[4];

			double df8dp3  = -1;
			double df8dQ34 = -2 * x[7] * k + 2 * k * x[5];
			double df8dQ13 =  2 * x[7] * k - 2 * k * x[5];


			//                                           ***     ������� �����    ***
			//signatureX' :: i :{ "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34",       "Q01", "Q02", "Q03", "Q04" }
			//                       0    1     2     3     4      5      6       7            8     9      10    11             // 
			//            :: j  :{f1 ..f12}
			// 

			vector<vector<double>> jm =
			{
										 {df1dp1,       0,        0,        0,     df1dQ12, df1dQ13,        0,        0},
										 {df2dp1,   df2dp2,       0,        0,     df2dQ12,       0,        0,        0},
										 {0,        df3dp2,       0,   df3dp4,       0,           0,  df3dQ24,        0},
										 {df4dp1,       0,     df4dp3,      0,       0,     df4dQ13,        0,        0},
										 {0,            0,     df5dp3, df5dp4,       0,           0,        0,  df5dQ34},
										 {0,            0,        0,   df6dp4,       0,           0,  df6dQ24,  df6dQ34},
										 {0,        df7dp2,       0,        0,     df7dQ12,       0,  df7dQ24,        0},
										 {0,            0,     df8dp3,      0,       0,     df8dQ13,        0,  df8dQ34}

			};
			vector<vector<double>>  inverseJ = InverseMatrix(jm);

			if (inverseJ[0][0] == error)
			{
				ERROR[0] = error;
				ERROR[1] = (double)ERROR_TYPE::INVERSE;                            //��� ������ ������ ��� ������ ����� � �������������� �������
				return ERROR;
			}

			vector<double> f = { -f1, -f2, -f3, -f4, -f5, -f6, -f7, -f8 };         //������ �������

			vector<double> dp(N, 0);                                               //Xk+1 = Xk - F * J^-1

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					dp[i] += stepK * inverseJ[i][j] * f[j];
				}
			}

			for (int i = 0; i < N; i++)                                            //���������� �������� ����������
			{
				x[i] += dp[i]; 
			}

			x[8] = x[4] + x[5];
			x[9] = x[6] - x[4];
			x[10] = x[7] - x[5];
			x[11] = x[6] + x[7];

			LOG("��������: # f1 # dx1", 3, (double)it, f[1], dp[1]);

			for (int i = 0; i < N; i++)
			{
				sumF += abs(f[i]);
			}

			if (sumF < eps)          //�������� ������� ��������� �������� �� �������
			{
				return x;
			}
		}
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
	void ShowCompressedMatrix(const void* matrix, DataType dataType)
	{
		auto m = (CompressedMatrix4V*)matrix;	
		cout << endl << "��������:";
		for (double value : m->values)
		{
			cout << "\t" << value <<" ";
		}
		cout << endl << endl;

		cout << "�������: ";
		for (int column : m->columns)
		{
			cout <<"\t" << column <<" ";
		}
		cout << endl;

		cout << "������: ";
		for (int rows : m->rows)
		{
			cout << "\t" << rows;
		}
		cout << endl;

		if (dataType == DataType::CSRAA)
		{
			cout << "C����� ������: ";
			for (int rows : m->nextRow)
			{
				cout << "\t" << rows << " ";
			}
			cout << endl;
		}
		cout << endl << line << endl<< endl;
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
				cout << right  << setprecision(5) << element->value <<", ";
				cout << left << "i=" << element->row;
				cout << left << " j=" << element->column << ",  ";
				dataType == DataType::KNUTHRowsList ? element = element->nextInRow : element = element->nextInColumn;
			}
		}
		cout << endl << endl << line << endl << endl;
	}

	//������� �� ������� �������������� ������ x
	void ShowResultX(vector<double> x)                                              //�������� ��������� � ������ x    
	{
		cout << "\t\t\t������ �������� � ����� � �������� ��������: \n\n";
		for (int i = 0; i < N; i++)
		{
			cout << i + 1 << ".\t" << signatureXN[i] << " = " << "\t" << x[i] << endl;
		}
		cout << endl << line << endl<< endl;
	}

	//������� ��� ������ ������� �� �����
	void ShowMatrix(const vector<vector<double>> matrix)
	{
		int N = (int)matrix.size();
		int M = (int)matrix[0].size();
		cout << endl;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				cout << setw(10) << fixed << setprecision(2) << matrix[i][j];
			}
			cout << endl;
		}
		cout << scientific;
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
			return k;
		default:
			cout << line << "\n\t�� ���������� ������ ������!" << endl << line;
			return (int)symbol;
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
			ERROR[0] = error;
			ERROR[1] = (double)ERROR_TYPE::FILE;
		}
		return matrix;
	}

	//��������� �������������� ������  � ������� JSON � ����, ���� - �������� ���������
	void WriteJSON(const string& filename,  const vector<double>& x)
	{
		string jsonData = ConvertVectorToJson(signatureXN, x);
		ofstream outputFile(filename);

		if (outputFile.is_open())
		{
			outputFile << jsonData;
			outputFile.close();
		}
		else
		{
			cout << line << "\n\t���������� ������� ���� ��� ������ : " << filename << endl << line;
			ERROR[0] = error;
			ERROR[1] = (double)ERROR_TYPE::FILE;
		}
	}

	//�������������� �������� � ������ JSON
	string ConvertVectorToJson(const vector<string>& title, const vector<double>& x)
	{
		string json = "{\n";

		for (size_t i = 0; i < title.size(); i++)
		{
			json += "    \"" + title[i] + "\": " + to_string(x[i]);                  //�������� �������� ��������� � ��� �������� � ������ JSON 

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
	HydraulicNet(vector<vector<double>>& aMatrix, int numOfParameters, double resistK) : a{ aMatrix }, N{ numOfParameters }, k{ resistK }
	{
		numNonZero = CountNonZeroElements(aMatrix);
	}
};



//��������; ��������� ����� ���������� �-���, � �������� ���������� ��������� ������ � ��� �����, ����� ���������� �������� ����������
double GetCalculationRuntime(HydraulicNet, vector<double>(HydraulicNet::* calculation)(vector<double>&), vector<double>, vector<double>*);

//��������; ��������� ��������� ����������� ��������
vector<double> Run(DataType, Task, HydraulicNet, const vector<vector<double>>,  vector<double>, const vector<double>, bool save = true);

//
//************************************* MAIN *******************************************
//
int main()
{
	setlocale(LC_ALL, "");

	Pipe pipe;
	double k = pipe.k;                          //��� ������� �������� ��������� ����� ��������� �����. ������������ - �������������� �������������

	//signatureX' :: i :{ "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34",       "Q01", "Q02", "Q03", "Q04" }
	//                       0    1     2     3     4      5      6       7            8     9      10    11             // 

	vector<vector<double>> a =
	{
		{1,  0, 0, 0,   0, 0, 0,  0,    k, 0, 0,   0},
		{0,  1, 0, 0,   0, 0, 0,  0,    0, k, 0,   0},
		{0,  0, 1, 0,   0, 0, 0,  0,    0, 0, k,   0},
		{0,  0, 0, 1,   0, 0, 0,  0,    0, 0, 0,  -k},

		{0,  0, 0, 0,   1, 1, 0,  0,   -1, 0, 0,   0},
		{0,  0, 0, 0,   0, 1, 0, -1,    0, 0, 1,   0},
		{0,  0, 0, 0,   0, 0, 1,  1,    0, 0, 0,  -1},
		{0,  0, 0, 0,   1, 0,-1,  0,    0, 1, 0,   0},

		{-1, 1, 0, 0,   k, 0, 0,  0,    0, 0, 0,   0},
		{0, -1, 0, 1,   0, 0, k,  0,    0, 0, 0,   0},
		{-1, 0, 1, 0,   0, k, 0,  0,    0, 0, 0,   0},
		{ 0, 0, 1,-1,   0, 0,  0, k,    0, 0, 0,   0},
	};

	HydraulicNet net(a, N, k);
	vector<double> b;
	vector<double> pInput(PInputCount, 0);                                                      //������ ������� ��������, 4 * P���

	Mode mode = Mode::DEBUG;

	//************************************ ��� ������ !!!!! ******************************      //��� ������ �������� Po 1..4

	(mode == Mode::RELEASE) ? pInput = net.InputValues(PInputCount) : pInput = { 150, 130, 100, 50 };


	b = { pInput[0], pInput[1], pInput[2], pInput[3], 0, 0, 0, 0, 0, 0, 0, 0 };
	vector<double> x(N, 0.0);                                                                   //�������� ��������� � ������ x

	Run(DataType::ALL, Task::All, net, a,b,x);																//������ ����������� ��������
	
	cout << endl << line << endl << "Error: " << ERROR[1];
	return 0;
}

//��������� ��������� ����������� ��������, save - ���� ���������� � JSON ����
vector<double> Run (DataType dataType, Task task, HydraulicNet net, const vector<vector<double>> a, vector<double> b, vector<double> x,  bool save)
{
	vector<double> funcRuntime((int)(Task::MaxTask) - 1, 0);
	bool runAll{ false };
	cout << endl << " �������������� �������� ������� � ���� ";
	if (dataType == DataType::ALL || task == Task::All)
	{
		runAll = true;																			//���������� �� ����� ������ �����, CLA ������ ���� ������ switch
	}

	switch (dataType) 
	{
	case DataType::CSRAA :
		cout << " CompressedSparseRowAdapted, 4 ������� CLA + CRSA: \n" << endl;											
		net.CSRAA = net.ConvertToCSRAAdapted(a);
		net.ShowCompressedMatrix(&net.CSRAA, DataType::CSRAA);
		if (!runAll)
			break;

	case DataType::ALL :

	case DataType::KNUTHRowsList :
		cout << "������ ����� : \n" << endl;													    //������ ������ ������ �������, ������ ������ ������� �������
		net.ConvertToKnuthLists(a);
		net.ShowKnuthList(net.KNUTHrows, DataType::KNUTHRowsList);
		net.ShowKnuthList(net.KNUTHcolumns, DataType::KNUTHColumnsList);
		if (!runAll)
			break;

	case DataType::CSRA :
		cout << " CompressedSparseRow, 3 �������, ���� - ������������ ����� �� 0 ���������: \n\n";	
		net.�SRA = net.ConvertToCSRA(a);
		net.ShowCompressedMatrix(&net.�SRA, DataType::CSRA);
		if (!runAll)
			break;

	case DataType::CLA :
		cout << " Coordinate List, 3 �������: \n\n";														
		net.CLA = net.ConvertToCLA(a);
		net.ShowCompressedMatrix(&net.CLA, DataType::CLA);
		if(!runAll) 
			break;
	}

	if (x[0] != error)
	{
		double z0{ 1.5 }, z{ 0 };																		//���������� ��� ������ ��������� � ������ �������
		
		switch (task)
		{
		case Task::All :

		case Task::GaussCLA:
			funcRuntime[(int)Task::GaussCLA] = GetCalculationRuntime(net, &HydraulicNet::GaussCLA, b, &x);
			cout << " ������� ������� �������� ��������� ������� ������.\n\n";
			net.ShowResultX(x);
			if (!runAll)
				break;

		case Task::Gauss:
			funcRuntime[(int)Task::Gauss] = GetCalculationRuntime(net, &HydraulicNet::Gauss, b, &x);
			cout << " ������� ������� �������� ��������� ������� ������.\n\n";
			net.ShowResultX(x);
			if (!runAll)
				break;

		case Task::Newton1 : 
			z = net.Newton(equation1, z0);
			cout << " ������� ���������   x^2 + x - 6   ������� �������.\n";
			cout << " ������ ���������:   " << z << endl << endl << line;
			if (!runAll)
				break;

		case Task::NewtonHardCode4 :
			x[8] = b[0];																				//x � ������ ������ - ������ ��������� �������� � p0
			x[9] = b[1];
			x[10] = b[2];
			x[11] = b[3];
			funcRuntime[(int)Task::NewtonHardCode4] = GetCalculationRuntime(net, &HydraulicNet::NewtonHardCode, x, &x);
			cout << " ������� ������� ���������� ��������� ������� �������.";
			x = net.NewtonHardCode(x);
			net.ShowResultX(x);
		}

		cout << "\n ����� �������� ������� � 3 ���������: " << funcRuntime[(int)Task::GaussCLA] << ", ���.";
		cout << "\n ����� �������� ������� � ����������� �������� ��� ������: " << funcRuntime[(int)Task::Gauss] << ", ���.";
		cout << "\n ����� �������� �������� ������� ���������: " << funcRuntime[(int)Task::NewtonHardCode4] << ", ���.";
		cout << endl << line << endl;
	}
		if (save)
		{
			net.WriteJSON(filenameResultJSON, x);
		}

	return x;
};

//
//��������� ����� ���������� �-���, � �������� ���������� ��������� ������ � ��� �����, ����� ���������� ������ ����������
//
double GetCalculationRuntime(HydraulicNet net, vector<double>(HydraulicNet::* calculation)(vector<double>&), vector<double> matrixB, vector<double>* x)
{
	auto start = chrono::high_resolution_clock::now();                                          // ���������� ����� ������ ���������� �������

	(*x) = (net.*calculation)(matrixB);

	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(end - start);                    // ��������� ����� ���������� ������� � �������������
	return duration.count();
}