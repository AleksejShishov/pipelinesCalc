#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdarg.h>

#include "gauss.h"

using namespace std;

const string line(100, '=');                    //строка разделитель для вывода на консоль
const string filenameM = "matrix.txt";          //матрицу можно загрузить из файла
const string filenameResultJSON = "result.json";//расчетные параметры сохраняются в формате JSON 

const double error = 1.7e308;					//флаг, что вернулась ошибка ERROR_TYPE в векторе ERRROR или ERROR 2
vector<double> ERROR = { 0 , 0 };
vector<vector<double>> ERROR2 = { {ERROR}, {0, 0} };

enum class ERROR_TYPE							
{
	GAUSS = 1,
	GAUSSCLA,
	GAUSSCSRAA,
	GAUSSKNUTH,
	NEWTON1,
	NEWTONHARDCODE,
	INVERSE,
	FILE,
	Max
};


//режимы для отладки
enum class Mode
{
	DEBUG = false,
	RELEASE = true
};

const string LOGfile = "param.log";

//логи
void LOG(string comment, const int num, ...)
{
	ofstream outputFile(LOGfile);

	if (outputFile.is_open())
	{
		string result{ comment };
		va_list factor;
		va_start(factor, num);
		for (int i = 0;i < num; i++)
		{
			result += to_string(va_arg(factor, double));						// получаем значение текущего параметра типа int
			result += " # ";
		}
		va_end(factor);

		outputFile << result << endl;
		outputFile.close();
	}
	else
	{
		cout << "\n\tНевозможно открыть файл для записи Log'ов: " << LOGfile << endl << line;
		ERROR[0] = error;
		ERROR[1] = (double)ERROR_TYPE::FILE;
	}
}

//enum задач
enum class Task
{
	Gauss = 0,
	GaussCLA,
	GaussCSRAA,
	GaussKnuth,
	Newton1,
	NewtonHardCode4,
	MaxTask = 6
};

const int PInputCount{ 4 };                     //число входных трубопроводов, 4 по умолчанию
const int N{ 12 };                              //число расчетных параметров, 12 по умолчанию, передаётся в конструкторе в модель

//
// СТРУКТУРЫ ДАННЫХ ДЛЯ ХРАНЕНИЯ РАЗРЕЖЕННОЙ МАТРИЦЫ
// 
//enum используемых структур данных
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

//матрица в виде связных списков по схеме Кнута, Node - элемент списка
struct Node {
	double value;
	int column;
	int row;
	Node* nextInRow;
	Node* nextInColumn;
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
	double k = ksi * lam * water.nu * reLam;    //функция проводимости - гидр. сопротивления для ламинарного режима
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

	vector<string> signatureXN { "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34", "Q01", "Q02", "Q03", "Q04" }; //расчетные параметры, подписи для Ньютона

	//
	//********************************* КОНВЕРТАЦИЯ РАЗРЕЖЕННОЙ МАТРИЦЫ *************
	//
	//
	void ConvetrtToVector(vector<vector<double>> matrix, vector<double>& values, vector<int>& rows, vector<int>& columns, vector<int>& nextRow)
	{
		values.resize(numNonZero);
		rows.resize(numNonZero);
		columns.resize(numNonZero);
		nextRow.resize(N + 1);
		int c = 0;                              //кумулятивное число не нулевых элементов в строке
		nextRow.push_back(c);
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
			nextRow.push_back(c);
		}	
		for (double c : values)
		cout << "   " << c;
	}

	//конвертирует разреженную матрицу в формат 3-х векторов: значения, координаты i, координаты j
	CompressedMatrix ConvertToCLA(const vector<vector<double>>& matrix)
	{
 		vector<double> values;
		vector<int> columns;
		vector<int> rows;
		CompressedMatrix CLAVectors;
		vector<int> none;

//		ConvetrtToVector(matrix, CLA.values,  CLA.rows, CLA.columns, none);
		

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
	//считает матрицу А и вектор В по методу Гаусса    
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

	//считает матрицу, представленную в виде Coordinate List и вектор В по методу Гауса 
	vector<double> GaussCLA(vector<double>& b)
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
				cout << "\n!\n";
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
				cout << endl << " change i to moain: " << i << " " << mainRow;
				mainRow = i;
			}

			bool zeroAtCurrentRow;                              //признак, что в строке есть 0 элемент, который после вычитания строк !=0
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
						if (CLA.rows[l] == i && CLA.columns[l] >= i)
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

					CLA.values.erase(CLA.values.begin() + j);                        //не обнуляем нижние элементы, а удаляем из списка
					CLA.rows.erase(CLA.rows.begin() + j);
					CLA.columns.erase(CLA.columns.begin() + j);
					--numNonZero;

					vector<int>::iterator it;                                        //добавляем новые не 0 элементы, равые 0 - coef*a[i][l]
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
		for (int i = N - 1; i >= 0; i--)                                             //обратный ход 
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
			x[i] = (b[i] - sum) / CLA.values[diagonal];                             //результирующий вектор 
		}
		return x;
	}

	//считает матрицу А и вектор В по методу Гаусса, матрица a в виде списка не 0 элементов с доп. вектором: кумулятивное кол-во не 0 элементов в строке   
	vector<double> GaussCRSAA(vector<double>& b)
	{
		ERROR[0] = error;
		ERROR[1] = (double)ERROR_TYPE::GAUSSCSRAA;
		return ERROR;
	}

	//считает матрицу, представленную в варианте Кнута - списки начала строк и столбцов
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

	//функция для нахождения обратной матрицы
	vector<vector<double>> InverseMatrix(vector<vector<double>> a)
	{
		int N{ (int)a.size() };
		vector<vector<double>> extended(N, vector<double>(N * 2, 0.0));
		vector<vector<double>>inverse(N, vector<double>(N, 0));
		for (int i = 0; i < N; i++)
		{													  //добавляем единичную матрицу справа
			for (int j = 0; j < N; j++)
			{
				extended[i][j] = a[i][j];
			}
			extended[i][i + N] = 1;
		}

		for (int i = 0; i < N; i++)						      //приводим расширенную матрицу к диагональному виду, проверим сводимость
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
					ERROR2[0][0] = error;
					ERROR2[0][1] = (double)ERROR_TYPE::INVERSE;                     
					return ERROR2;
				}
				swap(extended[i], extended[j]);
			}
		
			for (int j = i + 1; j < N; j++)                    //прямой ход Гаусса
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

		for (int i = 0; i < N; i++)								 //нормализуем строки расширенной матрицы, диагональ левой матрицы = 1
		{
			double coef = extended[i][i];
			for (int j = i; j < N * 2; j++)
			{
				extended[i][j] /= coef;
			}
		}

		for (int i = 0; i < N; i++)								 //извлекаем обратную матрицу из расширенной матрицы
		{
			for (int j = 0; j < N; j++)
			{
				inverse[i][j] = extended[i][j + N];
			}
		}

		return inverse;
	}

	//для произвольной функции можно использовать численный метод deltaY / deltaX
	double GetDiff(double (*equation)(double), double x, double dp = 1.0e-6)
	{
		double derivative = (equation(x + dp) - equation(x - dp)) / (2 * dp);
		return derivative;
	}
	double GetDifff(double (*equation)(double), double x, double dp = 1.0e-6)
	{
		return 2*x + 1;
	}


	//считает Ньютоном одно уравнение
	double Newton(double (*equation)(double), double x, double epsilon = 1e-06, int it = 100)
	{
		for (int i = 0; i < it; i++) 
		{
			double f = equation(x);
			double df = GetDiff(equation, x, epsilon);
			double dx = f / df;
			x = x - dx;
			LOG("шаг итерации Ньютона1",1,dx);
			if (abs(f) < epsilon)
			{
				return x;
			}
		}
		return x;
	}

	//TODO            typedef double (*FunctionArr[])(double);
	//считает СНАУ методом Ньютона; ограничено число итераций: 100, но можно перегрузить в параметрах эпсилон - приближение и число итераций
	vector<double> NewtonHardCode(const vector<double>& p0, vector<double>& x, double eps = 1e-6, int itMax = 100, int N = 8)
	{
		int iterations{ 0 };						             //ограничим число итераций: 100, но можно перегрузить в параметрах эпсилон - приближение и iter

		for (int it = 0; it < itMax; it++)
		{													     //УРАВНЕНИЯ В виде функций в Gauss.h
			x[8] = x[4] + x[5];
			x[9] = x[6] - x[4];
			x[10] = x[7] - x[5];
			x[11] = x[6] + x[7];

			double f1 = p0[0] - x[0] - x[4] * x[4] * k;			 //вычисление значений функций и их производных
			double f2 = x[0] - x[1] - x[5] * x[5] * k;
			double f3 = x[1] - x[3] - x[7] * x[7] * k;
			double f4 = x[0] - x[2] - x[11] * x[11] * k;
			double f5 = x[2] - x[3] - x[9] * x[9] * k;
			double f6 = x[3] - p0[3] - x[8] * x[8] * k;
			double f7 = p0[1] - x[1] - x[6] * x[6] * k;
			double f8 = p0[2] - x[2] - x[10] * x[10] * k;
			/*
						 double f9  = p0[0];                     //добавим уравнения, чтобы не менять размер матрица Якоби
						 double f10 = p0[1];
						 double f11 = p0[2];
						 double f12 = p0[3];
			*/
			double df1dp1 = -1;									 //производные
			double df1dQ01 = -2 * x[4] * k;
			double df2dp2 = -1;
			double df2dQ12 = -2 * x[5] * k;
			double df3dp2 = 1;
			double df3dp4 = -1;
			double df3dQ24 = -2 * x[7] * k;
			double df4dp1 = 1;
			double df4dp3 = -1;
			double df4dQ13 = -2 * x[11] * k;
			double df5dp3 = 1;
			double df5dp4 = -1;
			double df5dQ34 = -2 * x[9] * k;
			double df6dp4 = 1;
			double df6dQ04 = -2 * x[8] * k;
			double df7dp2 = -1;
			double df7dQ02 = -2 * x[6] * k;
			double df8dp3 = -1;
			double df8dQ03 = -2 * x[10] * k;

			//                                           ***     Матрица Якоби    ***
			//signatureX' :: i :{ "P1", "P2", "P3", "P4", "Q12", "Q13", "Q24", "Q34",       "Q01", "Q02", "Q03", "Q04" }
			//                       0    1     2     3     4      5      6       7            8     9      10    11             // 
			//            :: j  :{f1 ..f12}
			// 

			vector<vector<double>> jm =
			{
										 {df1dp1,       0,        0,        0,       0,       0,     0,        0},
										 {0,        df2dp2,       0,        0,       df2dQ12, 0,     0,        0},
										 {0,        df3dp2,       0,   df3dp4,       0,       0,     df3dQ24,  0},
										 {df4dp1,       0,     df4dp3,      0,       0, df4dQ13,     0,        0},
										 {0,            0,     df5dp3, df5dp4,       0,       0,     0,  df5dQ34},
										 {0,            0,        0,   df6dp4,       0,       0,     0,        0},
										 {0,        df7dp2,       0,        0,       0,       0,     0,        0},
										 {0,            0,     df8dp3,      0,       0,       0,     0,        0}

			};
			vector<vector<double>>  inverseJ = InverseMatrix(jm);

			if (inverseJ[0][0] == error)
			{
				ERROR[0] = error;
				ERROR[1] = (double)ERROR_TYPE::INVERSE;                            //код ошибки Гаусса для Ньтона зашит в результирующем векторе
				return ERROR;
			}

			vector<double> f = { -f1, -f2, -f3, -f4, -f5, -f6, -f7, -f8 };         //вектор функций

			vector<double> dp(N, 0);                                               //jMatrix * dp = f

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					dp[i] += inverseJ[i][j] * f[j];
				}
			}

			for (int i = 0; i < N; i++)                                            //обновление значений переменных
			{
				x[i] += dp[i]; 
			}

	//		for ()
			if (abs(dp[0]) < eps && abs(dp[1]) < eps && abs(dp[2]) < eps)          //проверка условия окончания итераций по эпсилон
			{
				break;
			}
		}
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
	void ShowCompressedMatrix(const void* matrix, DataType dataType)
	{
		auto m = (CompressedMatrix4V*)matrix;	
		cout << endl << "Значения:";
		for (double value : m->values)
		{
			cout << "\t" << value <<" ";
		}
		cout << endl << endl;

		cout << "Столбцы: ";
		for (int column : m->columns)
		{
			cout <<"\t" << column <<" ";
		}
		cout << endl;

		cout << "Строки: ";
		for (int rows : m->rows)
		{
			cout << "\t" << rows;
		}
		cout << endl;

		if (dataType == DataType::CSRAA)
		{
			cout << "Cжатые сторки: ";
			for (int rows : m->nextRow)
			{
				cout << "\t" << rows << " ";
			}
			cout << endl;
		}
		cout << endl << line << endl<< endl;
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
				cout << right  << setprecision(5) << element->value <<", ";
				cout << left << "i=" << element->row;
				cout << left << " j=" << element->column << ",  ";
				dataType == DataType::KNUTHRowsList ? element = element->nextInRow : element = element->nextInColumn;
			}
		}
		cout << endl << endl << line << endl << endl;
	}

	//выводит на консоль результирующий вектор x
	void ShowResultX(vector<double> x)                                              //собираем результат в вектор x    
	{
		cout << "\t\t\tРасчет давлений в узлах и объёмных расходов: \n\n";
		for (int i = 0; i < N; i++)
		{
			cout << i + 1 << ".\t" << signatureXN[i] << " = " << "\t" << x[i] << endl;
		}
		cout << endl << line << endl<< endl;
	}

	//функция для вывода матрицы на экран
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
			return k;
		default:
			cout << line << "\n\tНе правильный формат данных!" << endl << line;
			return (int)symbol;
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
			ERROR[0] = error;
			ERROR[1] = (double)ERROR_TYPE::FILE;
		}
		return matrix;
	}

	//сохраняет результирующий вектор  в формате JSON в файл, ключ - название параметра
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
			cout << line << "\n\tНевозможно открыть файл для записи : " << filename << endl << line;
			ERROR[0] = error;
			ERROR[1] = (double)ERROR_TYPE::FILE;
		}
	}

	//преобразования векторов в формат JSON
	string ConvertVectorToJson(const vector<string>& title, const vector<double>& x)
	{
		string json = "{\n";

		for (size_t i = 0; i < title.size(); i++)
		{
			json += "    \"" + title[i] + "\": " + to_string(x[i]);                  //собираем название параметра и его значение в строку JSON 

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
	HydraulicNet(vector<vector<double>>& aMatrix, int numOfParameters, double resistK) : a{ aMatrix }, N{ numOfParameters }, k{ resistK }
	{
		numNonZero = CountNonZeroElements(aMatrix);
	}
};



//прототип; вычисляет время выполнения ф-ции, в качестве параметров передатся объект и его метод, время выполнения которого измеряется
double GetCalculationRuntime(HydraulicNet, vector<double>(HydraulicNet::* calculation)(vector<double>&), vector<double>, vector<double>*);

//прототип; выполняет обработку результатов расчётов
vector<double> Run(DataType, HydraulicNet, const vector<vector<double>>,  vector<double>, const vector<double>, bool save = true);

//
//************************************* MAIN *******************************************
//
int main()
{
	setlocale(LC_ALL, "");

	Pipe pipe;
	double k = pipe.k;                          //для системы линейных уравнений сразу вычисляем коэфф. проводимости - гидравлическое сопротивление

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
	vector<double> pInput(PInputCount, 0);                                                      //вводим входные давления, 4 * Pвхю

	Mode mode = Mode::DEBUG;

	//************************************ ДЛЯ ТЕСТОВ !!!!! ******************************      //для тестов харкодим Po 1..4

	(mode == Mode::RELEASE) ? pInput = net.InputValues(PInputCount) : pInput = { 150, 130, 100, 50 };


	b = { pInput[0], pInput[1], pInput[2], pInput[3], 0, 0, 0, 0, 0, 0, 0, 0 };
	vector<double> x(N, 0.0);                                                                   //собираем результат в вектор x

	ERROR = { 0, 0 };																			//сброс ошибок
	ERROR2 = { {0, 0},{0,0} };																	

	Run(DataType::ALL, net, a,b,x);																//ЗАПУСК НЕОБХОДИМЫХ РАСЧЕТОВ
	
	cout << endl << line << endl << "Error: " << ERROR[1];
	return 0;
}

//выполняет обработку результатов расчётов, save - флаг сохранения в JSON файл
vector<double> Run (DataType dataType, HydraulicNet net, const vector<vector<double>> a, vector<double> b, vector<double> x,  bool save)
{
	vector<double> funcRuntime((int)(Task::MaxTask) - 1, 0);
	bool runAll{ false };
	cout << endl << " Преобразование исходной матрицы к виду ";
	if (dataType == DataType::ALL)
	{
		dataType = DataType::CLA;	
		runAll = true;																			//пробежимся по всему списку задач, CLA должна быть сверху switch
	}
	switch (dataType) 
	{
	case DataType::CLA :
		cout << " Coordinate List, 3 вектора: \n\n";														
		net.CLA = net.ConvertToCLA(a);
		funcRuntime [(int)Task::GaussCLA] = GetCalculationRuntime(net, &HydraulicNet::GaussCLA, b, &x);
		if(!runAll) 
			break;

	case DataType::CSRA :
		cout << " CompressedSparseRow, 3 вектора, ряды - кумулятивное число не 0 элементов: \n\n";	
		net.СSRA = net.ConvertToCSRA(a);
		net.ShowCompressedMatrix(&net.СSRA, DataType::CSRA);
		if (!runAll)
			break;

	case DataType::CSRAA :
		cout << " CompressedSparseRowAdapted, 4 вектора CLA + CRSA: \n" << endl;											
		net.CSRAA = net.ConvertToCSRAAdapted(a);
		net.ShowCompressedMatrix(&net.CSRAA, DataType::CSRAA);
//		funcRuntime [(int)Task::GaussCSRAA] = GetCalculationRuntime(net, &HydraulicNet::GaussCRSAA, b, &x); //выводим время выполнения метода Гауса
		net.ShowCompressedMatrix(&net.CSRAA, DataType::CSRAA);
		if (!runAll)
			break;

	case DataType::KNUTHRowsList :
		cout << "список Кнута : \n" << endl;													    //список начала каждой строчки, список начала каждого столбца
		net.ConvertToKnuthLists(a);
		net.ShowKnuthList(net.KNUTHrows, DataType::KNUTHRowsList);
		net.ShowKnuthList(net.KNUTHcolumns, DataType::KNUTHColumnsList);
//		funcRuntime[(int)Task::GaussKnuth] = GetCalculationRuntime(net, &HydraulicNet::GaussKnuth, b, &x);
		if (!runAll)
			break;
	case DataType::MATRIX :
		funcRuntime [(int)Task::Gauss] = GetCalculationRuntime(net, &HydraulicNet::Gauss, b, &x);	//матрица без сжатия	
	}

	if (x[0] != error) 
	{
		net.ShowResultX(x);
		cout << "\n Время расчетов Гауссом с 3 векторами: " << funcRuntime[(int)Task::GaussCLA] << ", мкс.";
		cout << "\n Время расчетов Гауссом с разреженной матрицей без сжатия: " << funcRuntime[(int)Task::Gauss] << ", мкс.";
		cout << endl << line << endl;
		if (save)
		{
			net.WriteJSON(filenameResultJSON, x);
		}

		cout << endl << line << endl;
		cout << " Расчет уравнения z^2 - 1 методом Ньютона.\n\n";

		double z0{1.5};
		double z = net.Newton(equation1, z0);
		cout << z << endl;

		vector<double> pInput = {b[0], b[1], b[2], b[3]};
		
		x = net.NewtonHardCode(pInput, x);

		net.ShowResultX(x);
	}
	else
	{
		cout << endl << line << "\nОшибка номер: " << x[1] << ", см. enum ERROR_TYPE.\n";
	}
	return b;
};

//
//вычисляет время выполнения ф-ции, в качестве параметров передатся объект и его метод, время выполнения метода измеряется
//
double GetCalculationRuntime(HydraulicNet net, vector<double>(HydraulicNet::* calculation)(vector<double>&), vector<double> matrixB, vector<double>* x)
{
	auto start = chrono::high_resolution_clock::now();                                          // запоминаем время начала выполнения функции

	(*x) = (net.*calculation)(matrixB);

	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(end - start);                    // вычисляем время выполнения функции в миллисекундах
	return duration.count();
}