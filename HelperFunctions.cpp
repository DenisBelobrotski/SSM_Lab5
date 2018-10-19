#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <float.h>


using namespace std;


double* calcConvertedSLE(double** matrixA, double* vectorF, int systemSize, int chainLength, int chainNum);
double calcNext(double** matrixCoefs, double* freeVector, double **markovMatrix, double* markovVector,
                double* helperVector, int* trajectory, int trajectorySize);
bool checkCoefs(double** coefs, int size);
double calcResidual(double* calcResult, double* exactResult, int size);

void printMatrix(double** mtr, int lines, int columns);
void printVector(double* vector, int size);


void calcSolutions() {
    const int chainLengthMin = 5;
    const int chainLengthStep = 1;
    const int chainLengthStepsNum = 45;
    const int chainNumMin = 1000;
    const int chainNumStep = 1000;
    const int chainNumStepsNum = 20;

    int systemSize = 0;
    double **matrixA = nullptr;
    double *vectorF = nullptr;
    double *exactSolution = nullptr;

    double tmp = 0;
    ifstream fin("input.txt");

    if (fin.is_open()) {
        fin >> systemSize;

        matrixA = new double *[systemSize];
        vectorF = new double[systemSize];
        exactSolution = new double[systemSize];

        for (int i = 0; i < systemSize; i++) {
            matrixA[i] = new double[systemSize];
            for (int j = 0; j < systemSize; j++) {
                fin >> tmp;
                matrixA[i][j] = tmp;
            }
        }

        for (int i = 0; i < systemSize; i++) {
            fin >> tmp;
            vectorF[i] = tmp;
        }

        for (int i = 0; i < systemSize; i++) {
            fin >> tmp;
            exactSolution[i] = tmp;
        }
    }

    fin.close();

    if (systemSize)
    {
        ofstream fout("output.txt");

        fout << "length;num;residual" << endl;

        auto newMatrixA = new double* [systemSize];
        for (int i = 0; i < systemSize; i++) {
            newMatrixA[i] = new double[systemSize];
        }
        for (int i = 0; i < systemSize; i++)
        {
            for (int j = 0; j < systemSize; j++)
            {
                newMatrixA[i][j] = (i == j) ? -(matrixA[i][j] - 1.0) : -matrixA[i][j];
            }
        }

        if (!checkCoefs(newMatrixA, systemSize))
        {
            printf("Incorrect matrix");
            return;
        }

        int currentChainLength = chainLengthMin;
        int currentChainNum;
        auto residual = DBL_MIN;
        double* result = nullptr;

        for (int i = 0; i < chainLengthStepsNum; i++)
        {
            currentChainNum = chainNumMin;

            for (int j = 0; j < chainNumStepsNum; j++)
            {
                result = calcConvertedSLE(newMatrixA, vectorF, systemSize, currentChainLength, currentChainNum);
                residual = max(residual, calcResidual(result, exactSolution, systemSize));
                fout << currentChainLength << ";" << currentChainNum << ";" << residual << endl;
                currentChainNum += chainNumStep;
                delete[] result;
            }
            currentChainLength += chainLengthStep;

            printf("Iteration: %d of %d\n", i + 1, chainLengthStepsNum);
        }
        printf("\nMax residual: %lf\n", residual);

        for (int i = 0; i < systemSize; i++)
        {
            delete[] newMatrixA[i];
        }
        delete[] newMatrixA;

        fout.close();
    }

    if (matrixA != nullptr)
    {
        for (int i = 0; i < systemSize; i++)
        {
            delete[] matrixA[i];
        }
    }
    delete[] matrixA;
    delete[] vectorF;
    delete[] exactSolution;
}


double* calcConvertedSLE(double** matrixA, double* vectorF, int systemSize, int chainLength, int chainNum)
{
    int size = systemSize;
    auto result = new double[size];
    auto markovMtr = new double* [size];
    for (int i = 0; i < size; i++)
    {
        markovMtr[i] = new double[size];
    }
    auto markovVector = new double[size];
    auto helperVector = new double[size];
    int** trajectories = new int* [chainNum];
    for (int i = 0; i < chainNum; i++)
    {
        trajectories[i] = new int[chainLength];
    }

    for (int i = 0; i < size; i++)
    {
        markovVector[i] = 1.0 / size;

        for (int j = 0; j < size; j++)
        {
            markovMtr[i][j] = 1.0 / size;
        }
    }

    for (int i = 0; i < chainNum; i++)
    {
        for (int j = 0; j < chainLength; j++)
        {
            trajectories[i][j] = rand() % size;
        }
    }

    for (int i = 0; i < chainNum; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                helperVector[k] = 0.0;
            }
            helperVector[j] = 1.0;
            result[j] += calcNext(matrixA, vectorF, markovMtr, markovVector, helperVector, trajectories[i], chainLength);
        }
    }

    for (int i = 0; i < size; i++)
    {
        result[i] /= chainNum;
    }

    return result;
}


double calcNext(double** matrixCoefs, double* freeVector, double **markovMatrix, double* markovVector,
                double* helperVector, int* trajectory, int trajectorySize)
{
    int indexPrev;
    int indexCur;
    double result;
    double weight;

    indexCur = trajectory[0];
    weight = (markovVector[indexCur] > 0.0) ? helperVector[indexCur] / markovVector[indexCur] : 0.0;
    result = weight * freeVector[indexCur];
    for (int i = 1; i < trajectorySize; i++)
    {
        indexPrev = trajectory[i - 1];
        indexCur = trajectory[i];
        weight *= (markovMatrix[indexPrev][indexCur] > 0.0) ? matrixCoefs[indexPrev][indexCur] / markovMatrix[indexPrev][indexCur] : 0.0;

        result += weight * freeVector[indexCur];
    }

    return result;
}


double calcResidual(double* calcResult, double* exactResult, int size)
{
    auto result = DBL_MIN;
    for (int i = 0; i < size; i++)
    {
        result = max(result, abs(calcResult[i] - exactResult[i]));
    }

    return result;
}


bool checkCoefs(double** coefs, int size)
{
    double sum;

    for (int i = 0; i < size; i++)
    {
        sum = 0.0;
        for (int j = 0; j < size; j++)
        {
            sum += abs(coefs[i][j]);
        }

        if (sum >= 1.0)
        {
            return false;
        }
    }

    return true;
}


void printMatrix(double** mtr, int lines, int columns)
{
    for (int i = 0; i < lines; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            printf("%lf ", mtr[i][j]);
        }
        printf("\n");
    }
}


void printVector(double* vector, int size)
{
    for (int i = 0; i < size; i++)
    {
        printf("%lf ", vector[i]);
    }
    printf("\n");
}
