#include "CMatrix.h"


//���캯������ʼ��ʱ�Զ�����
CMatrix::CMatrix(int m, int n)
{
	m_nRows = m;
	m_nColumns = n;
	//��Ҫ��ֹ�˷�������˴�����
	this->m_lpBuf = new double[m*n];
	if (this->m_lpBuf == 0)
	{
		printf("Error,Menory allocation failure!");
		exit(0);
	}
	memset(m_lpBuf, 0, m*n * sizeof(double));
}
//�������캯����������ָ���Աʱ��Ҫд�˺�����ʹ��ʱ�Զ�����
CMatrix::CMatrix(const CMatrix& A) //copy ructor
{
	this->m_nRows = A.m_nRows;
	this->m_nColumns = A.m_nColumns;
	if (m_nRows * m_nColumns != 0)
	{
		this->m_lpBuf = new double[m_nRows*m_nColumns];
		if (this->m_lpBuf == 0)
		{
			printf("Error,Menory allocation failure!");
			exit(0);
		}
		for (int i = 0; i < (m_nRows*m_nColumns); ++i)
			this->m_lpBuf[i] = A.m_lpBuf[i];
	}
	else
		m_lpBuf = NULL;
}
//��������������ʱ�Զ�����
CMatrix::~CMatrix()
{
	m_nRows = 0;
	m_nColumns = 0;
	if (m_lpBuf != NULL)
		delete[] this->m_lpBuf;
}
void CMatrix::resize(int m,int n)
{
	if(m_nRows!=m|| m_nColumns!=n)
		delete[] this->m_lpBuf;
	m_nRows = m;
	m_nColumns = n;
	//��Ҫ��ֹ�˷�������˴�����
	this->m_lpBuf = new double[m*n];
	if (this->m_lpBuf == 0)
	{
		printf("Error,Menory allocation failure!");
		exit(0);
	}
}
//����ĸ���
CMatrix& CMatrix::operator = (const CMatrix& A)
{
	if (this == &A)
		return *this;
	this->resize(A.m_nRows, A.m_nColumns);

	for (int i = 0; i<(m_nColumns*m_nRows); i++)
		m_lpBuf[i] = A.m_lpBuf[i];
	return *this;
}
//�����+=����
CMatrix& CMatrix::operator += (const CMatrix& A)
{
	if (!A.m_lpBuf) return *this;
	if ((this->m_nRows != A.m_nRows) || (this->m_nColumns != A.m_nColumns))
	{
		throw std::logic_error("Size mismatch in CMatrix addition");
	}
	for (int i = 0; i<A.m_nColumns*A.m_nRows; i++)
		this->m_lpBuf[i] += A.m_lpBuf[i];

	return *this;
}
//�����+����
CMatrix CMatrix::operator + (const CMatrix& A)
{
	if (!A.m_lpBuf) return *this;
	if ((this->m_nRows != A.m_nRows) || (this->m_nColumns != A.m_nColumns))
	{
		throw std::logic_error("Size mismatch in CMatrix addition");
	}
	for (int i = 0; i<A.m_nColumns*A.m_nRows; i++)
		this->m_lpBuf[i] += A.m_lpBuf[i];

	return *this;
}