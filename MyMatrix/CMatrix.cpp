#include "CMatrix.h"


//���캯������ʼ��ʱ�Զ�����
CMatrix::CMatrix(int rows, int cols)
{
	m_nRows = rows;
	m_nColumns = cols;
	//��Ҫ��ֹ�˷�������˴�����
	this->m_lpBuf = new double[rows*cols];
	if (this->m_lpBuf == 0)
	{
		printf("Error,Menory allocation failure!");
		exit(0);
	}
	memset(m_lpBuf, 0, rows*cols * sizeof(double));
}
//����һ��ֵȫ��Ϊvalue�ľ���
CMatrix::CMatrix(int rows, int cols, double value)
{
	m_nRows = rows;
	m_nColumns = cols;
	//��Ҫ��ֹ�˷�������˴�����
	this->m_lpBuf = new double[rows*cols];
	if (this->m_lpBuf == 0)
	{
		printf("Error,Menory allocation failure!");
		exit(0);
	}
	for (int i = 0; i < rows*cols; i++)
	{
		this->m_lpBuf[i] = value;
	}
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


//Բ����()����Ϊ��֤���±���������±�i��j(0<=i<row,0<=j<col)
//�������ĵ�i�е�j��Ԫ�أ�����һά�����е�λ����(i-1)*col+(j-1)
double& CMatrix::operator()(int i, int j)const
{
	if (i <= 0 || i > this->m_nRows || j <= 0 || j > this->m_nColumns)
		throw Bounds_error("�±�Խ��!");//�쳣
	return this->m_lpBuf[(i - 1)*this->m_nColumns + (j - 1)];
}

//��������С
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
//���������
std::istream& operator >> (std::istream &is, CMatrix &A)
{
	for (int i = 0; i<A.m_nColumns*A.m_nRows; i++)
	{
		is >> A.m_lpBuf[i];
	}
	return is;
}
//����Ĵ�ӡ����дstd��<<��������
std::ostream& operator << (std::ostream &out, CMatrix &A)
{
	for (int i = 0; i < A.m_nColumns*A.m_nRows; i++) 
	{
		out << A.m_lpBuf[i];
		if(i%A.m_nColumns == 0)
		{
			out << std::endl;
		}
		else out << "\t";
	}
	return out;
}

//����ĸ���
CMatrix& CMatrix::operator = (const CMatrix& A)
{
	if (this == &A)
		return *this;
	if ((this->m_nRows != A.m_nRows) || (this->m_nColumns != A.m_nColumns))
	{
		throw std::logic_error("Size mismatch in CMatrix");
	}
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
		throw std::logic_error("Size mismatch in CMatrix");
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
//�����-=����
CMatrix& CMatrix::operator -= (const CMatrix& A)
{
	if (!A.m_lpBuf) return *this;
	if ((this->m_nRows != A.m_nRows) || (this->m_nColumns != A.m_nColumns))
	{
		throw std::logic_error("Size mismatch in CMatrix");
	}
	for (int i = 0; i<A.m_nColumns*A.m_nRows; i++)
		this->m_lpBuf[i] -= A.m_lpBuf[i];

	return *this;
}
//�����-����
CMatrix CMatrix::operator - (const CMatrix& A)
{
	if (!A.m_lpBuf) return *this;
	if ((this->m_nRows != A.m_nRows) || (this->m_nColumns != A.m_nColumns))
	{
		throw std::logic_error("Size mismatch in CMatrix");
	}
	for (int i = 0; i<A.m_nColumns*A.m_nRows; i++)
		this->m_lpBuf[i] -= A.m_lpBuf[i];

	return *this;
}
//���������
CMatrix operator*(const CMatrix& lm, const CMatrix& rm) 
{
	if (lm.m_nColumns*lm.m_nRows == 0 || rm.m_nColumns*rm.m_nRows == 0 || lm.m_nColumns != rm.m_nRows)
	{
		CMatrix temp(0, 0);
		temp.m_lpBuf = NULL;
		throw std::logic_error("Size mismatch in CMatrix");
		return temp; //���ݲ��Ϸ�ʱ�򣬷��ؿվ���
	}
	CMatrix ret(lm.m_nRows, rm.m_nColumns);
	for (int i = 0; i<lm.m_nRows; i++)
	{
		for (int j = 0; j< rm.m_nColumns; j++)
		{
			for (int k = 0; k< lm.m_nColumns; k++)//lm.m_col == rm.m_row
			{
				ret.m_lpBuf[i*rm.m_nColumns + j] += lm.m_lpBuf[i*lm.m_nColumns + k] * rm.m_lpBuf[k*rm.m_nColumns + j];
			}
		}
	}
	return ret;
}
//������ ����
CMatrix operator*(double val, const CMatrix& rm)
{
	CMatrix ret(rm.m_nRows, rm.m_nColumns);
	for (int i = 0; i<ret.m_nRows*ret.m_nColumns; i++)
	{
		ret.m_lpBuf[i] = val * rm.m_lpBuf[i];
	}
	return ret;
}
//����� ����
CMatrix operator*(const CMatrix&lm, double val)  
{
	CMatrix ret(lm.m_nRows, lm.m_nColumns);
	for (int i = 0; i<ret.m_nRows*ret.m_nColumns; i++)
	{
		ret.m_lpBuf[i] = val * lm.m_lpBuf[i];
	}
	return ret;
}
//������� ����
CMatrix operator/(const CMatrix&lm, double val)  
{
	CMatrix ret(lm.m_nRows, lm.m_nColumns);
	if (val == 0)
	{
		throw std::logic_error("Run time error, divisor 0");
		return ret;
	}
	for (int i = 0; i<ret.m_nRows*ret.m_nColumns; i++)
	{
		ret.m_lpBuf[i] = lm.m_lpBuf[i] / val;
	}
	return ret;
}
// ��ӦԪ�����
CMatrix CMatrix::multi(const CMatrix&rm)
{
	if (m_nColumns != rm.m_nColumns || m_nRows != rm.m_nRows)
	{
		CMatrix temp(0, 0);
		temp.m_lpBuf = NULL;
		throw std::logic_error("Size mismatch in CMatrix");
		return temp; //���ݲ��Ϸ�ʱ�򣬷��ؿվ���
	}
	CMatrix ret(m_nRows, m_nColumns);
	for (int i = 0; i<ret.m_nRows*ret.m_nColumns; i++)
	{
		ret.m_lpBuf[i] = m_lpBuf[i] * rm.m_lpBuf[i];
	}
	return ret;

}
//��ȡ����ĵ�index������
CMatrix CMatrix::getrow(int index)
{
	CMatrix ret(1, m_nColumns); //һ�еķ���ֵ

	for (int i = 0; i< m_nColumns; i++)
	{
		ret[0][i] = m_lpBuf[(index-1)*m_nColumns + i];

	}
	return ret;
}
//��ȡ����ĵ�index������
CMatrix CMatrix::getcol(int index)
{
	CMatrix ret(m_nRows, 1); //һ�еķ���ֵ

	for (int i = 0; i< m_nRows; i++)
	{
		ret[i][0] = m_lpBuf[i *m_nColumns + (index - 1)];

	}
	return ret;
}
