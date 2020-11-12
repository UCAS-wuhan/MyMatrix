#include "CMatrix.h"


//构造函数，初始化时自动调用
CMatrix::CMatrix(int rows, int cols)
{
	m_nRows = rows;
	m_nColumns = cols;
	//需要防止乘法溢出，此处忽略
	this->m_lpBuf = new double[rows*cols];
	if (this->m_lpBuf == 0)
	{
		printf("Error,Menory allocation failure!");
		exit(0);
	}
	memset(m_lpBuf, 0, rows*cols * sizeof(double));
}
//声明一个值全部为value的矩阵
CMatrix::CMatrix(int rows, int cols, double value)
{
	m_nRows = rows;
	m_nColumns = cols;
	//需要防止乘法溢出，此处忽略
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
//拷贝构造函数，类中有指针成员时需要写此函数，使用时自动调用
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
//析构函数，销毁时自动调用
CMatrix::~CMatrix()
{
	m_nRows = 0;
	m_nColumns = 0;
	if (m_lpBuf != NULL)
		delete[] this->m_lpBuf;
}


//圆括号()重载为举证的下标操作，对下标i和j(0<=i<row,0<=j<col)
//输出矩阵的第i行第j列元素，它在一维数组中的位置是(i-1)*col+(j-1)
double& CMatrix::operator()(int i, int j)const
{
	if (i <= 0 || i > this->m_nRows || j <= 0 || j > this->m_nColumns)
		throw Bounds_error("下标越界!");//异常
	return this->m_lpBuf[(i - 1)*this->m_nColumns + (j - 1)];
}

//重设矩阵大小
void CMatrix::resize(int m,int n)
{
	if(m_nRows!=m|| m_nColumns!=n)
		delete[] this->m_lpBuf;
	m_nRows = m;
	m_nColumns = n;
	//需要防止乘法溢出，此处忽略
	this->m_lpBuf = new double[m*n];
	if (this->m_lpBuf == 0)
	{
		printf("Error,Menory allocation failure!");
		exit(0);
	}
}
//矩阵的输入
std::istream& operator >> (std::istream &is, CMatrix &A)
{
	for (int i = 0; i<A.m_nColumns*A.m_nRows; i++)
	{
		is >> A.m_lpBuf[i];
	}
	return is;
}
//矩阵的打印（重写std的<<操作符）
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

//矩阵的复制
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
//矩阵的+=操作
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
//矩阵的+操作
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
//矩阵的-=操作
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
//矩阵的-操作
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
//两矩阵相乘
CMatrix operator*(const CMatrix& lm, const CMatrix& rm) 
{
	if (lm.m_nColumns*lm.m_nRows == 0 || rm.m_nColumns*rm.m_nRows == 0 || lm.m_nColumns != rm.m_nRows)
	{
		CMatrix temp(0, 0);
		temp.m_lpBuf = NULL;
		throw std::logic_error("Size mismatch in CMatrix");
		return temp; //数据不合法时候，返回空矩阵
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
//单数乘 矩阵
CMatrix operator*(double val, const CMatrix& rm)
{
	CMatrix ret(rm.m_nRows, rm.m_nColumns);
	for (int i = 0; i<ret.m_nRows*ret.m_nColumns; i++)
	{
		ret.m_lpBuf[i] = val * rm.m_lpBuf[i];
	}
	return ret;
}
//矩阵乘 单数
CMatrix operator*(const CMatrix&lm, double val)  
{
	CMatrix ret(lm.m_nRows, lm.m_nColumns);
	for (int i = 0; i<ret.m_nRows*ret.m_nColumns; i++)
	{
		ret.m_lpBuf[i] = val * lm.m_lpBuf[i];
	}
	return ret;
}
//矩阵除以 单数
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
// 对应元素相乘
CMatrix CMatrix::multi(const CMatrix&rm)
{
	if (m_nColumns != rm.m_nColumns || m_nRows != rm.m_nRows)
	{
		CMatrix temp(0, 0);
		temp.m_lpBuf = NULL;
		throw std::logic_error("Size mismatch in CMatrix");
		return temp; //数据不合法时候，返回空矩阵
	}
	CMatrix ret(m_nRows, m_nColumns);
	for (int i = 0; i<ret.m_nRows*ret.m_nColumns; i++)
	{
		ret.m_lpBuf[i] = m_lpBuf[i] * rm.m_lpBuf[i];
	}
	return ret;

}
//获取矩阵的第index行数据
CMatrix CMatrix::getrow(int index)
{
	CMatrix ret(1, m_nColumns); //一行的返回值

	for (int i = 0; i< m_nColumns; i++)
	{
		ret[0][i] = m_lpBuf[(index-1)*m_nColumns + i];

	}
	return ret;
}
//获取矩阵的第index列数据
CMatrix CMatrix::getcol(int index)
{
	CMatrix ret(m_nRows, 1); //一列的返回值

	for (int i = 0; i< m_nRows; i++)
	{
		ret[i][0] = m_lpBuf[i *m_nColumns + (index - 1)];

	}
	return ret;
}
