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
	out << std::endl;
	for (int i = 0; i < A.m_nRows; i++) //打印逆矩阵
	{
		for (int j = 0; j < A.m_nColumns; j++)
		{
			out << (A[i][j]) << "\t";
		}
		out << std::endl;
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
//求矩阵的转置
CMatrix CMatrix::T()const
{
	CMatrix tem(m_nColumns, m_nRows);
	for (int i = 0; i<this->m_nRows; i++)
	{
		for (int j = 0; j<this->m_nColumns; j++)
		{
			tem[j][i] = m_lpBuf[i*m_nColumns + j];// (*this)[i][j]
		}
	}
	return tem;
}

//求协方差矩阵阵
//函数参数默认为true,对m*n矩阵，默认计算协方差方法为1 / (m - 1)，某些场合下协方差计算方法为1 / m，此时将参数设为flase。
//实现方法参考：https://blog.csdn.net/mr_hhh/article/details/78490576
CMatrix CMatrix::cov(bool flag)
{
	//m_row 样本数，column 变量数
	if (m_nColumns == 0)
	{
		CMatrix m(0,0);
		return m;
	}
	double *mean = new double[m_nColumns]; //均值向量

	for (int j = 0; j<m_nColumns; j++) //init
	{
		mean[j] = 0.0;
	}
	CMatrix ret(m_nColumns, m_nColumns);
	for (int j = 0; j<m_nColumns; j++) //mean
	{
		for (int i = 0; i<m_nRows; i++)
		{
			mean[j] += m_lpBuf[i*m_nColumns + j];
		}
		mean[j] /= m_nRows;
	}
	int i, k, j;
	for (i = 0; i<m_nColumns; i++) //第一个变量
	{
		for (j = i; j<m_nColumns; j++) //第二个变量
		{
			for (k = 0; k<m_nRows; k++) //计算
			{
				ret[i][j] += (m_lpBuf[k*m_nColumns + i] - mean[i])*(m_lpBuf[k*m_nColumns + j] - mean[j]);

			}
			if (flag == true)
			{
				ret[i][j] /= (m_nRows - 1);
			}
			else
			{
				ret[i][j] /= (m_nRows);
			}
		}
	}
	for (i = 0; i<m_nColumns; i++) //补全对应面
	{
		for (j = 0; j<i; j++)
		{
			ret[i][j] = ret[j][i];
		}
	}
	return ret;
}

//递归调用
double calcDet(int n, double *&aa)
{
	if (n == 1)
		return aa[0];
	double *bb = new double[(n - 1)*(n - 1)];//创建n-1阶的代数余子式阵bb
	double sum = 0.0;
	for (int Ai = 0; Ai<n; Ai++)
	{
		for (int Bi = 0; Bi < n - 1; Bi++)//把aa阵第一列各元素的代数余子式存到bb
		{
			int offset = Bi < Ai ? 0 : 1; //bb中小于Ai的行，同行赋值，等于的错过，大于的加一
			for (int j = 0; j<n - 1; j++)
			{
				bb[Bi*(n - 1) + j] = aa[(Bi + offset)*n + j + 1];
			}
		}
		int flag = (Ai % 2 == 0 ? 1 : -1);//因为列数为0，所以行数是偶数时候，代数余子式为1.
		sum += flag* aa[Ai*n] * calcDet(n - 1, bb);//aa第一列各元素与其代数余子式积的和即为行列式
	}
	delete[]bb;
	return sum;
}
//求行列式
double CMatrix::det()
{
	if (m_nColumns == m_nRows)
		return calcDet(m_nRows, m_lpBuf);
	else
	{
		throw RCdismatch("行列不相等无法计算!");
		return 0;
	}
}

//取对角线元素
CMatrix CMatrix::diag()
{
	if (m_nRows != m_nColumns)
	{
		CMatrix m(0,0);
		throw RCdismatch("行列不相等无法计算!");
		return m;
	}
	CMatrix m(m_nRows,m_nRows);
	for (int i = 0; i<m_nRows; i++)
	{
		m.m_lpBuf[i*m_nRows + i] = m_lpBuf[i*m_nRows + i];
	}
	return m;
}

//返回代数余子式
double CalcAlgebraicCofactor(CMatrix& srcMat, int ai, int aj)
{
	int temMatLen = srcMat.row() - 1;
	CMatrix temMat(temMatLen, temMatLen);
	for (int bi = 0; bi < temMatLen; bi++)
	{
		for (int bj = 0; bj < temMatLen; bj++)
		{
			int rowOffset = bi < ai ? 0 : 1;
			int colOffset = bj < aj ? 0 : 1;
			temMat[bi][bj] = srcMat[bi + rowOffset][bj + colOffset];
		}
	}
	int flag = (ai + aj) % 2 == 0 ? 1 : -1;
	return flag * temMat.det();
}
//返回伴随阵
CMatrix CMatrix::adjoint()
{
	if (m_nRows != m_nColumns)
	{
		throw RCdismatch("行列不相等无法计算!");
		return CMatrix(0,0);
	}

	CMatrix adjointMat(m_nRows, m_nRows);
	for (int ai = 0; ai < m_nRows; ai++)
	{
		for (int aj = 0; aj < m_nRows; aj++)
		{
			adjointMat.m_lpBuf[aj*m_nRows + ai] = CalcAlgebraicCofactor(*this, ai, aj);
		}
	}
	return adjointMat;
}

//求逆矩阵
CMatrix CMatrix::inverse()
{
	double detOfMat = det();
	if (detOfMat == 0)
	{
		throw RCdismatch("行列不相等无法计算!");
		CMatrix ret(0, 0);
		return ret;
	}
	return adjoint() / detOfMat;
}

//QR分解
void  CMatrix::QR(CMatrix &Q, CMatrix &R) const
{
	//如果A不是一个二维方阵，则提示错误，函数计算结束
	if (m_nRows != m_nColumns)
	{
		throw std::logic_error("ERROE: QR() parameter A is not a square matrix!");
		return;
	}
	const int N = m_nRows;
	double *a = new double[N];
	double *b = new double[N];

	for (int j = 0; j < N; ++j)  //(Gram-Schmidt) 正交化方法
	{
		for (int i = 0; i < N; ++i)  //第j列的数据存到a，b
			a[i] = b[i] = m_lpBuf[i * N + j];

		for (int i = 0; i<j; ++i)  //第j列之前的列
		{
			R.m_lpBuf[i * N + j] = 0;  //
			for (int m = 0; m < N; ++m)
			{
				R.m_lpBuf[i * N + j] += a[m] * Q.m_lpBuf[m *N + i]; //R[i,j]值为Q第i列与A的j列的内积
			}
			for (int m = 0; m < N; ++m)
			{
				b[m] -= R.m_lpBuf[i * N + j] * Q.m_lpBuf[m * N + i]; //
			}
		}

		double norm = 0;
		for (int i = 0; i < N; ++i)
		{
			norm += b[i] * b[i];
		}
		norm = (double)sqrt(norm);

		R.m_lpBuf[j*N + j] = norm; //向量b[]的2范数存到R[j,j]

		for (int i = 0; i < N; ++i)
		{
			Q.m_lpBuf[i * N + j] = b[i] / norm; //Q 阵的第j列为单位化的b[]
		}
	}
	delete[]a;
	delete[]b;
}
//求特征值(使用QR分解方法)
CMatrix CMatrix::eig_val(_In_opt_ int _iters)
{
	if (m_nRows*m_nColumns == 0 || m_nRows != m_nColumns)
	{
		throw RCdismatch("矩阵为空或者非方阵");
		CMatrix rets(0,0);
		return rets;
	}
	//if (det() == 0)
	//{
	//  cout << "非满秩矩阵没法用QR分解计算特征值！" << endl;
	//  Matrix rets(0);
	//  return rets;
	//}
	const int N = m_nRows;
	CMatrix matcopy(*this);//备份矩阵
	CMatrix Q(N,N), R(N,N);
	/*当迭代次数足够多时,A 趋于上三角矩阵，上三角矩阵的对角元就是A的全部特征值。*/
	for (int k = 0; k < _iters; ++k)
	{
		QR(Q, R);
		*this = R*Q;
	}
	CMatrix val = diag();
	*this = matcopy;//恢复原始矩阵；
	return val;
}

//求特征向量
CMatrix CMatrix::eig_vect(_In_opt_ int _iters)
{
	if (m_nRows*m_nColumns == 0 || m_nRows != m_nColumns)
	{
		throw RCdismatch("矩阵为空或者非方阵");
		CMatrix rets(0,0);
		return rets;
	}
	if (det() == 0)
	{
		throw std::logic_error("ERROE: 非满秩矩阵没法用QR分解计算特征向量！");
		CMatrix rets(0,0);
		return rets;
	}
	CMatrix matcopy(*this);//备份矩阵
	CMatrix eigenValue = eig_val(_iters);
	CMatrix ret(m_nRows, m_nRows);
	const int NUM = m_nColumns;
	double eValue;
	double sum, midSum, diag;
	CMatrix copym(*this);
	for (int count = 0; count < NUM; ++count)
	{
		//计算特征值为eValue，求解特征向量时的系数矩阵
		*this = copym;
		eValue = eigenValue[count][count];

		for (int i = 0; i < m_nColumns; ++i)//A-lambda*I
		{
			m_lpBuf[i * m_nColumns + i] -= eValue;
		}
		//cout<<*this<<endl;
		//将 this为阶梯型的上三角矩阵
		for (int i = 0; i < m_nRows - 1; ++i)
		{
			diag = m_lpBuf[i*m_nColumns + i];  //提取对角元素
			for (int j = i; j < m_nColumns; ++j)
			{
				m_lpBuf[i*m_nColumns + j] /= diag; //【i,i】元素变为1
			}
			for (int j = i + 1; j<m_nRows; ++j)
			{
				diag = m_lpBuf[j *  m_nColumns + i];
				for (int q = i; q < m_nColumns; ++q)//消去第i+1行的第i个元素
				{
					m_lpBuf[j*m_nColumns + q] -= diag*m_lpBuf[i*m_nColumns + q];
				}
			}
		}
		//cout<<*this<<endl;
		//特征向量最后一行元素置为1
		midSum = ret.m_lpBuf[(ret.m_nRows - 1) * ret.m_nColumns + count] = 1;
		for (int m = m_nRows - 2; m >= 0; --m)
		{
			sum = 0;
			for (int j = m + 1; j < m_nColumns; ++j)
			{
				sum += m_lpBuf[m *  m_nColumns + j] * ret.m_lpBuf[j * ret.m_nColumns + count];
			}
			sum = -sum / m_lpBuf[m *  m_nColumns + m];
			midSum += sum * sum;
			ret.m_lpBuf[m * ret.m_nColumns + count] = sum;
		}
		midSum = sqrt(midSum);
		for (int i = 0; i < ret.m_nRows; ++i)
		{
			ret.m_lpBuf[i * ret.m_nColumns + count] /= midSum; //每次求出一个列向量
		}
	}
	*this = matcopy;//恢复原始矩阵；
	return ret;
}

//对矩阵元素排序，true为从小到大
void CMatrix::sort(bool flag)
{
	double tem;
	for (int i = 0; i<m_nRows*m_nColumns; i++)
	{
		for (int j = i + 1; j<m_nRows*m_nColumns; j++)
		{
			if (flag == true)
			{
				if (m_lpBuf[i]>m_lpBuf[j])
				{
					tem = m_lpBuf[i];
					m_lpBuf[i] = m_lpBuf[j];
					m_lpBuf[j] = tem;
				}
			}
			else
			{
				if (m_lpBuf[i]<m_lpBuf[j])
				{
					tem = m_lpBuf[i];
					m_lpBuf[i] = m_lpBuf[j];
					m_lpBuf[j] = tem;
				}
			}

		}
	}
}

//使用Cramer法则求解非齐次线性方程组Ax=b，返回矩阵x
CMatrix CMatrix::solveAb(CMatrix &obj)
{
	CMatrix ret(m_nRows, 1);
	if (m_nRows*m_nColumns == 0 || obj.m_nRows*obj.m_nColumns == 0)
	{
		throw std::logic_error("solveAb(Matrix &obj):this or obj is null");
		return ret;
	}
	if (m_nRows != obj.m_nRows*obj.m_nColumns)
	{
		throw RCdismatch("solveAb(Matrix &obj):the row of two matrix is not equal!");
		return ret;
	}

	double *Dx = new double[m_nRows*m_nRows];
	for (int i = 0; i<m_nRows; i++)
	{
		for (int j = 0; j<m_nRows; j++)
		{
			Dx[i*m_nRows + j] = m_lpBuf[i*m_nRows + j];
		}
	}
	double D = calcDet(m_nRows, Dx);
	if (D == 0)
	{
		throw RCdismatch("Cramer法则只能计算系数矩阵为满秩的矩阵");
		return  ret;
	}

	for (int j = 0; j<m_nRows; j++)
	{
		for (int i = 0; i<m_nRows; i++)
		{
			for (int j = 0; j<m_nRows; j++)
			{
				Dx[i*m_nRows + j] = m_lpBuf[i*m_nRows + j];
			}
		}
		for (int i = 0; i<m_nRows; i++)
		{
			Dx[i*m_nRows + j] = obj.m_lpBuf[i]; //obj赋值给第j列
		}

		//for( int i=0;i<m_row;i++) //print
		//{
		//    for(int j=0; j<m_row;j++)
		//    {
		//        cout<< Dx[i*m_row+j]<<"\t";
		//    }
		//    cout<<endl;
		//}
		ret[j][0] = calcDet(m_nRows, Dx) / D;

	}

	delete[]Dx;
	return ret;
}

//求矩阵的1范数
double CMatrix::norm1()
{
	double sum = 0;
	for (int i = 0; i < m_nRows*m_nColumns; ++i)
	{
		sum += abs(m_lpBuf[i]);
	}
	return sum;
}
//||matrix||_2  求A矩阵的2范数
double CMatrix::norm2()
{
	double norm = 0;
	for (int i = 0; i < m_nRows*m_nColumns; ++i)
	{
		norm += m_lpBuf[i] * m_lpBuf[i];
	}
	return (double)sqrt(norm);
}
//求矩阵均值
double CMatrix::mean()
{
	double sum = 0;
	for (int i = 0; i < m_nRows*m_nColumns; ++i)
	{
		sum += (m_lpBuf[i]);
	}
	return sum / (m_nRows*m_nColumns);
}

//中心化（零均值化）
void CMatrix::zeromean(_In_opt_  bool flag)
{
	if (flag == true) //计算列均值
	{
		double *mean = new double[m_nColumns];
		for (int j = 0; j < m_nColumns; j++)
		{
			mean[j] = 0.0;
			for (int i = 0; i < m_nRows; i++)
			{
				mean[j] += m_lpBuf[i*m_nColumns + j];
			}
			mean[j] /= m_nRows;
		}
		for (int j = 0; j < m_nColumns; j++)
		{

			for (int i = 0; i < m_nRows; i++)
			{
				m_lpBuf[i*m_nColumns + j] -= mean[j];
			}
		}
		delete[]mean;
	}
	else //计算行均值
	{
		double *mean = new double[m_nRows];
		for (int i = 0; i< m_nRows; i++)
		{
			mean[i] = 0.0;
			for (int j = 0; j < m_nColumns; j++)
			{
				mean[i] += m_lpBuf[i*m_nColumns + j];
			}
			mean[i] /= m_nColumns;
		}
		for (int i = 0; i < m_nRows; i++)
		{
			for (int j = 0; j < m_nColumns; j++)
			{
				m_lpBuf[i*m_nColumns + j] -= mean[i];
			}
		}
		delete[]mean;
	}
}
//标准化（归一化）
void CMatrix::normalize(_In_opt_  bool flag)
{
	if (flag == true) //计算列均值
	{
		double *mean = new double[m_nColumns];

		for (int j = 0; j < m_nColumns; j++)
		{
			mean[j] = 0.0;
			for (int i = 0; i < m_nRows; i++)
			{
				mean[j] += m_lpBuf[i*m_nColumns + j];
			}
			mean[j] /= m_nRows;
		}
		for (int j = 0; j < m_nColumns; j++)
		{

			for (int i = 0; i < m_nRows; i++)
			{
				m_lpBuf[i*m_nColumns + j] -= mean[j];
			}
		}
		//计算标准差
		for (int j = 0; j < m_nColumns; j++)
		{
			mean[j] = 0;
			for (int i = 0; i < m_nRows; i++)
			{
				mean[j] += pow(m_lpBuf[i*m_nColumns + j], 2);//列平方和
			}
			mean[j] = sqrt(mean[j] / m_nRows); // 开方
		}
		for (int j = 0; j < m_nColumns; j++)
		{
			for (int i = 0; i < m_nRows; i++)
			{
				m_lpBuf[i*m_nColumns + j] /= mean[j];//列平方和
			}
		}
		delete[]mean;
	}
	else //计算行均值
	{
		double *mean = new double[m_nRows];
		for (int i = 0; i< m_nRows; i++)
		{
			mean[i] = 0.0;
			for (int j = 0; j < m_nColumns; j++)
			{
				mean[i] += m_lpBuf[i*m_nColumns + j];
			}
			mean[i] /= m_nColumns;
		}
		for (int i = 0; i < m_nRows; i++)
		{
			for (int j = 0; j < m_nColumns; j++)
			{
				m_lpBuf[i*m_nColumns + j] -= mean[i];
			}
		}
		//计算标准差
		for (int i = 0; i< m_nRows; i++)
		{
			mean[i] = 0.0;
			for (int j = 0; j < m_nColumns; j++)
			{
				mean[i] += pow(m_lpBuf[i*m_nColumns + j], 2);//列平方和
			}
			mean[i] = sqrt(mean[i] / m_nColumns); // 开方
		}
		for (int i = 0; i < m_nRows; i++)
		{
			for (int j = 0; j < m_nColumns; j++)
			{
				m_lpBuf[i*m_nColumns + j] /= mean[i];
			}
		}
		delete[]mean;
	}
}
//每个元素x次幂
CMatrix CMatrix::exponent(double x)
{
	CMatrix ret(m_nRows, m_nColumns);
	for (int i = 0; i< m_nRows; i++)
	{
		for (int j = 0; j < m_nColumns; j++)
		{
			ret[i][j] = pow(m_lpBuf[i*m_nColumns + j], x);
		}
	}
	return ret;
}
//对角阵
CMatrix CMatrix::eye()
{

	for (int i = 0; i< m_nRows; i++)
	{
		for (int j = 0; j < m_nColumns; j++)
		{
			if (i == j)
			{
				m_lpBuf[i*m_nColumns + j] = 1.0;
			}
		}
	}
	return *this;
}
