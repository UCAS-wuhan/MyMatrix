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
	out << std::endl;
	for (int i = 0; i < A.m_nRows; i++) //��ӡ�����
	{
		for (int j = 0; j < A.m_nColumns; j++)
		{
			out << (A[i][j]) << "\t";
		}
		out << std::endl;
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
//������ת��
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

//��Э���������
//��������Ĭ��Ϊtrue,��m*n����Ĭ�ϼ���Э�����Ϊ1 / (m - 1)��ĳЩ������Э������㷽��Ϊ1 / m����ʱ��������Ϊflase��
//ʵ�ַ����ο���https://blog.csdn.net/mr_hhh/article/details/78490576
CMatrix CMatrix::cov(bool flag)
{
	//m_row ��������column ������
	if (m_nColumns == 0)
	{
		CMatrix m(0,0);
		return m;
	}
	double *mean = new double[m_nColumns]; //��ֵ����

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
	for (i = 0; i<m_nColumns; i++) //��һ������
	{
		for (j = i; j<m_nColumns; j++) //�ڶ�������
		{
			for (k = 0; k<m_nRows; k++) //����
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
	for (i = 0; i<m_nColumns; i++) //��ȫ��Ӧ��
	{
		for (j = 0; j<i; j++)
		{
			ret[i][j] = ret[j][i];
		}
	}
	return ret;
}

//�ݹ����
double calcDet(int n, double *&aa)
{
	if (n == 1)
		return aa[0];
	double *bb = new double[(n - 1)*(n - 1)];//����n-1�׵Ĵ�������ʽ��bb
	double sum = 0.0;
	for (int Ai = 0; Ai<n; Ai++)
	{
		for (int Bi = 0; Bi < n - 1; Bi++)//��aa���һ�и�Ԫ�صĴ�������ʽ�浽bb
		{
			int offset = Bi < Ai ? 0 : 1; //bb��С��Ai���У�ͬ�и�ֵ�����ڵĴ�������ڵļ�һ
			for (int j = 0; j<n - 1; j++)
			{
				bb[Bi*(n - 1) + j] = aa[(Bi + offset)*n + j + 1];
			}
		}
		int flag = (Ai % 2 == 0 ? 1 : -1);//��Ϊ����Ϊ0������������ż��ʱ�򣬴�������ʽΪ1.
		sum += flag* aa[Ai*n] * calcDet(n - 1, bb);//aa��һ�и�Ԫ�������������ʽ���ĺͼ�Ϊ����ʽ
	}
	delete[]bb;
	return sum;
}
//������ʽ
double CMatrix::det()
{
	if (m_nColumns == m_nRows)
		return calcDet(m_nRows, m_lpBuf);
	else
	{
		throw RCdismatch("���в�����޷�����!");
		return 0;
	}
}

//ȡ�Խ���Ԫ��
CMatrix CMatrix::diag()
{
	if (m_nRows != m_nColumns)
	{
		CMatrix m(0,0);
		throw RCdismatch("���в�����޷�����!");
		return m;
	}
	CMatrix m(m_nRows,m_nRows);
	for (int i = 0; i<m_nRows; i++)
	{
		m.m_lpBuf[i*m_nRows + i] = m_lpBuf[i*m_nRows + i];
	}
	return m;
}

//���ش�������ʽ
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
//���ذ�����
CMatrix CMatrix::adjoint()
{
	if (m_nRows != m_nColumns)
	{
		throw RCdismatch("���в�����޷�����!");
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

//�������
CMatrix CMatrix::inverse()
{
	double detOfMat = det();
	if (detOfMat == 0)
	{
		throw RCdismatch("���в�����޷�����!");
		CMatrix ret(0, 0);
		return ret;
	}
	return adjoint() / detOfMat;
}

//QR�ֽ�
void  CMatrix::QR(CMatrix &Q, CMatrix &R) const
{
	//���A����һ����ά��������ʾ���󣬺����������
	if (m_nRows != m_nColumns)
	{
		throw std::logic_error("ERROE: QR() parameter A is not a square matrix!");
		return;
	}
	const int N = m_nRows;
	double *a = new double[N];
	double *b = new double[N];

	for (int j = 0; j < N; ++j)  //(Gram-Schmidt) ����������
	{
		for (int i = 0; i < N; ++i)  //��j�е����ݴ浽a��b
			a[i] = b[i] = m_lpBuf[i * N + j];

		for (int i = 0; i<j; ++i)  //��j��֮ǰ����
		{
			R.m_lpBuf[i * N + j] = 0;  //
			for (int m = 0; m < N; ++m)
			{
				R.m_lpBuf[i * N + j] += a[m] * Q.m_lpBuf[m *N + i]; //R[i,j]ֵΪQ��i����A��j�е��ڻ�
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

		R.m_lpBuf[j*N + j] = norm; //����b[]��2�����浽R[j,j]

		for (int i = 0; i < N; ++i)
		{
			Q.m_lpBuf[i * N + j] = b[i] / norm; //Q ��ĵ�j��Ϊ��λ����b[]
		}
	}
	delete[]a;
	delete[]b;
}
//������ֵ(ʹ��QR�ֽⷽ��)
CMatrix CMatrix::eig_val(_In_opt_ int _iters)
{
	if (m_nRows*m_nColumns == 0 || m_nRows != m_nColumns)
	{
		throw RCdismatch("����Ϊ�ջ��߷Ƿ���");
		CMatrix rets(0,0);
		return rets;
	}
	//if (det() == 0)
	//{
	//  cout << "�����Ⱦ���û����QR�ֽ��������ֵ��" << endl;
	//  Matrix rets(0);
	//  return rets;
	//}
	const int N = m_nRows;
	CMatrix matcopy(*this);//���ݾ���
	CMatrix Q(N,N), R(N,N);
	/*�����������㹻��ʱ,A ���������Ǿ��������Ǿ���ĶԽ�Ԫ����A��ȫ������ֵ��*/
	for (int k = 0; k < _iters; ++k)
	{
		QR(Q, R);
		*this = R*Q;
	}
	CMatrix val = diag();
	*this = matcopy;//�ָ�ԭʼ����
	return val;
}

//����������
CMatrix CMatrix::eig_vect(_In_opt_ int _iters)
{
	if (m_nRows*m_nColumns == 0 || m_nRows != m_nColumns)
	{
		throw RCdismatch("����Ϊ�ջ��߷Ƿ���");
		CMatrix rets(0,0);
		return rets;
	}
	if (det() == 0)
	{
		throw std::logic_error("ERROE: �����Ⱦ���û����QR�ֽ��������������");
		CMatrix rets(0,0);
		return rets;
	}
	CMatrix matcopy(*this);//���ݾ���
	CMatrix eigenValue = eig_val(_iters);
	CMatrix ret(m_nRows, m_nRows);
	const int NUM = m_nColumns;
	double eValue;
	double sum, midSum, diag;
	CMatrix copym(*this);
	for (int count = 0; count < NUM; ++count)
	{
		//��������ֵΪeValue�������������ʱ��ϵ������
		*this = copym;
		eValue = eigenValue[count][count];

		for (int i = 0; i < m_nColumns; ++i)//A-lambda*I
		{
			m_lpBuf[i * m_nColumns + i] -= eValue;
		}
		//cout<<*this<<endl;
		//�� thisΪ�����͵������Ǿ���
		for (int i = 0; i < m_nRows - 1; ++i)
		{
			diag = m_lpBuf[i*m_nColumns + i];  //��ȡ�Խ�Ԫ��
			for (int j = i; j < m_nColumns; ++j)
			{
				m_lpBuf[i*m_nColumns + j] /= diag; //��i,i��Ԫ�ر�Ϊ1
			}
			for (int j = i + 1; j<m_nRows; ++j)
			{
				diag = m_lpBuf[j *  m_nColumns + i];
				for (int q = i; q < m_nColumns; ++q)//��ȥ��i+1�еĵ�i��Ԫ��
				{
					m_lpBuf[j*m_nColumns + q] -= diag*m_lpBuf[i*m_nColumns + q];
				}
			}
		}
		//cout<<*this<<endl;
		//�����������һ��Ԫ����Ϊ1
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
			ret.m_lpBuf[i * ret.m_nColumns + count] /= midSum; //ÿ�����һ��������
		}
	}
	*this = matcopy;//�ָ�ԭʼ����
	return ret;
}

//�Ծ���Ԫ������trueΪ��С����
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

//ʹ��Cramer��������������Է�����Ax=b�����ؾ���x
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
		throw RCdismatch("Cramer����ֻ�ܼ���ϵ������Ϊ���ȵľ���");
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
			Dx[i*m_nRows + j] = obj.m_lpBuf[i]; //obj��ֵ����j��
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

//������1����
double CMatrix::norm1()
{
	double sum = 0;
	for (int i = 0; i < m_nRows*m_nColumns; ++i)
	{
		sum += abs(m_lpBuf[i]);
	}
	return sum;
}
//||matrix||_2  ��A�����2����
double CMatrix::norm2()
{
	double norm = 0;
	for (int i = 0; i < m_nRows*m_nColumns; ++i)
	{
		norm += m_lpBuf[i] * m_lpBuf[i];
	}
	return (double)sqrt(norm);
}
//������ֵ
double CMatrix::mean()
{
	double sum = 0;
	for (int i = 0; i < m_nRows*m_nColumns; ++i)
	{
		sum += (m_lpBuf[i]);
	}
	return sum / (m_nRows*m_nColumns);
}

//���Ļ������ֵ����
void CMatrix::zeromean(_In_opt_  bool flag)
{
	if (flag == true) //�����о�ֵ
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
	else //�����о�ֵ
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
//��׼������һ����
void CMatrix::normalize(_In_opt_  bool flag)
{
	if (flag == true) //�����о�ֵ
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
		//�����׼��
		for (int j = 0; j < m_nColumns; j++)
		{
			mean[j] = 0;
			for (int i = 0; i < m_nRows; i++)
			{
				mean[j] += pow(m_lpBuf[i*m_nColumns + j], 2);//��ƽ����
			}
			mean[j] = sqrt(mean[j] / m_nRows); // ����
		}
		for (int j = 0; j < m_nColumns; j++)
		{
			for (int i = 0; i < m_nRows; i++)
			{
				m_lpBuf[i*m_nColumns + j] /= mean[j];//��ƽ����
			}
		}
		delete[]mean;
	}
	else //�����о�ֵ
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
		//�����׼��
		for (int i = 0; i< m_nRows; i++)
		{
			mean[i] = 0.0;
			for (int j = 0; j < m_nColumns; j++)
			{
				mean[i] += pow(m_lpBuf[i*m_nColumns + j], 2);//��ƽ����
			}
			mean[i] = sqrt(mean[i] / m_nColumns); // ����
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
//ÿ��Ԫ��x����
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
//�Խ���
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
