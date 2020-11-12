#pragma once
#include <string.h>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
//��д��ͨ�ļ��㺯��
//�ٰѼ��㺯�����ɺ���ģ��
//�Ѽ��㺯���ŵ����У���һ����Ĺ��졢�����������
//��������ģ�壨��һ����ģ��ͺ���ģ�������÷���

//�����࣬��������Ķ��弰������
class CMatrix
{
private:
	int  m_nRows; // ���������
	int  m_nColumns; // ���������
	double* m_lpBuf; // ��̬���������������Ŀռ�
public:
	CMatrix(int, int);
	CMatrix(int, int, double);	//����ʱ��ʼ��Ϊָ��ֵ
	CMatrix(const CMatrix& A);//������캯��
	virtual ~CMatrix();//����Ϊ����ʱӦȥ���ؼ���

	double& CMatrix::operator()(int i, int j)const;
	double* operator[](int i) { return m_lpBuf + (i-1)*m_nColumns; }//ע��this�����ţ� (*this)[i][j]
	void resize(int, int);//���·���ռ�
	friend std::istream &operator >> (std::istream&,CMatrix&);
	friend std::ostream &operator << (std::ostream&, CMatrix&);     // �������Ļ

	CMatrix& operator = (const CMatrix&);//����ĸ���
	CMatrix& operator+=(const CMatrix& );//�����+=����
	CMatrix operator + (const CMatrix& A);//+
	CMatrix& operator-=(const CMatrix&);//-=
	CMatrix operator - (const CMatrix& A);//-
	friend CMatrix  operator*(const CMatrix&, const CMatrix&);  //���������
	friend CMatrix  operator*(double, const CMatrix&);  //�����˾���
	friend CMatrix  operator*(const CMatrix&, double);  //����˵���
	friend CMatrix  operator/(const CMatrix&, double);  //���� ���Ե���

	CMatrix multi(const CMatrix&); // ��ӦԪ�����
	int row()const { return m_nRows; }	//�����к�
	int col()const { return m_nColumns; } //�����к�
	CMatrix getrow(int index); // ���ص�index ��,������1����
	CMatrix getcol(int index); // ���ص�index ��,������1����

	
	//CMatrix cov(_In_opt_ bool flag = true);   //Э������ ������������
	//double det();   //����ʽ
	//CMatrix solveAb(CMatrix &obj);  // b������������������
	//CMatrix diag();  //���ضԽ���Ԫ��
	//				//Matrix asigndiag();  //�Խ���Ԫ��
	//CMatrix T()const;   //ת��
	//void sort(bool);//trueΪ��С����
	//CMatrix adjoint();
	//CMatrix inverse();
	//void QR(_Out_ CMatrix&, _Out_ CMatrix&)const;
	//CMatrix eig_val(_In_opt_ Index_T _iters = 1000);
	//CMatrix eig_vect(_In_opt_ Index_T _iters = 1000);

	//double norm1();//1����
	//double norm2();//2����
	//double mean();// �����ֵ
	//double*operator[](int i) { return m_lpBuf + i*m_nColumns; }//ע��this�����ţ� (*this)[i][j]
	//void zeromean(_In_opt_  bool flag = true);//Ĭ�ϲ���Ϊtrue������
	//void normalize(_In_opt_  bool flag = true);//Ĭ�ϲ���Ϊtrue������
	//CMatrix exponent(double x);//ÿ��Ԫ��x����
	//CMatrix  eye();//�Խ���
	//void  maxlimit(double max, double set = 0);//�Խ���

};

class Bounds_error                                          //�±�Խ���쳣
{
public:
	Bounds_error(char* s) :str(s) {}
	void display()
	{
		std::cout << str << std::endl;
	}
private:
	std::string str;
};

class RCdismatch                                           //���в�ƥ���쳣
{
public:
	RCdismatch(char* s) :str(s) {}
	void display()
	{
		std::cout << str << std::endl;
	}
private:
	std::string str;
};
