#pragma once
#include <string.h>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
//��д��ͨ�ļ��㺯��
//�ٰѼ��㺯�����ɺ���ģ��
//�Ѽ��㺯���ŵ����У���һ����Ĺ��졢�����������
//��������ģ�壨��һ����ģ��ͺ���ģ�������÷���

//�����࣬��������Ķ��弰������
class CMatrix
{
public:
	int  m_nRows; // ���������
	int  m_nColumns; // ���������
	double* m_lpBuf; // ��̬���������������Ŀռ�
public:
	CMatrix(int, int);
	CMatrix::CMatrix(const CMatrix& A);
	virtual ~CMatrix();//����Ϊ����ʱӦȥ���ؼ���

	void resize(int, int);//���·���ռ�
	CMatrix& operator = (const CMatrix& A);//����ĸ���
	CMatrix& operator+=(const CMatrix&);//�����+=����
	CMatrix operator + (const CMatrix& A);//�����+����
};

