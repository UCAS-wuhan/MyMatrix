#pragma once
#include <string.h>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
//先写普通的计算函数
//再把计算函数做成函数模板
//把计算函数放到类中（查一下类的构造、析构、深拷贝）
//把类做成模板（查一下类模板和函数模板区别用法）

//矩阵类，声明矩阵的定义及处理函数
class CMatrix
{
public:
	int  m_nRows; // 矩阵的行数
	int  m_nColumns; // 矩阵的列数
	double* m_lpBuf; // 动态分配用来存放数组的空间
public:
	CMatrix(int, int);
	CMatrix::CMatrix(const CMatrix& A);
	virtual ~CMatrix();//不作为基类时应去掉关键字

	void resize(int, int);//重新分配空间
	CMatrix& operator = (const CMatrix& A);//矩阵的复制
	CMatrix& operator+=(const CMatrix&);//矩阵的+=操作
	CMatrix operator + (const CMatrix& A);//矩阵的+操作
};

