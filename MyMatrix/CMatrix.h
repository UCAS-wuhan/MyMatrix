#pragma once
#include <string.h>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
//先写普通的计算函数
//再把计算函数做成函数模板
//把计算函数放到类中（查一下类的构造、析构、深拷贝）
//把类做成模板（查一下类模板和函数模板区别用法）

//矩阵类，声明矩阵的定义及处理函数
class CMatrix
{
private:
	int  m_nRows; // 矩阵的行数
	int  m_nColumns; // 矩阵的列数
	double* m_lpBuf; // 动态分配用来存放数组的空间
public:
	CMatrix(int, int);
	CMatrix(int, int, double);	//构造时初始化为指定值
	CMatrix(const CMatrix& A);//深拷贝构造函数
	virtual ~CMatrix();//不作为基类时应去掉关键字

	double& CMatrix::operator()(int i, int j)const;
	double* operator[](int i) { return m_lpBuf + (i-1)*m_nColumns; }//注意this加括号， (*this)[i][j]
	void resize(int, int);//重新分配空间
	friend std::istream &operator >> (std::istream&,CMatrix&);
	friend std::ostream &operator << (std::ostream&, CMatrix&);     // 输出到屏幕

	CMatrix& operator = (const CMatrix&);//矩阵的复制
	CMatrix& operator+=(const CMatrix& );//矩阵的+=操作
	CMatrix operator + (const CMatrix& A);//+
	CMatrix& operator-=(const CMatrix&);//-=
	CMatrix operator - (const CMatrix& A);//-
	friend CMatrix  operator*(const CMatrix&, const CMatrix&);  //两矩阵相乘
	friend CMatrix  operator*(double, const CMatrix&);  //单数乘矩阵
	friend CMatrix  operator*(const CMatrix&, double);  //矩阵乘单数
	friend CMatrix  operator/(const CMatrix&, double);  //矩阵 除以单数

	CMatrix multi(const CMatrix&); // 对应元素相乘
	int row()const { return m_nRows; }	//返回行号
	int col()const { return m_nColumns; } //返回列号
	CMatrix getrow(int index); // 返回第index 行,索引从1算起
	CMatrix getcol(int index); // 返回第index 列,索引从1算起

	
	//CMatrix cov(_In_opt_ bool flag = true);   //协方差阵 或者样本方差
	//double det();   //行列式
	//CMatrix solveAb(CMatrix &obj);  // b是行向量或者列向量
	//CMatrix diag();  //返回对角线元素
	//				//Matrix asigndiag();  //对角线元素
	//CMatrix T()const;   //转置
	//void sort(bool);//true为从小到大
	//CMatrix adjoint();
	//CMatrix inverse();
	//void QR(_Out_ CMatrix&, _Out_ CMatrix&)const;
	//CMatrix eig_val(_In_opt_ Index_T _iters = 1000);
	//CMatrix eig_vect(_In_opt_ Index_T _iters = 1000);

	//double norm1();//1范数
	//double norm2();//2范数
	//double mean();// 矩阵均值
	//double*operator[](int i) { return m_lpBuf + i*m_nColumns; }//注意this加括号， (*this)[i][j]
	//void zeromean(_In_opt_  bool flag = true);//默认参数为true计算列
	//void normalize(_In_opt_  bool flag = true);//默认参数为true计算列
	//CMatrix exponent(double x);//每个元素x次幂
	//CMatrix  eye();//对角阵
	//void  maxlimit(double max, double set = 0);//对角阵

};

class Bounds_error                                          //下标越界异常
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

class RCdismatch                                           //行列不匹配异常
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
