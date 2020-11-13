#include <iostream>
#include "CMatrix.h"

using namespace std;

int main()
{
	CMatrix A(2, 1,2);//
	CMatrix B(1, 2,3);
	CMatrix C(2, 1);
	//C=B + A;
	//cin >> B;
	cout << "A=" << A << endl;
	B=A.T();
	cout <<"B="<< B << endl;
	system("pause");
}