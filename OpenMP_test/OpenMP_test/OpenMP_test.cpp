// OpenMP_test.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "DSsimulation.h"


int main()
{
	DSsimulation ds(20, 20, 0.8, 0.1, true);
	for (int i = 0;i < 1;i++)
	{
		ds.StartSimulation();
		ds.SaveSimulation(i);
	}
	return 0;
}