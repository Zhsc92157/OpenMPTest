#include <time.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <functional>
#include <string>
#include<omp.h>
#include <windows.h>
using namespace std;

#pragma once

typedef struct C3DPointInformation
{
	int x;
	int y;
	int z;
	double values;
	double distances;

	C3DPointInformation() {}
	C3DPointInformation(int xx, int yy, int zz, double valuess, double dis)
	{
		x = xx;
		y = yy;
		z = zz;
		values = valuess;
		distances = dis;
	}
}C3DPoint;


class DSsimulation
{
public:
	vector<vector<vector<double>>> m_Ti;		//ѵ��ͼ��			˳��Z X Y
	vector<vector<vector<double>>> m_Sim;		//����������		˳��Z X Y
	vector<C3DPoint> m_Samples;					//��ʼ��Ʒ����

	bool IsHaveTi;						//Ŀǰ�Ƿ���ѵ��ͼ��		
	bool IsHaveSamples;					//Ŀǰ�Ƿ�����Ʒ

	bool IsUseSamples;					//�Ƿ�ʹ����Ʒ

	bool IsSimulation;						//Ŀǰ�Ƿ��������		

											//	int *p_PathSim;					//ģ��·��
											//vector<int> p_PathSim;         //ģ��·��
	int *p_PathSim;					//ģ��·��
	int m_SimX, m_SimY, m_SimZ;				//��������ĳߴ�
	int m_TiX, m_TiY, m_TiZ;				//ѵ��ͼ��ĳߴ�

	int m_SamplesMinX, m_SamplesMaxX;		//��ʼ��Ʒ��Χ
	int m_SamplesMinY, m_SamplesMaxY;
	int m_SamplesMinZ, m_SamplesMaxZ;

	double m_SearchRadius;  //�����뾶�������ҵ�������������Ч��n������Ϊ�����¼�
	int m_MaxPoint;                //������Ч�����

	float m_F;					//ÿ����ѵ��ͼ����ɨ��İٷֱ�,(0,1]
	double m_T;				//�ԱȾ������ֵ

	int *p_ScanCouts;			//ÿ�ν��յ�ʱ��ɨ��ѵ��ͼ��Ĵ���////////////////////////////////������Ҫ�޸�
	int *p_Bestmin;				//ÿ�ν��յ�ʱ�����С����////////////////////////////////////////������Ҫ�޸�
	double m_Time;				//������ʱ
	double m_Time_omp;
	//search_radius�������뾶��f��ÿ����ѵ��ͼ����ɨ��İٷֱ�(0,1]��thr�ǶԱȾ������ֵ��isUseSamples�Ƿ�����Ʒ
	/*int simx, int simy, int simz, double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples = true*/
	DSsimulation(double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples);
	//search_radius�������뾶��f��ÿ����ѵ��ͼ����ɨ��İٷֱ�(0,1]��thr�ǶԱȾ������ֵ��isUseSamples�Ƿ�����Ʒ
	void SetParameter(double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples);
	/*int simx, int simy, int simz, double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples*/
	bool LoadTi(string FilePath);
	bool LoadSamples(string FilePath);
	bool SaveSimulation(int i);
	string StartSimulation();
	bool IsConflict(int *current_simulate_position, int my_rank, int num_threads);
	void ResetSim();
	//	void swap(int i, int j);
	~DSsimulation();

private:
	//	int _xSim, _ySim, _zSim;		//��ǰģ���
	//	int _x0, _x1;			//�����¼���X��Χ��x0����x1��ֹ��
	//	int _y0, _y1;			//�����¼���Y��Χ��y0����y1��ֹ��
	//	int _z0, _z1;			//�����¼���Z��Χ��z0����z1��ֹ��

	string InsertSamples();
	int *SetPath(int X, int Y, int Z);
	void SetPath(vector<int> & Path, int X, int Y, int Z);
	int Compare(C3DPoint a, C3DPoint b);
	void InsertEffectivePoint(vector<C3DPoint> &EffectivePoint, C3DPoint pointz, int &_xSim, int &_ySim, int &_zSim, int &_x0, int &_x1, int &_y0, int &_y1, int &_z0, int &_z1);
	//�õ���Ч�㣬��ά�˷�
	//	void GetEffectivePoint(vector<C3DPoint> &EffectivePoint, int &_xSim, int &_ySim, int &_zSim);
	//�����¼��ľ���
	void GetEffectivePoint(vector<C3DPoint>& EffectivePoint, int &_xSim, int &_ySim, int &_zSim, int &_x0, int &_x1, int &_y0, int &_y1, int &_z0, int &_z1);
	double GetDistances(vector<C3DPoint> EffectivePoint, int x, int y, int z, int &_xSim, int &_ySim, int &_zSim);
};