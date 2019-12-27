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
	vector<vector<vector<double>>> m_Ti;		//训练图像			顺序Z X Y
	vector<vector<vector<double>>> m_Sim;		//保存采样结果		顺序Z X Y
	vector<C3DPoint> m_Samples;					//初始样品数据

	bool IsHaveTi;						//目前是否有训练图像		
	bool IsHaveSamples;					//目前是否有样品

	bool IsUseSamples;					//是否使用样品

	bool IsSimulation;						//目前是否采样结束		

											//	int *p_PathSim;					//模拟路径
											//vector<int> p_PathSim;         //模拟路径
	int *p_PathSim;					//模拟路径
	int m_SimX, m_SimY, m_SimZ;				//采样结果的尺寸
	int m_TiX, m_TiY, m_TiZ;				//训练图像的尺寸

	int m_SamplesMinX, m_SamplesMaxX;		//初始样品范围
	int m_SamplesMinY, m_SamplesMaxY;
	int m_SamplesMinZ, m_SamplesMaxZ;

	double m_SearchRadius;  //搜索半径，用于找到数据样板中有效的n个点作为数据事件
	int m_MaxPoint;                //最大的有效点个数

	float m_F;					//每次在训练图像上扫描的百分比,(0,1]
	double m_T;				//对比距离的阈值

	int *p_ScanCouts;			//每次接收的时候扫描训练图像的次数////////////////////////////////可能需要修改
	int *p_Bestmin;				//每次接收的时候的最小距离////////////////////////////////////////可能需要修改
	double m_Time;				//采样耗时
	double m_Time_omp;
	//search_radius是搜索半径，f是每次在训练图像上扫描的百分比(0,1]，thr是对比距离的阈值，isUseSamples是否有样品
	/*int simx, int simy, int simz, double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples = true*/
	DSsimulation(double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples);
	//search_radius是搜索半径，f是每次在训练图像上扫描的百分比(0,1]，thr是对比距离的阈值，isUseSamples是否有样品
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
	//	int _xSim, _ySim, _zSim;		//当前模拟点
	//	int _x0, _x1;			//数据事件的X范围（x0是起，x1是止）
	//	int _y0, _y1;			//数据事件的Y范围（y0是起，y1是止）
	//	int _z0, _z1;			//数据事件的Z范围（z0是起，z1是止）

	string InsertSamples();
	int *SetPath(int X, int Y, int Z);
	void SetPath(vector<int> & Path, int X, int Y, int Z);
	int Compare(C3DPoint a, C3DPoint b);
	void InsertEffectivePoint(vector<C3DPoint> &EffectivePoint, C3DPoint pointz, int &_xSim, int &_ySim, int &_zSim, int &_x0, int &_x1, int &_y0, int &_y1, int &_z0, int &_z1);
	//得到有效点，三维八分
	//	void GetEffectivePoint(vector<C3DPoint> &EffectivePoint, int &_xSim, int &_ySim, int &_zSim);
	//数据事件的距离
	void GetEffectivePoint(vector<C3DPoint>& EffectivePoint, int &_xSim, int &_ySim, int &_zSim, int &_x0, int &_x1, int &_y0, int &_y1, int &_z0, int &_z1);
	double GetDistances(vector<C3DPoint> EffectivePoint, int x, int y, int z, int &_xSim, int &_ySim, int &_zSim);
};