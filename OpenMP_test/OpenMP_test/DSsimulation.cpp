#include "stdafx.h"
#include "DSsimulation.h"
#define MAX_THREAD_NUM 10

/*int simx, int simy, int simz, double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples /*= true*/
DSsimulation::DSsimulation(double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples /*= true*/)
{
	/*	m_SimX = simx;
	m_SimY = simy;
	m_SimZ = simz;*/

	m_SearchRadius = search_radius;
	m_MaxPoint = max_pointsCount;
	m_F = f;
	m_T = thr;
	IsUseSamples = isUseSamples;

	IsSimulation = false;
	m_TiX = m_TiY = m_TiZ = 0;
	IsHaveTi = false;
	IsHaveSamples = false;
	m_Sim.clear();

	//3D
	LoadTi("D:\\VSProject\\OpenMP_test\\data\\fold_categorical_180x150x120.SGEMS");
	LoadSamples("D:\\VSProject\\OpenMP_test\\data\\543_xyz.sgems");

	for (int i = 0; i < m_SimZ; i++)
	{
		vector<vector<double>> temp;
		for (int j = 0; j < m_SimX; j++)
		{
			vector<double> temp0;
			for (int k = 0; k < m_SimY; k++)
			{
				temp0.push_back(-1);
			}
			temp.push_back(temp0);
		}
		m_Sim.push_back(temp);
	}
	p_ScanCouts = new int[m_SimX*m_SimY*m_SimZ];
	p_Bestmin = new int[m_SimX*m_SimY*m_SimZ];
}


DSsimulation::~DSsimulation()
{
	if (p_ScanCouts != NULL)
	{
		delete p_ScanCouts;
		p_ScanCouts = NULL;
	}
	if (p_Bestmin != NULL)
	{
		delete p_Bestmin;
		p_Bestmin = NULL;
	}
}

/*int simx, int simy, int simz, double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples*/
void DSsimulation::SetParameter(double search_radius, int max_pointsCount, float f, double thr, bool isUseSamples)
{
	/*
	m_SimX = simx;
	m_SimY = simy;
	m_SimZ = simz;
	*/
	m_SearchRadius = search_radius;
	m_MaxPoint = max_pointsCount;
	m_F = f;
	m_T = thr;
	IsUseSamples = isUseSamples;

	m_Sim.clear();
	for (int i = 0; i < m_SimZ; i++)
	{
		vector<vector<double>> temp;
		for (int j = 0; j < m_SimX; j++)
		{
			vector<double> temp0;
			for (int k = 0; k < m_SimY; k++)
			{
				temp0.push_back(-1);
			}
			temp.push_back(temp0);
		}
		m_Sim.push_back(temp);
	}
	p_ScanCouts = new int[m_SimX*m_SimY*m_SimZ];
	p_Bestmin = new int[m_SimX*m_SimY*m_SimZ];
}

bool DSsimulation::LoadTi(string FilePath)
{
	m_Ti.clear();
	m_TiX = m_TiY = m_TiZ = 0;
	IsHaveTi = false;

	ifstream ifile;
	ifile.open(FilePath);
	if (!ifile)
	{
		cout << "训练图像载入错误！" << endl;
		return false;
	}

	string temp;
	ifile >> m_TiY >> m_TiX >> m_TiZ;
	ifile >> temp;
	if (temp != "1")
	{
		return false;
	}
	ifile >> temp;
	/*	if (temp != "facies")
	{
	return false;
	}*/

	for (int k = 0; k < m_TiZ; k++)
	{
		vector<vector<double>> xyTi;
		for (int i = 0; i < m_TiX; i++)
		{
			vector<double> yTi;
			for (int j = 0; j < m_TiY; j++)
			{
				double values;
				ifile >> values;

				yTi.push_back(values);
			}
			xyTi.push_back(yTi);
		}
		m_Ti.push_back(xyTi);
	}
	ifile.close();
	IsHaveTi = true;
	cout << "训练图像载入成功！" << endl;
	return true;
}

bool DSsimulation::LoadSamples(string FilePath)
{
	m_Samples.clear();
	IsHaveSamples = false;
	m_SamplesMinX = m_SamplesMinY = m_SamplesMinZ = INT_MAX;
	m_SamplesMaxX = m_SamplesMaxY = m_SamplesMaxZ = -1;

	ifstream ifile;
	ifile.open(FilePath);
	if (!ifile)
	{
		cout << "样品图像载入失败！" << endl;
		return false;
	}

	string temp;
	int Couts;
	ifile >> Couts;
	ifile >> temp;
	if (temp != "4")
	{
		return false;
	}
	ifile >> temp;
	if (temp != "x")
	{
		return false;
	}
	ifile >> temp;
	if (temp != "y")
	{
		return false;
	}
	ifile >> temp;
	if (temp != "z")
	{
		return false;
	}
	ifile >> temp;
	/*	if (temp != "sample")
	{
	return false;
	}*/

	for (int i = 0; i < Couts; i++)
	{
		C3DPoint tempPoint;
		ifile >> tempPoint.x >> tempPoint.y >> tempPoint.z >> tempPoint.values;
		/*
		tempPoint.x--;
		tempPoint.y--;
		tempPoint.z--;*/
		m_Samples.push_back(tempPoint);

		m_SamplesMinX = m_SamplesMinX < tempPoint.x ? m_SamplesMinX : tempPoint.x;
		m_SamplesMaxX = m_SamplesMaxX > tempPoint.x ? m_SamplesMaxX : tempPoint.x;
		m_SamplesMinY = m_SamplesMinY < tempPoint.y ? m_SamplesMinY : tempPoint.y;
		m_SamplesMaxY = m_SamplesMaxY > tempPoint.y ? m_SamplesMaxY : tempPoint.y;
		m_SamplesMinZ = m_SamplesMinZ < tempPoint.z ? m_SamplesMinZ : tempPoint.z;
		m_SamplesMaxZ = m_SamplesMaxZ > tempPoint.z ? m_SamplesMaxZ : tempPoint.z;
	}
	/////////////////////////////////////////////////////////////////////////////////////
	m_SimX = m_SamplesMaxX + 1;
	m_SimY = m_SamplesMaxY + 1;
	m_SimZ = m_SamplesMaxZ + 1;
	/////////////////////////////////////////////////////////////////////////////////////
	ifile.close();
	IsHaveSamples = true;
	cout << "样品图像载入成功！" << endl;
	return true;
}

bool DSsimulation::SaveSimulation(int i)
{
	ofstream ofile;
	string FilePath = "D:\\VSProject\\OpenMP_test\\";
	string FileName = "openmpNEW2D.txt";
	FileName[6] += i;
	string File = FilePath + FileName;
	ofile.open(FileName);
	if (!ofile)
	{
		cout << "文件打开失败！" << endl;
		return false;
	}
	ofile << m_SimY << " " << m_SimX << " " << m_SimZ << endl;
	ofile << "1" << endl;
	ofile << "code" << endl;
	for (int k = 0; k < m_SimZ; k++)
	{
		for (int i = 0; i < m_SimX; i++)
		{
			for (int j = 0; j < m_SimY; j++)
			{
				ofile << m_Sim[k][i][j] << endl;
			}
		}
	}
	cout << "保存成功！" << endl;
	ofile.close();
	cout << "time:" << m_Time << endl;
	cout << "omp time" << m_Time_omp << endl;
	return true;
}

bool DSsimulation::IsConflict(int *current_simulate_position, int my_rank, int num_threads)
{
	int my_position = p_PathSim[current_simulate_position[my_rank]];
	int my_x, my_y, my_z;
	my_z = my_position / (m_SimX * m_SimY);
	int temp = my_position % (m_SimX * m_SimY);
	my_x = temp / m_SimY;
	my_y = temp % m_SimY;
	for (int i = 0;i < num_threads;i++)
	{
		if (i == my_rank || current_simulate_position[i] == -1)
			continue;
		int i_position = p_PathSim[current_simulate_position[i]];
		int i_x, i_y, i_z;
		i_z = i_position / (m_SimX * m_SimY);
		int temp = i_position % (m_SimX * m_SimY);
		i_x = temp / m_SimY;
		i_y = temp % m_SimY;
		if (abs(my_x - i_x) <= m_SearchRadius && abs(my_y - i_y) <= m_SearchRadius && abs(my_z - i_z) <= m_SearchRadius)
		{
			return true;
		}
	}
	return false;
}
/*
void DSsimulation::swap(int i, int j)
{
if (i == j)
return;
int temp = p_PathSim[i];
p_PathSim[i] = p_PathSim[j];
p_PathSim[j] = temp;
return;
}
*/
std::string DSsimulation::StartSimulation()
{
	ResetSim();
	if (!IsHaveTi)
	{
		cout << "没有训练图像，请先加载训练图像！！！" << endl;;
		return NULL;
	}
	if (!IsHaveSamples)
	{
		cout << "没有加载样品，请先加载样品数据！！！" << endl;;
		return NULL;
	}

	clock_t begin, end;


	if (IsUseSamples)
	{
		string results = InsertSamples();
		if (results != "加载成功！")
		{
			cout << results;
			return "加载样品失败," + results;
		}
	}
	//	SYSTEM_INFO SystemInfo;
	//	GetSystemInfo(&SystemInfo);
	//	cout << "cpu num:" << SystemInfo.dwNumberOfProcessors << endl;

	begin = clock();
	double begin_omp = omp_get_wtime();

	int sim_size = m_SimX*m_SimY*m_SimZ;//模拟图像的总尺寸
	p_PathSim = SetPath(m_SimX, m_SimY, m_SimZ);
	//SetPath(p_PathSim, m_SimX, m_SimY, m_SimZ);
	/*	double distances = 0;
	double sum = 0;
	double Thr = 0;*/

	srand((unsigned)time(NULL));
	//////////////////////////////////////////////////////////////////设置线程数
	int m_thread_num = 10;
	int current_sim_position[MAX_THREAD_NUM];
	for (int i = 0;i < m_thread_num;i++)
		current_sim_position[i] = -1;
	int local_sim_num = sim_size / m_thread_num;

#pragma omp parallel num_threads(m_thread_num)
	{
		int _xSim, _ySim, _zSim;//当前模拟点
		int _x0 = -1, _x1 = -1, _y0 = -1, _y1 = -1, _z0 = -1, _z1 = -1;//数据事件的起始范围
		double distances = 0;
		double sum = 0;
		double Thr = 0;
		int my_rank = omp_get_thread_num();
		//		SetThreadAffinityMask(GetCurrentThread(),pow(2,my_rank) );
		vector<int> sim_points;//需要模拟的点，值为模拟路径的下标
							   //		vector<bool> sim_flag;//标记点是否需要模拟，当为true时需要重新模拟或已模拟
		vector<int> resim_points;//需要重新模拟的点
		int n = 5;//冲突次数超过n时需要重新模拟
		sim_points.clear();
		resim_points.clear();
		int start0 = my_rank*local_sim_num;//每个线程开始位置
		int end0 = (my_rank + 1)*local_sim_num;//每个线程结束位置
#pragma omp parallel for
		for (int i = start0;i < end0;i++)
		{

			{
				sim_points.push_back(i);
				//				sim_flag.push_back(false);
				//				cout << "线程" << my_rank << "：" << i << endl;
			}
		}
#pragma omp critical
		cout << "thread" << my_rank << endl;
#pragma omp barrier
		{
#pragma omp critical
			cout << "thread" << my_rank << "start sim" << endl;
			for (int i = 0;i < sim_points.size();i++)
			{
				if (i % 100 == 0)
					cout << "sim线程号：" << my_rank << " 当前第" << i << endl;
				int conflict_count = 0;
				{
					for (conflict_count = 0;conflict_count < n&&conflict_count < (sim_points.size() - i);conflict_count++)
					{
						current_sim_position[my_rank] = sim_points[i + conflict_count];
						if (!IsConflict(current_sim_position, my_rank, m_thread_num))
						{
							break;
						}
					}
				}

				if (conflict_count == n || conflict_count == (sim_points.size() - i))
				{

					//					resim_points.push_back(sim_points[i + conflict_count - 1]);
					/////////////////////////////////////////////////////
					swap(sim_points[i], sim_points[i + conflict_count - 1]);
				}
				else
				{
					int position_in_PathSim = sim_points[i + conflict_count];
					{
						double mindist = DBL_MAX;
						int bestX, bestY, bestZ;
						{
							_zSim = p_PathSim[position_in_PathSim] / (m_SimX * m_SimY);
							int tempSim = p_PathSim[position_in_PathSim] % (m_SimX * m_SimY);
							_xSim = tempSim / m_SimY;
							_ySim = tempSim % m_SimY;
						}

						if (m_Sim[_zSim][_xSim][_ySim] == -1)
						{
							//					cout << "sim线程号：" << my_rank << "else if" << endl;
							vector<C3DPoint> EffectivePoint;

							GetEffectivePoint(EffectivePoint, _xSim, _ySim, _zSim, _x0, _x1, _y0, _y1, _z0, _z1);


							if (EffectivePoint.size() == 0)
							{
								//						cout << "sim线程号：" << my_rank << "else if if" << endl;
								m_Sim[_zSim][_xSim][_ySim] = m_Ti[rand() % m_TiZ][rand() % m_TiX][rand() % m_TiY];

								p_ScanCouts[position_in_PathSim] = 0;
								p_Bestmin[position_in_PathSim] = 0;

								Thr = sum / (position_in_PathSim + 1)*m_T;
							}
							else
							{
								//						cout << "sim线程号：" << my_rank << "else if else" << endl;
								int scanX = m_TiX - _x0 - _x1;
								int scanY = m_TiY - _y0 - _y1;
								int scanZ = m_TiZ - _z0 - _z1;
								//						cout << "sim线程号：" << my_rank << " 1111111111111111111111111111111111111" << endl;
								//						vector<int> pathTi;
								//						pathTi.clear();
								//						int *pathTi = SetPath(scanX, scanY, scanZ);
								int *pathTi = SetPath(scanX, scanY, scanZ);
								//						SetPath(pathTi, scanX, scanY, scanZ);//vector为什么会有问题？
								int j = 0;
								//						cout << "sim线程号：" << my_rank << " 2222222222222222222222222222222222222" << endl;
								for (j = 0; j < scanX*scanY*scanZ*m_F; j++)
								{
									int xTi = pathTi[j] / (scanY * scanZ) + _x0;
									int tempTi = pathTi[j] % (scanY * scanZ);
									int yTi = tempTi / scanZ + _y0;
									int zTi = tempTi % scanZ + _z0;
									distances = GetDistances(EffectivePoint, xTi, yTi, zTi,_xSim,_ySim,_zSim);
									if (distances < mindist)
									{
										mindist = distances;
										bestX = xTi;
										bestY = yTi;
										bestZ = zTi;
									}

									if (mindist <= Thr*(j / scanX*scanY*scanZ*m_F))
									{
										j++;
										//								cout << "sim线程号：" << my_rank << "break at" <<j<< endl;
										break;
									}
								}
								//						cout << "sim线程号：" << my_rank << " 333333333333333333333333333333" << endl;
								if (pathTi)
									delete pathTi;


								m_Sim[_zSim][_xSim][_ySim] = m_Ti[bestZ][bestX][bestY];

								p_ScanCouts[position_in_PathSim] = j;
								p_Bestmin[position_in_PathSim] = mindist;

								sum = sum + mindist;

								Thr = sum / (position_in_PathSim + 1)*m_T;
								//						cout << "sim线程号：" << my_rank << " 444444444444444444444444444444" << endl;
							}

							//						sim_flag[i + conflict_count] = true;
							//											cout << "thread" << my_rank << "end simulation" << endl;
						}
					}
					swap(sim_points[i], sim_points[i + conflict_count]);
					//						sim_flag[i + conflict_count] = true;
					//						cout << "thread" << my_rank << "end simulation" << endl;
				}
			}
		}

#pragma omp critical
		cout << "start resim" << endl;
#pragma omp barrier
#pragma omp parallel for
		for (int i = 0;i < sim_points.size();i++)
		{
#pragma omp critical
			cout << "resim thread：" << my_rank << " ：" << sim_points[i] << endl;
			int position_in_PathSim = sim_points[i];
			double mindist = DBL_MAX;
			int bestX, bestY, bestZ;
			{
				_zSim = p_PathSim[position_in_PathSim] / (m_SimX * m_SimY);
				int tempSim = p_PathSim[position_in_PathSim] % (m_SimX * m_SimY);
				_xSim = tempSim / m_SimY;
				_ySim = tempSim % m_SimY;
			}

			if (m_Sim[_zSim][_xSim][_ySim] != -1)
			{
				continue;
			}

			vector<C3DPoint> EffectivePoint;
			GetEffectivePoint(EffectivePoint, _xSim, _ySim, _zSim, _x0, _x1, _y0, _y1, _z0, _z1);


			if (EffectivePoint.size() == 0)
			{
				//				omp_set_lock(&lock);
				m_Sim[_zSim][_xSim][_ySim] = m_Ti[rand() % m_TiZ][rand() % m_TiX][rand() % m_TiY];
				//				omp_unset_lock(&lock);

				p_ScanCouts[position_in_PathSim] = 0;
				p_Bestmin[position_in_PathSim] = 0;

				Thr = sum / (position_in_PathSim + 1)*m_T;
			}
			else
			{
				int scanX = m_TiX - _x0 - _x1;
				int scanY = m_TiY - _y0 - _y1;
				int scanZ = m_TiZ - _z0 - _z1;

				int *pathTi = SetPath(scanX, scanY, scanZ);
				int j = 0;

				for (j = 0; j < scanX*scanY*scanZ*m_F; j++)
				{
					int xTi = pathTi[j] / (scanY * scanZ) + _x0;
					int tempTi = pathTi[j] % (scanY * scanZ);
					int yTi = tempTi / scanZ + _y0;
					int zTi = tempTi % scanZ + _z0;
				
					distances = GetDistances(EffectivePoint, xTi, yTi, zTi, _xSim, _ySim, _zSim);
					//bestX bestY bestZ都没有获取到？？？？？？？？？？？？
					if (distances < mindist)
					{
						mindist = distances;
						bestX = xTi;
						bestY = yTi;
						bestZ = zTi;
					}

					if (mindist <= Thr*(j / scanX*scanY*scanZ*m_F))
					{
						j++;
						break;
					}
				}
				if (pathTi)
					delete pathTi;

				//				omp_set_lock(&lock);
				m_Sim[_zSim][_xSim][_ySim] = m_Ti[bestZ][bestX][bestY];
				//				omp_unset_lock(&lock);

				p_ScanCouts[position_in_PathSim] = j;
				p_Bestmin[position_in_PathSim] = mindist;

				sum = sum + mindist;

				Thr = sum / (position_in_PathSim + 1)*m_T;
			}
		}
	}

	//p_PathSim.clear();
	if (p_PathSim)
	{
		delete p_PathSim;
		p_PathSim = NULL;
	}
	end = clock();
	double end_omp = omp_get_wtime();
	m_Time = end - begin;
	m_Time /= 1000;
	m_Time_omp = end_omp - begin_omp;
	IsSimulation = true;
	cout << "模拟完成！" << endl;
	return "成功";
}

void DSsimulation::ResetSim()
{
	IsSimulation = false;
	for (int i = 0; i < m_SimZ; i++)
	{
		for (int j = 0; j < m_SimX; j++)
		{
			for (int k = 0; k < m_SimY; k++)
			{
				m_Sim[i][j][k] = -1;
			}
		}
	}
}

std::string DSsimulation::InsertSamples()
{
	if (m_Sim.size() < 1)
	{
		return "sim的Z层初始化失败!";
	}
	else if (m_Sim[0].size() < 1)
	{
		return "sim的X层初始化失败!";
	}
	else if (m_Sim[0][0].size() < 1)
	{
		return "sim的Y层初始化失败!";
	}

	for (int i = 0; i < m_Samples.size(); i++)
	{
		if (m_Samples[i].z < m_Sim.size() && m_Samples[i].x < m_Sim[0].size() && m_Samples[i].y < m_Sim[0][0].size())
		{
			m_Sim[m_Samples[i].z][m_Samples[i].x][m_Samples[i].y] = m_Samples[i].values;
		}
		else
		{
			cout << i << "样品位置越界！" << endl;
			return "样品位置越界!";
		}
	}
	return "加载成功！";
}

int * DSsimulation::SetPath(int X, int Y, int Z)
{
	srand((unsigned)time(NULL));

	int sizeXYZ = X * Y * Z;
	int *p_Path = new int[sizeXYZ];
	for (int i = 0; i < sizeXYZ; i++)
		p_Path[i] = i;
	for (int i = 0; i < sizeXYZ; i++)
	{
		if (i < sizeXYZ)
		{
			int m = rand() % sizeXYZ;
			int temp = p_Path[i];
			p_Path[i] = p_Path[m];
			p_Path[m] = temp;
		}
	}
	return p_Path;
}

void DSsimulation::SetPath(vector<int> & Path, int X, int Y, int Z)
{
	for (long i = 0; i < X * Y * Z; i++)
	{
		Path.push_back(i);
	}

	srand((unsigned)time(NULL));
	std::random_shuffle(Path.begin(), Path.end());
}


int DSsimulation::Compare(C3DPoint a, C3DPoint b)
{
	if (a.distances < b.distances)
		return 1; //升序排列，如果改为 a >b，则为降序,要注意sort()中cmp()的返值只有1和0，不像qsort中存在－1！！！！
	else
		return 0;
}

void DSsimulation::InsertEffectivePoint(vector<C3DPoint> &EffectivePoint, C3DPoint point, int &_xSim, int &_ySim, int &_zSim, int &_x0, int &_x1, int &_y0, int &_y1, int &_z0, int &_z1)
{
	int x0 = _x0;
	int x1 = _x1;
	int y0 = _y0;
	int y1 = _y1;
	int z0 = _z0;
	int z1 = _z1;
	if (point.x <= _xSim)
	{
		if (_x0 < _xSim - point.x)
		{
			if (_xSim - point.x + _x1 < m_TiX)
			{
				x0 = _xSim - point.x;
			}
			else
			{
				return;
			}
		}
	}
	else
	{
		if (_x1 < point.x - _xSim)
		{
			if (point.x - _xSim + _x0 < m_TiX)
			{
				x1 = point.x - _xSim;
			}
			else
			{
				return;
			}
		}
	}

	if (point.y <= _ySim)
	{
		if (_y0 < _ySim - point.y)
		{
			if (_ySim - point.y + _y1 < m_TiY)
			{
				y0 = _ySim - point.y;
			}
			else
			{
				return;
			}
		}
	}
	else
	{
		if (_y1 < point.y - _ySim)
		{
			if (point.y - _ySim + _y0 < m_TiY)
			{
				y1 = point.y - _ySim;
			}
			else
			{
				return;
			}
		}
	}

	if (point.z <= _zSim)
	{
		if (_z0 < _zSim - point.z)
		{
			if (_zSim - point.z + _z1 < m_TiZ)
			{
				z0 = _zSim - point.z;
			}
			else
			{
				return;
			}
		}
	}
	else
	{
		if (_z1 < point.z - _zSim)
		{
			if (point.z - _zSim + _z0 < m_TiZ)
			{
				z1 = point.z - _zSim;
			}
			else
			{
				return;
			}
		}
	}
	{
		_x0 = x0;
		_x1 = x1;
		_y0 = y0;
		_y1 = y1;
		_z0 = z0;
		_z1 = z1;
	}
	EffectivePoint.push_back(point);
}

void DSsimulation::GetEffectivePoint(vector<C3DPoint>& EffectivePoint, int &_xSim, int &_ySim, int &_zSim, int &_x0, int &_x1, int &_y0, int &_y1, int &_z0, int &_z1)
{
	_x0 = _x1 = _y0 = _y1 = _z0 = _z1 = 0;
	int Radius = m_SearchRadius;
	int x0 = _xSim - Radius > 0 ? _xSim - Radius : 0;
	int x1 = _xSim + Radius < m_SimX - 1 ? _xSim + Radius : m_SimX - 1;
	int y0 = _ySim - Radius > 0 ? _ySim - Radius : 0;
	int y1 = _ySim + Radius < m_SimY - 1 ? _ySim + Radius : m_SimY - 1;
	int z0 = _zSim - Radius > 0 ? _zSim - Radius : 0;
	int z1 = _zSim + Radius < m_SimZ - 1 ? _zSim + Radius : m_SimZ - 1;
	for (int iRadius = 1; iRadius < Radius&&EffectivePoint.size() < m_MaxPoint; iRadius++)
	{
		int ix0 = _xSim - iRadius;
		int ix1 = _xSim + iRadius;
		int iy0 = _ySim - iRadius;
		int iy1 = _ySim + iRadius;
		int iz0 = _zSim - iRadius;
		int iz1 = _zSim + iRadius;
		//上面
		if (iz1 >= z0&&iz1 <= z1)
		{
			_z1 = iRadius;
			for (int i = ix0;i <= ix1;i++)
			{
				if (i >= x0&&i <= x1)
				{
					for (int j = iy0;j <= iy1;j++)
					{
						if (j >= y0&&j <= y1&&m_Sim[iz1][i][j] != -1 && EffectivePoint.size() < m_MaxPoint)
						{
							double distances = sqrt(double((i - _xSim)*(i - _xSim) + (j - _ySim)*(j - _ySim) + (iz1 - _zSim)*(iz1 - _zSim)));
							EffectivePoint.push_back(C3DPoint(i, j, iz1, m_Sim[iz1][i][j], distances));
						}
					}
				}
			}
		}
		//下面
		if (iz0 >= z0&&iz0 <= z1)
		{
			_z0 = iRadius;
			for (int i = ix0;i <= ix1;i++)
			{
				if (i >= x0&&i <= x1)
				{
					for (int j = iy0;j <= iy1;j++)
					{
						if (j >= y0&&j <= y1&&m_Sim[iz0][i][j] != -1 && EffectivePoint.size() < m_MaxPoint)
						{
							double distances = sqrt(double((i - _xSim)*(i - _xSim) + (j - _ySim)*(j - _ySim) + (iz0 - _zSim)*(iz0 - _zSim)));
							EffectivePoint.push_back(C3DPoint(i, j, iz0, m_Sim[iz0][i][j], distances));
						}
					}
				}
			}
		}
		//左面
		if (iy0 >= y0&&iy0 <= y1)
		{
			_y0 = iRadius;
			for (int k = iz0 + 1;k < iz1;k++)
			{
				if (k >= z0&&k <= z1)
				{
					for (int i = ix0;i < ix1;i++)
					{
						if (i >= x0&&i <= x1&&m_Sim[k][i][iy0] != -1 && EffectivePoint.size() < m_MaxPoint)
						{
							double distances = sqrt(double((i - _xSim)*(i - _xSim) + (iy0 - _ySim)*(iy0 - _ySim) + (k - _zSim)*(k - _zSim)));
							EffectivePoint.push_back(C3DPoint(i, iy0, k, m_Sim[k][i][iy0], distances));
						}
					}
				}
			}
		}
		//右面
		if (iy1 >= y0&&iy1 <= y1)
		{
			_y1 = iRadius;
			for (int k = iz0 + 1;k < iz1;k++)
			{
				if (k >= z0&&k <= z1)
				{
					for (int i = ix1;i > ix0;i--)
					{
						if (i >= x0&&i <= x1&&m_Sim[k][i][iy1] != -1 && EffectivePoint.size() < m_MaxPoint)
						{
							double distances = sqrt(double((i - _xSim)*(i - _xSim) + (iy1 - _ySim)*(iy1 - _ySim) + (k - _zSim)*(k - _zSim)));
							EffectivePoint.push_back(C3DPoint(i, iy1, k, m_Sim[k][i][iy1], distances));
						}
					}
				}
			}
		}
		//后面
		if (ix0 >= x0&&ix0 <= x1)
		{
			_x0 = iRadius;
			for (int k = iz0 + 1;k < iz1;k++)
			{
				if (k >= z0&&k <= z1)
				{
					for (int j = iy0;j < iy1;j++)
					{
						if (j >= y0&&j <= y1&&m_Sim[k][ix0][j] != -1 && EffectivePoint.size() < m_MaxPoint)
						{
							double distances = sqrt(double((ix0 - _xSim)*(ix0 - _xSim) + (j - _ySim)*(j - _ySim) + (k - _zSim)*(k - _zSim)));
							EffectivePoint.push_back(C3DPoint(ix0, j, k, m_Sim[k][ix0][j], distances));
						}
					}
				}
			}
		}
		//前面
		if (ix1 >= x0&&ix1 <= x1)
		{
			_x1 = iRadius;
			for (int k = iz0 + 1;k < iz1;k++)
			{
				if (k >= z0&&k <= z1)
				{
					for (int j = iy1;j > iy0;j--)
					{
						if (j >= y0&&j <= y1&&m_Sim[k][ix1][j] != -1 && EffectivePoint.size() < m_MaxPoint)
						{
							double distances = sqrt(double((ix1 - _xSim)*(ix1 - _xSim) + (j - _ySim)*(j - _ySim) + (k - _zSim)*(k - _zSim)));
							EffectivePoint.push_back(C3DPoint(ix1, j, k, m_Sim[k][ix1][j], distances));
						}
					}
				}
			}
		}
	}
	return;
}


double DSsimulation::GetDistances(vector<C3DPoint> EffectivePoint, int x, int y, int z, int &_xSim, int &_ySim, int& _zSim)
{
	double sum = 0;
	for (vector<C3DPoint>::iterator it = EffectivePoint.begin(); it != EffectivePoint.end(); it++)
	{
		if (abs(it->values - m_Ti[z - _zSim + it->z][x - _xSim + it->x][y - _ySim + it->y]) > 1e-6)
		{
			sum += 1.0;
		}
	}
	return sum / EffectivePoint.size();
}
