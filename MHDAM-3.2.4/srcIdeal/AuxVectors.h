#ifndef _AUX_VECTORS_H_
#define _AUX_VECTORS_H_


struct int_vector3d
{
	int m_iData[3];
	int_vector3d(int* data)
	{
		m_iData[0]=data[0];m_iData[1]=data[1];m_iData[2]=data[2];
		//std::sort(m_iData,m_iData+3);
	}

	int_vector3d(int v1,int v2,int v3 = 0)
	{
		m_iData[0]=v1;m_iData[1]=v2;m_iData[2]=v3;
		//std::sort(m_iData,m_iData+3);
	}

	void set_data(int v1,int v2,int v3 = 0)
	{
		m_iData[0]=v1;m_iData[1]=v2;m_iData[2]=v3;
		//std::sort(m_iData,m_iData+3);
	}
};

bool operator < (const int_vector3d& vec1,const int_vector3d& vec2)
{
	if (vec1.m_iData[0]<vec2.m_iData[0]) return true;
	if (vec1.m_iData[0]>vec2.m_iData[0]) return false;

	if (vec1.m_iData[1]<vec2.m_iData[1]) return true;
	if (vec1.m_iData[1]>vec2.m_iData[1]) return false;

	if (vec1.m_iData[2]<vec2.m_iData[2]) return true;
	if (vec1.m_iData[2]>vec2.m_iData[2]) return false;
	
	return false;
}

bool operator == (const int_vector3d& vec1,const int_vector3d& vec2)
{
	return 
		(vec1.m_iData[0]==vec2.m_iData[0])&&
		(vec1.m_iData[1]==vec2.m_iData[1])&&
		(vec1.m_iData[2]==vec2.m_iData[2]);	
}

#endif

