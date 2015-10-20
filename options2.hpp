#ifndef OPTION_HEADER_RKROTATE
#define OPTION_HEADER_RKROTATE
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <algorithm>

using std::cerr;
using std::endl;

class Options
{
private:
	std::vector<int> integers;
	std::vector<int> angleSums;
	std::vector<int> ringAtoms;
	std::vector<int> angAtoms;
	std::vector<double> lengths;
	std::vector<double> equiang;
	bool workingStatus;
	bool ringCalculating;
	bool allIn;
	double k_value;
	double a_value;
        
public:
	Options(std::istream &is)
		: workingStatus(false)
		, ringCalculating(false)
		, allIn(false)
	{
		if (!is)
		{
			cerr << "Open file Error!!!" << __LINE__ << endl;
			exit(0);
		}
		std::string str;
		if (!(is >> str))
		{
			cerr << "Format Error!!!" << __LINE__ << endl;
			exit(0);
		}
		bool error_flag = true;
		if (str == "0")
		{
			workingStatus = true;
			error_flag = false;
		}
		if (str == "1")
		{
			workingStatus = false;
			error_flag = false;
		}
		if (str == "2")
		{
			workingStatus = true;
			form(is);
			error_flag = false;
		}
		if (str == "3")
		{
			workingStatus = false;
			form(is);
			error_flag = false;
		}
		if (error_flag)
		{
			cerr << "Format Error!!!" << __LINE__ << endl;
			cerr << "Accepted:" << str << endl;
			exit(0);
		}
		int count;
		int step=workingStatus?4:1;
		bool isFormingAtoms=true;
		int turningK=step;
		while (is>>count)
		{
			if (isFormingAtoms)
			{
				if (count<0)
				{
					allIn = true;
					break;
				}
				integers.push_back(count);
				if (!(--turningK))
				{
					turningK=2;
					isFormingAtoms=false;
				}
			}
			else
			{
				angleSums.push_back(count);
				if (!(--turningK))
				{
					turningK=step;
					isFormingAtoms=true;
				}
			}
			
		}
	}
	inline const std::vector<int>& getResult()
	{
		return integers;
	}
	inline const std::vector<int>& getAngles()
	{
		return angleSums;
	}
	inline const std::vector<int>& getAtoms()
	{
		return ringAtoms;
	}
	inline const std::vector<int>& getangAtoms()
	{
		return angAtoms;
	}
	inline const std::vector<double>& getLengths()
	{
		return lengths;
	}
	inline const std::vector<double>& getequiang()
	{
		return equiang;
	}
	inline const double getK()
	{
		return k_value;
	}
	inline const double getA()
	{
		return a_value;
	}
	inline bool getRingCalculating()
	{
		return ringCalculating;
	}
	inline bool getWorkingStatus()
	{
		return workingStatus;
	}
	inline bool getAllInStatus()
	{
		return allIn;
	}
private:
	void form(std::istream &ost)
	{
		ringCalculating = true;
		ringAtoms.clear();
                angAtoms.clear();
		int size = 0;
		ost >> k_value; // constant of bond
                ost >> a_value; // constan of angle
		ost >> size;  // (size/3) is the number of broken bond 
              
              // read in the atoms of broken bond
		for (int i=0;i<size;++i)
		{
			if (i % 3 == 2)
			{
				double db;
				ost >> db;
				lengths.push_back(db); // l0 of bond
			}
			else
			{
				int id;
				ost >> id;
				ringAtoms.push_back(id);
			}
		}
               //-------------------------------------
              // read in the atoms forming the angles 
             for (int i=0;i<size;++i)
                {
                        if (i % 3 == 2)
                        {
                                double ag;
                                ost >> ag;
                                equiang.push_back(ag); // v0 of bond
                        }
                        else
                        {
                                int tt;
                                ost >> tt;
                                angAtoms.push_back(tt);
                        }
                }

	}
};
#endif
