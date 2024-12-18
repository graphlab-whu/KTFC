#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
#include <chrono>
#include <array>
#include <queue>
#include <cmath>
// #include<cstring>
#include <string>
#include <iomanip>
#include <stack>
using namespace std;
using namespace chrono;
template <typename T>
void hash_combine(size_t &seed, T const &v)
{
	seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct pair_hash
{
	template <typename T1, typename T2>
	size_t operator()(std::pair<T1, T2> const &p) const
	{
		size_t seed = 0;
		hash_combine(seed, p.first);
		hash_combine(seed, p.second);
		return seed;
	}
};

// const int VMAX = 2610000;
const int VMAX = 4000000;
const int EMAX = 70000000;
const int TMAX = 239800000;
const int feps = 100000;
int vern = 0;
int arcn = 0;
int tmax = 0;
int mxt = 0;
int mxk = 0;
vector<vector<vector<pair<int, int>>>> CFI;

vector<map<int, int>> skyline;

int getCF(int v, int k, int t)
{
	if (CFI.size() < v + 1 || CFI[v].size() < k + 1 || CFI[v][k].empty())
		return -feps - 2;
	auto it = upper_bound(CFI[v][k].begin(), CFI[v][k].end(), t,
						  [](int value, const pair<int, int> &element)
						  {
							  return value < element.first;
						  });
	return (--it)->second;
}
int *isRes = new int[VMAX];
int lb1[feps + 1];
int main(int argc, char *argv[])
{

	char *graph = argv[1];
	freopen((string("../result/CFI_") + string(graph)).c_str(), "r", stdin);
	string ii, num12;
	string kk;
	CFI.push_back(vector<vector<pair<int, int>>>());
	while (getline(cin, ii))
	{
		CFI.push_back(vector<vector<pair<int, int>>>());
		CFI[stoi(ii)].push_back(vector<pair<int, int>>());
		while (getline(cin, kk))
		{
			CFI[stoi(ii)].push_back(vector<pair<int, int>>());
			if (kk.length() == 0)
				break;
			while (getline(cin, num12))
			{
				size_t pos = num12.find(' ');
				string firstNumStr = num12.substr(0, pos);
				string secondNumStr = num12.substr(pos + 1);
				int num1 = stoi(firstNumStr);
				int num2 = stoi(secondNumStr);
				int k = stoi(kk);
				mxk = max(mxk, k);
				int i = stoi(ii);
				vern = max(vern, i);
				CFI[i][k].push_back({num1, num2});
				mxt = max(mxt, num1 - 1);
				if (num2 == -feps)
					break;
			}
		}
	}
	//		cout<<mxk<<endl;
	for (int i = 1; i <= mxk + 1; ++i)
	{
		skyline.push_back(map<int, int>());
	}
	auto start = system_clock::now();
	for (int j = mxt; j; --j)
	{
		int lb2 = 0;
		if(j%10==0) 
		cout<<j<<endl; 
		for (int f = feps; f >= 0; --f)
		{
			int klb = max(lb1[f], lb2) + 1;
			int k = klb;
			int flag = 0;
			if (k > mxk)
				break;
			for (int kk = klb; kk <= mxk; ++kk)
			{
				flag = 0;
				for (int i = 1; i <= vern; ++i)
				{
					if (getCF(i, kk, j) >= f)
					{
						flag = 1;
						break;
					}
				}
				if (flag == 0)
				{
					k = kk - 1;
					break;
				}
			}

			lb2 = max(lb2, k);
			lb1[f] = max(lb1[f], k);
			if (k >= klb)
			{
				skyline[k][j] = f;
			}
		}
	}
	auto end = system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	cout << fixed << setprecision(6) << static_cast<double>(duration.count()) * 1.0 << "ms";
	freopen(("../result/" + string("Skyline_baseline_") + string(argv[1])).c_str(), "w", stdout);
	for (int k = 1; k <= mxk; ++k)
	{
		for (auto x : skyline[k])
		{
			cout << k << " " << x.first << " " << fixed << setprecision(6) << x.second * 1.0 / feps << endl;
		}
	}
	return 0;
}
