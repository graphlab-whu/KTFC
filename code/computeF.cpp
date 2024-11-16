#include <bits/stdc++.h>
using namespace std;
using namespace chrono;
const int TMAX = 32000000;
int T[32000000];
int tn;
int qt;
void generateRandomNumbers(int tn)
{
	std::set<int> uniqueNumbers; 

	std::random_device rd;							  
	std::mt19937 gen(rd());							  
	std::uniform_int_distribution<int> dist(1, TMAX); 

	while (uniqueNumbers.size() < tn)
	{

		int randomNum = dist(gen); 
		uniqueNumbers.insert(randomNum);
	}

	// 将唯一的随机数存入数组 T
	int index = 0;
	for (int num : uniqueNumbers)
	{
		T[index++] = num;
	}
}

double checkf_baseline()
{
	double ans = 0.0;
	for (int i = 0; i < tn; ++i)
	{
		for (int j = i + qt - 1; j < tn; ++j)
		{
			ans = max(ans, (j - i + 1) * 1.0 / (T[j] - T[i] + 1) * 1.0);
		}
	}
	return ans;
}
int dq[TMAX];
int isLCP(int ida, int idb, int idc, long long a, long long b, long long c)
{
	return (idb - ida) * (c - a) - (idc - ida) * (b - a);
}
int isT(int ida, int idb, int idc, long long a, long long b, long long c)
{
	return (idc - idb + 1) * (c - a + 1) - (idc - ida + 1) * (c - b + 1);
}
double checkf_opt()
{

	int l = 0;
	int r = 0;
	double ans = 0x7f7f7f7f * 1.0;
	for (int i = qt - 1; i < tn; ++i)
	{
		while (r - l > 1 && isLCP(dq[r - 2], dq[r - 1], i - qt + 1, T[dq[r - 2]], T[dq[r - 1]], T[i - qt + 1]) >= 0)
		{
			r--;
		}
		dq[r++] = i - qt + 1;
		while (r - l > 1 && isT(dq[l], dq[l + 1], i, T[dq[l]], T[dq[l + 1]], T[i]) > 0)
		{
			++l;
		}
		double res = 1.0 * (T[i] - T[dq[l]] + 1) / (i - dq[l] + 1);
		ans = min(ans, res);
	}

	return 1.0 / ans;
}
int main(int argc, char *argv[])
{
	tn = stoi(argv[1]);
	generateRandomNumbers(tn);
	qt = stoi(argv[2]);
	auto start1 = system_clock::now();
	double ans1;
	for(int i=1;i<=1000;++i) ans1=checkf_baseline();
	auto end1 = system_clock::now();
	auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
	cout << fixed << setprecision(6) << static_cast<double>(duration1.count()) * 1.0 << "μs" << endl;
	auto start2 = system_clock::now();
	double ans2;
	for(int i=1;i<=1000;++i)ans2 = checkf_opt();
	auto end2 = system_clock::now();
	auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
	cout << fixed << setprecision(6) << static_cast<double>(duration2.count()) * 1.0 << "μs" << endl;
	cout << ans1 << " " << ans2 << endl;
	return 0;
}
