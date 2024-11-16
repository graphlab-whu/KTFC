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

const int VMAX = 4000000;
const int EMAX = 70000000;
const int TMAX = 239800000;
const int feps = 100000;
int vern = 0;
int arcn = 0;
long long tmax = 0ll;
int mxk = 0;
int mxt = 0;
int choice = 0;
//unordered_map<int,int>node;
struct arc
{
	int src, dst;
	long long t;
} arcs[EMAX];

unordered_map<pair<int, int>, vector<long long>, pair_hash> tss;
int offset[VMAX];
int mxoffset[VMAX];
vector<vector<int>> G2(VMAX);
void loadgraph(const char *name)
{
	ifstream fin(name, ios::in);
	vector<int> v;
	long long tmin = 0x7fffffffll;

	string l;
	while (getline(fin, l))
	{
		//		cout<<l<<endl;
		int uvt[2] = {0};
		long long uvt3 = 0ll;
		int p = -1, np = -1;
		for (int i = 0; i < 2; ++i)
		{
			p = np + 1, np = l.find(' ', np + 1);
			if (np == -1)
				np = l.size();
			uvt[i] = stoi(l.substr(p, np - p));
		}
		p = np + 1, np = l.find(' ', np + 1);
		if (np == -1)
			np = l.size();
		uvt3 = stoll(l.substr(p, np - p));
//		uvt3=uvt3/(24*3600);
		if (uvt[0] == uvt[1])
			continue;
		v.push_back(uvt[0]), v.push_back(uvt[1]);
		tmin = min(tmin, uvt3);
		arcs[arcn++] = {uvt[0], uvt[1], uvt3};
		//		if(arcn%10000000==0)cout<<arcn<<endl;
	}
	sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
	vern = v.size();
	//	cout<<vern<<endl;
	auto get = [&](int k)
	{
		return (lower_bound(v.begin(), v.end(), k) - v.begin()) + 1;
	};

	for (int i = 0; i < arcn; ++i)
	{
		arcs[i].t -= tmin;
		arcs[i].src = get(arcs[i].src), arcs[i].dst = get(arcs[i].dst);
		if (arcs[i].src > arcs[i].dst)
			swap(arcs[i].src, arcs[i].dst);
		tmax = max(tmax, arcs[i].t);
		//		cout<<tmax<<endl;
		int u = arcs[i].src;
		int v = arcs[i].dst;
		if (u > v)
		{
			u = u ^ v;
			v = u ^ v;
			u = u ^ v;
		}

		if (!(find(tss[make_pair(u, v)].begin(), tss[make_pair(u, v)].end(), arcs[i].t) != tss[make_pair(u, v)].end()))
		{

			tss[make_pair(u, v)].push_back(arcs[i].t);
		}
	}
}

void tsort()
{
	for (auto &x : tss)
	{
		sort(x.second.begin(), x.second.end());
	}
}

unordered_map<pair<int, int>, vector<pair<int, double>>, pair_hash> TF;
int isLCP(int ida, int idb, int idc, long long a, long long b, long long c)
{
	return (idb - ida) * (c - a) - (idc - ida) * (b - a);
}
int isT(int ida, int idb, int idc, long long a, long long b, long long c)
{
	return (idc - idb + 1) * (c - a + 1) - (idc - ida + 1) * (c - b + 1);
}
int dq[VMAX];
void buildTF_baseline()
{
	for (int u = 1; u <= vern; ++u)
	{
//		cout<<u<<endl;
		for (int v : G2[u])
		{
			if (u > v)
				continue;
			vector<long long> tuv(tss[make_pair(u, v)]);
			int tsz = tuv.size();
			int t = 1;
			while (1)
			{
				double ans = 0x7f7f7f7f * 1.0;
				double eps = 1e-6;
				int length = 0;
				for (int i = 0; i < tsz; ++i)
				{
					for (int j = i + t - 1; j < tsz; ++j)
					{
						double res = (tuv[j] - tuv[i] + 1) * 1.0 / (j - i + 1) * 1.0;
						if (res + eps < ans)
						{
							ans = res;
							length = j - i + 1;
						}
						else if (fabs(ans - res) < eps)
						{
							ans = res;
							length = max(length, j - i + 1);
						}
					}
				}
				if (length >= 1)
					TF[make_pair(u, v)].push_back(make_pair(length, ans));
				t = length + 1;
				if (t > tsz)
					break;
			}
		}
	}
}
long long ttt=0ll;
long long ttt2=0ll; 
void buildTF()
{
	for (int u = 1; u <= vern; ++u)
	{
//					if(u%100==0)cout<<u<<endl;
		for (int v : G2[u])
		{
			if (u > v)
				continue;
			vector<long long> tuv(tss[make_pair(u, v)]);

			int tsz = tuv.size();
			ttt+=tsz*1ll;
			int t = 1;
			while (1)
			{
				int l = 0;
				int r = 0;
				double ans = 0x7f7f7f7f * 1.0;
				double eps = 1e-6;
				//				double targ=1.0/qf;
				int length = 0;
				for (int i = t - 1; i < tsz; ++i)
				{

					while (r - l > 1 && isLCP(dq[r - 2], dq[r - 1], i - t + 1, tuv[dq[r - 2]], tuv[dq[r - 1]], tuv[i - t + 1]) >= 0)
					{
						r--;
					}

					dq[r++] = i - t + 1;
					while (r - l > 1 && isT(dq[l], dq[l + 1], i, tuv[dq[l]], tuv[dq[l + 1]], tuv[i]) > 0)
					{
						++l;
					}
					double res = 0.0;
					res = (tuv[i] - tuv[dq[l]] + 1) * 1.0 / ((i - dq[l] + 1) * 1.0);
					if (res + eps < ans)
					{
						ans = res;
						length = i - dq[l] + 1;
					}
					else if (fabs(ans - res) < eps)
					{
						ans = res;
						length = max(length, i - dq[l] + 1);
					}
				}
				if (length >= 1)
					TF[make_pair(u, v)].push_back(make_pair(length, ans));
				//				cout<<t<<" "<<1.0/ans<<endl;
				t = length + 1;
				ttt2++;
				if (t > tsz)
					break;
			}
		}
	}
}
int gbu = 0;
bool cmpt(int v1, int v2)
{
	int x1 = gbu;
	int x2 = v1;
	int x3 = gbu;
	int x4 = v2;
	if (x1 > x2)
	{
		x1 = x1 ^ x2;
		x2 = x1 ^ x2;
		x1 = x1 ^ x2;
	}
	if (x3 > x4)
	{
		x3 = x3 ^ x4;
		x4 = x3 ^ x4;
		x3 = x3 ^ x4;
	}
	int sz1 = tss[make_pair(x1, x2)].size();
	int sz2 = tss[make_pair(x3, x4)].size();
	return sz1 < sz2;
}
void tsort2()
{
	int id = 0;
	for (auto &u : G2)
	{
		gbu = id;
		sort(u.begin(), u.end(), cmpt);
		mxoffset[id] = u.size();
		id++;
	}
}
vector<vector<pair<int, int>>> tarcs;
void preCalculateTedge()
{
	for (int i = 0; i <= mxt + 3; ++i)
	{
		tarcs.push_back(vector<pair<int, int>>());
	}
	for (auto x : TF)
	{
		int u = x.first.first;
		int v = x.first.second;
		for (auto y : x.second)
		{
			int t = y.first;
			tarcs[t].push_back(make_pair(u, v));
		}
	}
}
void buildG2()
{
	for (int i = 0; i < arcn; ++i)
	{
		int src = arcs[i].src;
		int dst = arcs[i].dst;
		long long t = arcs[i].t;
		if (src > dst)
		{
			src = src ^ dst;
			dst = src ^ dst;
			src = src ^ dst;
		}
		if (tss[make_pair(src, dst)].size() >= 1)
		{
			if (!(find(G2[src].begin(), G2[src].end(), dst) != G2[src].end()))
			{
				G2[src].push_back(dst);
				G2[dst].push_back(src);
			}
		}
		else
		{
			tss.erase(make_pair(src, dst));
		}
	}
}

unordered_map<int, int> Mv;
unordered_map<pair<int, int>, int, pair_hash> Mc;
set<pair<int, int>> Hv;
int dmax = 0;

int *core = new int[VMAX];
int *_core = new int[VMAX];

bool cAdd(int src, int dst)
{
	if (Mc.count({src, dst}) == 0)
	{
		Mc[{src, dst}] = 1;
		return true;
	}
	Mc[{src, dst}]++;
	return false;
}

void vAdd(int v)
{
	if (Mv.count(v) == 0)
		Mv[v] = 0;
	Hv.erase({Mv[v], v});
	Mv[v]++;
	Hv.insert({Mv[v], v});
}

void initMH()
{
	Mv.clear();
	Mc.clear();
	Hv.clear();

	for (int i = 1; i <= vern; ++i)
	{
		for (int j = 0; j < G2[i].size(); ++j)
		{
			int src = i;
			int dst = G2[i][j];
			if (src > dst)
			{
				src = src ^ dst;
				dst = src ^ dst;
				src = src ^ dst;
			}
			if (cAdd(src, dst))
			{
				vAdd(src);
				vAdd(dst);
			}
		}
	}
}

void vUpd(int v)
{
	if (Mv.count(v) == 0)
		return;
	int d = Mv[v];
	Hv.erase({d, v});
	Hv.insert({d - 1, v});
	Mv[v]--;
}
void coredecomp()
{
	memset(core, 0, sizeof(int) * VMAX);
	memset(_core, 0, sizeof(int) * VMAX);
	initMH();
	if (Hv.empty())
		return;
	dmax = (--Hv.end())->first;
	for (int k = 2; k <= dmax + 1; ++k)
	{
		while (Hv.size() && (Hv.begin()->first < k))
		{
			auto nv = *(Hv.begin());
			Hv.erase(Hv.begin());
			int n = nv.first, v = nv.second;
			Mv.erase(v);
			core[v] = k - 1;
			mxk = max(mxk, core[v]);
			_core[v] = k - 1;
			for (int i = 0; i < G2[v].size(); ++i)
			{
				int u = G2[v][i];
				int src = v;
				int dst = u;
				if (src > dst)
				{
					src = src ^ dst;
					dst = src ^ dst;
					src = src ^ dst;
				}
				vUpd(u);
			}
		}
	}
}

vector<vector<unordered_map<int, int>>> CFN;
unordered_map<pair<int, int>, int, pair_hash> TFoffset;
vector<vector<vector<pair<long long, int>>>> CFI;
vector<vector<int>> CF;

int getf(int u, int v)
{
	int uu = u;
	int vv = v;
	if (uu > vv)
		swap(uu, vv);
	int id = TFoffset[make_pair(uu, vv)];
	if (id >= TF[make_pair(uu, vv)].size())
		return -feps;
	double f = 1.0 / TF[make_pair(uu, vv)][id].second;
	int fi = (int)(f * feps);
	return fi;
}
void initCFN()
{
	initMH();
	CFN.push_back(vector<unordered_map<int, int>>());
	for (int i = 1; i <= vern; ++i)
	{
		CFN.push_back(vector<unordered_map<int, int>>());
		CFN[i].push_back(unordered_map<int, int>());
		for (int j = 1; j <= _core[i]; ++j)
		{
			CFN[i].push_back(unordered_map<int, int>());
		}
	}
	for (int i = 1; i <= vern; ++i)
		for (int j = 1; j <= _core[i]; ++j)
			CFN[i][j].clear();
	for (auto vn : Mv)
	{
		int v = vn.first;
		for (int kk = 1; kk <= core[v]; kk++)
		{
			for (int i = 0; i < G2[v].size(); ++i)
			{
				int u = G2[v][i];
				if (core[u] < kk)
					continue;
				int qd = getf(v, u);
				if (CF[u][kk] >= CF[v][kk] && qd >= CF[v][kk])
				{
					CFN[v][kk][u] = 1;
				}
			}
		}
	}
}
void updCFN(int u, int k)
{
	CFN[u][k].clear();
	for (int i = offset[u]; i < mxoffset[u]; ++i)
	{

		int v = G2[u][i];
		if (core[v] < k)
			continue;
		if (CF[v][k] >= CF[u][k] && CF[u][k] <= getf(u, v))
			CFN[u][k][v] = 1;
	}
}

void clear(queue<pair<int, int>> &q)
{
	queue<pair<int, int>> empty;
	swap(empty, q);
}
void initCF()
{
	CF.push_back(vector<int>());
	for (int v = 1; v <= vern; ++v)
	{
		CF.push_back(vector<int>());
		CF[v].push_back(0);
		for (int k = 1; k <= core[v]; ++k)
		{
			if (!(CFI[v][k].empty()))
			{
				CF[v].push_back(CFI[v][k][0].second);
			}
		}
	}
}

queue<pair<int, int>> pq;
unordered_map<pair<int, int>, int, pair_hash> XPQ;
void check(int u, int v, int k, int qd1, int qd2)
{
	if (CF[v][k] >= CF[u][k] && qd1 >= CF[u][k] && qd2 < CF[u][k])
	{
		CFN[u][k].erase(v);
	}

	if (CFN[u][k].size() < k)
	{
		if (!XPQ[make_pair(u, k)])
		{
			pq.push(make_pair(u, k));
			XPQ[make_pair(u, k)] = 1;
		}
	}
}
int findk(vector<int> &a, int s, int t, int k)
{
	int i = s, j = t;
	int temp = a[s];
	if (s < t)
	{
		while (i != j)
		{
			while (j > i && a[j] >= temp)
				j--;
			a[i] = a[j];
			while (i < j && a[i] <= temp)
				i++;
			a[j] = a[i];
		}
		a[i] = temp;
		if (k - 1 == i)
			return a[i];
		else if (k - 1 < i)
		{
			return findk(a, s, i - 1, k);
		}
		else
			return findk(a, i + 1, t, k);
	}

	else
		return a[k - 1];
}
int localCF(int u, int t, int k)
{
	vector<int> T;
	int ofu = offset[u];
	while (ofu < mxoffset[u])
	{
		int v = G2[u][ofu];
		int uu = u;
		if (uu > v)
			swap(uu, v);
		int sz = tss[make_pair(uu, v)].size();
		if (sz < t)
			ofu++;
		else
			break;
	}
	offset[u] = ofu;
	for (int i = offset[u]; i < mxoffset[u]; ++i)
	{
		int v = G2[u][i];
		if (core[v] < k)
			continue;
		int x = getf(u, v);
		if (x <= 0)
			continue;
		int ct = min(CF[v][k], x);
		T.push_back(-ct);
	}
	if (T.size() >= k)
	{
		return -findk(T, 0, T.size() - 1, k);
	}
	core[u] = min(core[u], k - 1);
	return -feps;
}

void buildCFI()
{
	coredecomp();
		cout<<mxk<<endl;
	CFI.push_back(vector<vector<pair<long long, int>>>());
	for (int i = 1; i <= vern; ++i)
	{
		CFI.push_back(vector<vector<pair<long long, int>>>());
		CFI[i].push_back(vector<pair<long long, int>>());
		for (int k = 1; k <= core[i]; ++k)
		{
			CFI[i].push_back(vector<pair<long long, int>>());
			CFI[i][k].push_back({1ll, 1 * feps});
		}
	}
	initCF();
	initCFN();
	memset(offset, 0, sizeof(offset));
	tsort2();
	for (int t = 2; t <= mxt + 1; ++t)
	{	
//		cout<<t<<endl;
		int hh, tt;
		hh = 0, tt = -1;
		XPQ.clear();
		clear(pq);
		for (auto p : tarcs[t - 1])
		{
			int u = p.first;
			int v = p.second;
			int uu = u;
			int vv = v;
			if (uu > vv)
				swap(uu, vv);
			int qd1 = getf(uu, vv);
			TFoffset[make_pair(uu, vv)]++;
			int qd2 = getf(uu, vv);
			for (int k = 1; k <= min(core[uu], core[vv]); ++k)
			{
				check(u, v, k, qd1, qd2);
				check(v, u, k, qd1, qd2);
			}
		}
		while (pq.size())
		{
			auto uk = pq.front();
			pq.pop();
			XPQ.erase(uk);
			int u = uk.first;
			int k = uk.second;

			int oldCF = CF[u][k];
			CF[u][k] = localCF(u, t, k);

			if (CFI[u][k].back().first < t)
				CFI[u][k].push_back({t, CF[u][k]});
			else
				CFI[u][k].back().second = min(CFI[u][k].back().second, CF[u][k]);
			if (CF[u][k] < 0)
				CFN[u][k].clear();
			else
				updCFN(u, k);
			for (int i = offset[u]; i < mxoffset[u]; ++i)
			{
				int v = G2[u][i];
				if (_core[v] < k)
					continue;
				if (CF[v][k] < 0)
					continue;
				if (CFN[v][k].count(u) && CF[u][k] < CF[v][k])
				{
					CFN[v][k].erase(u);
					if (CFN[v][k].size() < k)
					{
						if (!XPQ[make_pair(v, k)])
						{
							pq.push(make_pair(v, k));
							XPQ[make_pair(v, k)] = 1;
						}
					}
				}
			}
		}
	}
}

int getCT(int v, int k, int t)
{
	if (CFI.size() < v + 1 || CFI[v].size() < k + 1 || CFI[v][k].empty())
		return -feps;
	auto it = upper_bound(CFI[v][k].begin(), CFI[v][k].end(), t,
						  [](int value, const pair<long long, int> &element)
						  {
							  return value < element.first;
						  });
	return (--it)->second;
}
int *isRes = new int[VMAX];
// 计算基本数据类型的占用空间
template <typename T>
size_t calculate_vector_size(const std::vector<T> &vec)
{
	size_t size = vec.size() * sizeof(T); // 只计算元素占用的空间
	return size;
}

// 特化模板，用于处理 vector 容器类型的元素
template <typename T>
size_t calculate_vector_size(const std::vector<std::vector<T>> &vec)
{
	size_t size = sizeof(int); // 不计算vec本身的指针占用，只计算内容
	for (const auto &item : vec)
	{
		size += calculate_vector_size(item); // 递归计算内部 vector 的大小
	}
	return size;
}

// 继续特化处理三维vector
template <typename T>
size_t calculate_vector_size(const std::vector<std::vector<std::vector<T>>> &vec)
{
	size_t size = sizeof(int); // 同样，只计算内容
	for (const auto &item : vec)
	{
		size += calculate_vector_size(item); // 递归计算
	}
	return size;
}
int main(int argc, char *argv[])
{
	loadgraph((string("../data/") + string(argv[1])).c_str());
	choice = stoi(argv[2]);
	buildG2();
	tsort();
	auto start = system_clock::now();
	if (choice == 0)
	{
		buildTF();
		cout<<arcn<<endl;
		cout<<ttt<<" "<<ttt2<<endl; 
	}
	else
		buildTF_baseline();
	auto end2 = system_clock::now();
		auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start);
	cout << fixed << setprecision(6) << static_cast<double>(duration2.count()) * 1.0 << "ms"<<endl;
	for (int i = 1; i <= vern; ++i)
	{
		for (int j = 0; j < G2[i].size(); ++j)
		{
			if (i > G2[i][j])
				continue;
			int sz = tss[make_pair(i, G2[i][j])].size();
			mxt = max(mxt, sz);
		}
	}
	preCalculateTedge();
		cout<<mxt<<endl;
		
	buildCFI();
	
	cout << fixed << setprecision(2) << 1.0 * sizeof(arc) * arcn / (1024 * 1024) << endl;
//cout << fixed << setprecision(2) << 1.0 * ttt*3*4 / (1024 * 1024) << endl;
	cout << fixed << setprecision(2) << (1.0 * calculate_vector_size(CFI)+mxk*vern*4.0) / (1024 * 1024) << endl;
	auto end = system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	cout << fixed << setprecision(6) << static_cast<double>(duration.count()) * 1.0 << "ms";
	for (int k = 1; k <= mxk; ++k)
	{
		freopen((std::string("../result/CFI_") + to_string(k) + std::string(argv[1])).c_str(), "w", stdout);
		for (int i = 1; i <= vern; ++i)
		{
			cout << i << endl;
			if (CFI[i].size() >= k + 1)
			{
				for (auto p : CFI[i][k])
				{
					cout << p.first << " " << p.second << endl;
				}
			}
		}
	}
	freopen((std::string("../result/CFI_") + std::string(argv[1])).c_str(), "w", stdout);
	for (int i = 1; i <= vern; ++i)
	{
		cout << i << endl;
		for (int k = 1; k <= _core[i]; ++k)
		{
			if (CFI[i][k].empty())
				break;
			cout << k << endl;
			for (auto p : CFI[i][k])
			{
				cout << p.first << " " << p.second << endl;
			}
		}
		cout << endl;
	}
	return 0;
}
