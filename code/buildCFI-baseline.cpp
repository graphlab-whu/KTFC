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

// const int VMAX = 2610000;
const int VMAX = 4000000;
const int EMAX = 70000000;
const int TMAX = 239800000;
const int feps = 100000;
int vern = 0;
int qt;
int qk;
double qf;
int arcn = 0;
long long tmax = 0ll;
int mxt = 0;
int choice = 0;
struct arc
{
	int src, dst;
	long long t;
} arcs[EMAX];
int isRes[VMAX];
int deg[VMAX];
int _deg[VMAX];
int pos[VMAX];
int bin[VMAX];
int vert[VMAX];
int core[VMAX];
int _core[VMAX];
int mxk = 0;
unordered_map<pair<int, int>, vector<long long>, pair_hash> tss;
unordered_map<pair<int, int>, int, pair_hash> last;
unordered_map<pair<int, int>, int, pair_hash> _last;
unordered_map<int, set<int>> fft;
vector<vector<int>> G2(VMAX);
vector<vector<int>> lk(VMAX);
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
		if (uvt[0] == uvt[1])
			continue;
		v.push_back(uvt[0]), v.push_back(uvt[1]);
		tmin = min(tmin, uvt3);
		arcs[arcn++] = {uvt[0], uvt[1], uvt3};
	}
	sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
	vern = v.size();

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
int isLCP(int ida, int idb, int idc, long long a, long long b, long long c)
{
	return (idb - ida) * (c - a) - (idc - ida) * (b - a);
}
int isT(int ida, int idb, int idc, long long a, long long b, long long c)
{
	return (idc - idb + 1) * (c - a + 1) - (idc - ida + 1) * (c - b + 1);
}
int dq[VMAX];
bool checkf_opt(int u, int v)
{
	vector<long long> tuv(tss[make_pair(u, v)]);
	int tsz = tuv.size();
	int l = 0;
	int r = 0;
	double ans = 0x7f7f7f7f * 1.0;
	double eps = 1e-6;
	double targ = 1.0 / qf;
	for (int i = qt - 1; i < tsz; ++i)
	{
		while (r - l > 1 && isLCP(dq[r - 2], dq[r - 1], i - qt + 1, tuv[dq[r - 2]], tuv[dq[r - 1]], tuv[i - qt + 1]) >= 0)
		{
			r--;
		}
		dq[r++] = i - qt + 1;
		while (r - l > 1 && isT(dq[l], dq[l + 1], i, tuv[dq[l]], tuv[dq[l + 1]], tuv[i]) > 0)
		{
			++l;
		}
		double res = 1.0 * (tuv[i] - tuv[dq[l]] + 1) / (i - dq[l] + 1);
		if (res + eps < targ || fabs(res - targ) < eps)
		{
			return 1;
		}
		ans = min(ans, res);
	}
	return ans + eps < targ || fabs(ans - targ) < eps;
}
bool checkf_baseline(int u, int v)
{
	vector<long long> tuv(tss[make_pair(u, v)]);
	int tsz = tuv.size();
	double targ = 1.0 / qf;
	double ans = 0.0;
	double eps = 1e-6;
	for (int i = 0; i < tsz; ++i)
	{
		for (int j = i + qt - 1; j < tsz; ++j)
		{
			ans = max(ans, (j - i + 1) * 1.0 / (tuv[j] - tuv[i] + 1) * 1.0);
		}
	}
	if (1.0 / ans + eps < targ || fabs(1.0 / ans - targ) < eps)
		return 1;
	else
		return 0;
}
void checkedge()
{
	double eps = 1e-6;
	for (int u = 1; u <= vern; ++u)
	{
		if (!isRes[u])
			continue;
		for (auto v : G2[u])
		{
			if (u > v)
				continue;
			if (!isRes[v])
				continue;
			vector<long long> tuv(tss[make_pair(u, v)]);
			if (choice == 0)
			{
				if (checkf_opt(u, v))
				{
					lk[u].push_back(v);
					deg[u]++;
					lk[v].push_back(u);
					deg[v]++;
				}
			}
			else
			{
				if (checkf_baseline(u, v))
				{
					lk[u].push_back(v);
					lk[v].push_back(u);
					deg[u]++;
					deg[v]++;
				}
			}
		}
	}
}
void buildTF_baseline()
{
	for (int u = 1; u <= vern; ++u)
	{
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
				{

					for (int ll = t; ll <= length; ++ll)
					{
						fft[ll].insert((int)(1.0 / ans * feps));
					}
				}
				t = length + 1;
				if (t > tsz)
					break;
			}
		}
	}
}
void buildTF()
{
	for (int u = 1; u <= vern; ++u)
	{
		for (int v : G2[u])
		{
			if (u > v)
				continue;
			vector<long long> tuv(tss[make_pair(u, v)]);

			int tsz = tuv.size();
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
				{

					for (int ll = t; ll <= length; ++ll)
					{
						fft[ll].insert((int)(1.0 / ans * feps));
					}
				}
				t = length + 1;
				if (t > tsz)
					break;
			}
		}
	}
}
void coredecomp(vector<vector<int>> nb, int ccore[])
{
	int md = -1;

	for (int i = 1; i <= vern; ++i)
	{
		md = max(md, deg[i]);
	}

	//    int bin[vern + 1];
	for (int i = 0; i <= vern; ++i)
	{
		bin[i] = 0;
	}
	for (int i = 1; i <= vern; ++i)
	{
		++bin[deg[i]];
	}
	int start = 0;
	for (int i = 0; i <= md; ++i)
	{
		int num = bin[i];
		bin[i] = start;
		start += num;
	}

	for (int i = 1; i <= vern; ++i)
	{
		int v = i;
		pos[v] = bin[deg[v]];
		vert[pos[v]] = v;
		++bin[deg[v]];
	}

	for (int i = md; i > 0; --i)
	{
		bin[i] = bin[i - 1];
	}
	bin[0] = 0;
	for (int i = 0; i < vern; ++i)
	{
		int v = vert[i];
		ccore[v] = deg[v];
		mxk = max(mxk, ccore[v]);
		for (auto u : nb[v])
		{
			if (deg[u] == deg[v])
				continue;

			int w = vert[bin[deg[u]]];
			if (u != w)
			{
				swap(vert[pos[u]], vert[pos[w]]);
				swap(pos[u], pos[w]);
			}
			++bin[deg[u]];
			--deg[u];
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
				deg[src]++;
				deg[dst]++;
			}
		}
		else
		{
			tss.erase(make_pair(src, dst));
		}
	}
}
vector<vector<vector<pair<long long, int>>>> CFI;
void buildCFI_baseline()
{
	coredecomp(G2, core);
	fft[mxt + 1].insert(feps * 2);
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
	for (int k = 1; k <= mxk; ++k)
	{
		cout << "k" << " " << k << endl;
		for (int i = 1; i <= mxt + 1; ++i)
		{
			//		cout<<i<<endl;
			int jcnt = 0;
			for (auto j : fft[i])
			{
				//			cout<<i<<" "<<j<<endl;
				jcnt++;
				qt = i;
				qf = 1.0 * j / feps;
				for (int v = 1; v <= vern; ++v)
				{
					deg[v] = 0;
				}
				for (int v = 1; v <= vern; ++v)
				{
					lk[v].clear();
				}
				for (int v = 1; v <= vern; ++v)
				{
					isRes[v] = 1;
				}
				qk = k;

				for (int v = 1; v <= vern; ++v)
				{
					if (core[v] < qk)
					{
						isRes[v] = 0;
					}
				}
				checkedge();
				coredecomp(lk, _core);
				for (int v = 1; v <= vern; ++v)
				{
					if (_core[v] < qk)
					{
						isRes[v] = 0;
					}
				}
				for (int v = 1; v <= vern; ++v)
				{
					if (isRes[v] == 0 && qk <= core[v] && jcnt == 1 && _last[make_pair(v, k)] > 0)
					{

						CFI[v][k].push_back({i, -feps});
						last[make_pair(v, k)] = -1;
					}
					if (isRes[v] == 1)
					{
						if (_last[make_pair(v, k)] != j)
						{

							if (CFI[v][k].back().first == i)
								CFI[v][k].back().second = j;
							else
								CFI[v][k].push_back({i, j});

							last[make_pair(v, k)] = j;
						}
						else
						{

							if (CFI[v][k].back().first == i)
								CFI[v][k].pop_back();
							last[make_pair(v, k)] = j;
						}
					}
				}
			}
			for (int v = 1; v <= vern; ++v)
			{
				for (int k = 1; k <= core[v]; ++k)
					_last[make_pair(v, k)] = last[make_pair(v, k)];
			}
		}
	}
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
	}
	else
		buildTF_baseline();
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
	//	cout<<mxt<<endl;
	buildCFI_baseline();
	auto end = system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	cout << fixed << setprecision(6) << static_cast<double>(duration.count()) * 1.0 << "ms";
	freopen((std::string("../result/CFI_baseline_") + std::string(argv[1])).c_str(), "w", stdout);
	for (int i = 1; i <= vern; ++i)
	{
		cout << i << endl;
		for (int k = 1; k <= core[i]; ++k)
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
