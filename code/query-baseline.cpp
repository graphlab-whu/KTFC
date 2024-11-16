#include <bits/stdc++.h>
using namespace std;
using namespace chrono;
const int VMAX = 4000000;
const int EMAX = 70000000;
int qt;
int qk;
double qf;
double feps = 100000;
int choice;
vector<vector<pair<int, long long>>> G(VMAX);
vector<vector<int>> G2(VMAX);
vector<vector<int>> lk(VMAX);
vector<int> vs;
int isRes[VMAX];
int deg[VMAX];
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
unordered_map<pair<int, int>, vector<long long>, pair_hash> mp;
bool operator!=(const pair<int, int> &p1, const pair<int, int> &p2)
{
	return (p1.first != p2.first) || (p1.second != p2.second);
}
bool operator==(const pair<int, int> &p1, const pair<int, int> &p2)
{
	return (p1.first == p2.first) && (p1.second == p2.second);
}
bool operator<(const pair<long long, int> &p1, const pair<long long, int> &p2)
{
	if (p1.first < p2.first)
		return true;
	else if (p1.first == p2.first)
	{
		return p1.second < p2.second;
	}
	else
		return false;
}
struct arc
{
	int src, dst;
	long long t;
	arc()
	{
	}
	arc(int u, int v, long long t)
	{
		this->src = u;
		this->dst = v;
		this->t = t;
	}
} arcs[EMAX];

int ide = 0;
int vn = 0;
int core[VMAX];
int pos[VMAX];
int bin[VMAX];
int vert[VMAX];
// int cntcheckf;
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
	vector<long long> tuv(mp[make_pair(u, v)]);
	//	sort(tuv.begin(),tuv.end());
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
	vector<long long> tuv(mp[make_pair(u, v)]);
	int tsz = tuv.size();
	double targ = 1.0 / qf;
	double ans = 0x7f7f7f7f * 1.0;
	double eps = 1e-6;
	for (int i = 0; i < tsz; ++i)
	{
		for (int j = i + qt - 1; j < tsz; ++j)
		{
			ans = min(ans, 1.0 * (tuv[j] - tuv[i] + 1) / (j - i + 1));
		}
	}
	return ans + eps < targ || fabs(ans - targ) < eps;
}
void addEdge(int u, int v, long long t)
{
	if (!(find(G[u].begin(), G[u].end(), make_pair(v, t)) != G[u].end()))
	{
		G[u].push_back(make_pair(v, t));
		G[v].push_back(make_pair(u, t));
		mp[make_pair(u, v)].push_back(t);
	}
	if (!(find(G2[u].begin(), G2[u].end(), v) != G2[u].end()))
	{
		G2[u].push_back(v);
		G2[v].push_back(u);
		deg[u]++;
		deg[v]++;
	}
}
void checkedge()
{
	double eps = 1e-6;
	for (int u = 0; u < vn; ++u)
	{
		if (!isRes[u])
			continue;
		for (auto v : G2[u])
		{
			if (u > v)
				continue;
			if (!isRes[v])
				continue;
			vector<long long> tuv(mp[make_pair(u, v)]);
			sort(tuv.begin(), tuv.end());
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
					deg[u]++;
					lk[v].push_back(u);
					deg[v]++;
				}
			}
		}
	}
}
void loadgraph(const char *filename)
{
	ifstream fin(filename, ios::in);

	long long tmin = 0x7fffffffll;
	string l;
	while (getline(fin, l))
	{
		int uvt[3] = {0};
		long long uvt3=0ll;
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
		vs.push_back(uvt[0]), vs.push_back(uvt[1]);
		tmin = min(tmin, uvt3);
		arcs[++ide] = {uvt[0], uvt[1], uvt[2]};
	}
	sort(vs.begin(), vs.end());
	vs.erase(unique(vs.begin(), vs.end()), vs.end());
	vn = vs.size();

	auto get = [&](int k)
	{
		return (lower_bound(vs.begin(), vs.end(), k) - vs.begin());
	};
	for (int i = 1; i <= ide; ++i)
	{
		arcs[i].t -= tmin;
		int u = get(arcs[i].dst);
		int v = get(arcs[i].src);
		//		addEdge(u,v,arcs[i].t);
		if (u > v)
		{
			u = u ^ v;
			v = u ^ v;
			u = u ^ v;
		}
		addEdge(u, v, arcs[i].t);
	}
}
void coredecomp(vector<vector<int>> nb)
{
	int md = -1;

	for (int i = 0; i < vn; ++i)
	{
		md = max(md, deg[i]);
	}

	//    int bin[vn + 1];
	for (int i = 0; i <= vn; ++i)
	{
		bin[i] = 0;
	}
	for (int i = 0; i < vn; ++i)
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

	for (int i = 0; i < vn; ++i)
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
	for (int i = 0; i < vn; ++i)
	{
		int v = vert[i];
		core[v] = deg[v];
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
	for (int i = 0; i < vn; ++i)
	{
		if (core[i] < qk)
		{
			isRes[i] = 0;
		}
	}
}
int main(int argc, char *argv[])
{
	loadgraph((string("../data/") + string(argv[1])).c_str());
	qk = stoi(argv[2]);
	qt = stoi(argv[3]);
	qf = stod(argv[4]);
	choice = stoi(argv[5]);
	for (int i = 0; i < vn; ++i)
	{
		isRes[i] = 1;
	}
	auto start = system_clock::now();
	coredecomp(G2);
	for (int i = 0; i <= vn; ++i)
		deg[i] = 0;
	checkedge();

	coredecomp(lk);
	auto end = system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	cout << fixed << setprecision(6) << static_cast<double>(duration.count()) * 1.0 << "¦Ìs" << endl;

	return 0;
}
