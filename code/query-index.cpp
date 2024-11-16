#include <bits/stdc++.h>
using namespace std;
using namespace chrono;
const int VMAX = 4000000;
const int EMAX = 70000000;
const int feps = 100000;
int arcn = 0;
vector<vector<pair<int, int>>> CFI;
int vern = 0;
int getCF(int v, int t)
{
	if (CFI.size() < v + 1 || CFI[v].empty())
		return -feps;
	auto it = upper_bound(CFI[v].begin(), CFI[v].end(), t,
						  [](int value, const pair<int, int> &element)
						  {
							  return value < element.first;
						  });
	return (--it)->second;
}
vector<int> res;
vector<vector<int>> G2(VMAX);
int isRes[VMAX + 3];
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
unordered_map<pair<int, int>, vector<long long>, pair_hash> tss;
unordered_map<int,int>node;
struct arc
{
	int src, dst;
	long long t;
} arcs[EMAX];
void loadgraph(const char *name)
{
	ifstream fin(name, ios::in);
	vector<int> v;
	long long tmin = 0x7fffffffll;

	string l;
	while (getline(fin, l))
	{
		//		cout<<l<<endl;
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
		node[get(arcs[i].src)]=arcs[i].src;
		node[get(arcs[i].dst)]=arcs[i].dst;
		arcs[i].src = get(arcs[i].src), arcs[i].dst = get(arcs[i].dst);
		if (arcs[i].src > arcs[i].dst)
			swap(arcs[i].src, arcs[i].dst);
		int u = arcs[i].src;
		int v = arcs[i].dst;
		if (u > v)
		{
			u = u ^ v;
			v = u ^ v;
			u = u ^ v;
		}
		if (!(find(G2[u].begin(), G2[u].end(), v) != G2[u].end()))
		{
			G2[u].push_back(v);
			G2[v].push_back(u);
		}
		if (!(find(tss[make_pair(u, v)].begin(), tss[make_pair(u, v)].end(), arcs[i].t) != tss[make_pair(u, v)].end()))
		{
			tss[make_pair(u, v)].push_back(arcs[i].t);
		}
	}
}
int main(int argc, char *argv[])
{
	loadgraph((string("../data/") + string(argv[1])).c_str());
	char *graph = argv[1];
	int qk = stoi(argv[2]);
	int qt = stoi(argv[3]);
	double qf = stod(argv[4]);
	if (freopen((string("../result/CFI_") + to_string(qk) + string(graph)).c_str(), "r", stdin) == nullptr)
		cout << "k is too large" << endl;
	string ii, num12;
	int i;
	CFI.push_back(vector<pair<int, int>>());
	while (getline(cin, ii))
	{
		CFI.push_back(vector<pair<int, int>>());
		i = stoi(ii);
		vern = max(vern, i);
		while (getline(cin, num12))
		{
			size_t pos = num12.find(' ');

			if (pos >= num12.length())
			{
				i++;
				CFI.push_back(vector<pair<int, int>>());
				vern = max(vern, i);
				continue;
			}
			string firstNumStr = num12.substr(0, pos);
			string secondNumStr = num12.substr(pos + 1);
			int num1 = stoi(firstNumStr);
			int num2 = stoi(secondNumStr);
			CFI[i].push_back({num1, num2});
			if (num2 == -feps)
				break;
		}
	}
	int targ = (int)(qf * feps);
	auto start = system_clock::now();
	for (int i = 1; i <= vern; ++i){
		if (getCF(i, qt) >= targ){
			res.push_back(i);
			isRes[i] = 1;
		}
	}

 

	auto end = system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	cout << fixed << setprecision(6) << static_cast<double>(duration.count()) * 1.0 << "μs";
//	freopen(("../result/" + string(argv[1]) + "_" + string(argv[2]) + "_" + string(argv[3]) +
//			 "_" + string(argv[4]) + ".txt")
//				.c_str(),
//			"w", stdout);

	long long bs=0;		
	int countEdges = 0;  // 用于计算边的总数
    int countVertices = 0;  // 用于计算在图中的顶点数量
    
	for (int i = 0; i < res.size(); ++i)
	{
		 ++countVertices;  // 计算在图中的顶点数量
		for (int j : G2[res[i]])
		{
			if (!isRes[j])
				continue;
			if (j <= res[i])
				continue;
			++countEdges; 
		}
	}
//	outputFile.close();
	for (int i = 0; i < res.size(); ++i)
	{
//		 ++countVertices;  // 计算在图中的顶点数量
		for (int j : G2[res[i]])
		{
			if (!isRes[j]){
				bs++;
				continue;
			}
				
		
		}
	}
	double averageDegree = 0.0;
    double density = 0.0;
	 double cr=0.0;
	 double ct=0.0;
    if (countVertices > 0) {
        averageDegree = 2.0 * countEdges / countVertices;
        density = (countVertices > 1) ? (2.0 * countEdges / (countVertices * (countVertices - 1))) : 0.0;
    }
    if(countVertices==0){
    	cr=0;
    	ct=0;
	}
    else {
    	cr=1.0-bs*1.0/(countVertices*(vern-countVertices));
    ct=1-bs*1.0/(2.0*countEdges+bs);
	}
//    freopen("CON","w",stdout); 
//	 freopen(("../result/ex6.1ad"+ string(argv[1])+"_" + string(argv[2])+".ad.txt").c_str(), "a", stdout);
//    cout << averageDegree <<"\t";
//     freopen(("../result/ex6.1d"+ string(argv[1])+"_" + string(argv[2])+".d.txt").c_str(), "a", stdout);
//    cout << density << "\t";
//    freopen(("../result/ex6.1cr"+ string(argv[1])+"_" + string(argv[2])+".cr.txt").c_str(), "a", stdout);
//    cout <<cr<< "\t";
//    freopen(("../result/ex6.1ct"+ string(argv[1])+"_" + string(argv[2])+".ct.txt").c_str(), "a", stdout);
//    cout <<ct<< "\t";
    return 0;
}
