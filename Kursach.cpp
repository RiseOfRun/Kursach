#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
typedef vector<vector<double>> Matrix;
class Net
{
public:
	Net()
	{
		
	}
	vector<vector<double>> Node;
	vector<vector<int>> Elements;
	double lambda;

	void BuildNet(double xmin, double xmax, double ymin, double ymax, int nx,int ny)
	{
		double hx = (xmax - xmin) / nx;
		double hy = (ymax - ymin) / ny;
		Node = vector<vector<double>>((nx + 1) * (ny + 1));
		Node[0] = vector<double>{xmin, ymin};
		double y = ymin;

		for (int i = 0; i < ny; i++)
		{
			double x = xmin;
			for (int j = 0; j < nx; j++)
			{
				Node[j + 1] = { x + hx, y};
				Node[i + nx + 1] = { x,y + hy };
				Elements.push_back({ j+i*(nx+1),j + 1+i*(nx+1), j + (nx+1)*(i+1)});
				x = x + hx;	
			}
			y = y + hy;
		}
		Node[Node.size()-1]={ xmax,ymax };

		for (int i = ny; i>0; i--)
		{
			for (int j = nx; j >0; j--)
			{
				Elements.push_back({ j + i * (nx + 1) - nx - 1,j - 1 + i * (nx + 1), j + i * (nx + 1)});
			}
		}
	}
private:
};


class Eq
{
public:
	Net FuckingNet;
	double Lambda=10, gamma;
	vector<int> ia, ja;
	vector<double> al, au,di;

	Eq(Net Net)
	{
		FuckingNet = Net;
	}

	vector<vector<double>> BuildG(vector<vector<double>>& D_1,  double DetD)
	{
		vector<vector<double>> G(3);
		double multix = Lambda * abs(DetD) / 2.;
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				double alpha1i = D_1[i][1];
				double alpha1j = D_1[j][1];
				double alpha2i = D_1[i][2];
				double alpha2j = D_1[j][2];
				double sum = alpha1i * alpha1j + alpha2i * alpha2j;
				G[i].push_back(multix*sum);
			}
		}
		return G;
	}

	vector<vector<double>> BuildM(double gamma, double DetD)
	{
		vector<vector<double>> M = {
		{2,1,1},
		{1,2,1},
		{1,1,2}
		};

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				M[i][j] = M[i][j] * gamma * abs(DetD) / 24;
			}
		}
		return M;
	}

	vector<vector<double>> BuildLocal(int i, int j, int k)
	{
		double x1 = FuckingNet.Node[i][0];
		double x2 = FuckingNet.Node[j][0];
		double x3 = FuckingNet.Node[k][0];
		double y1 = FuckingNet.Node[i][1];
		double y2 = FuckingNet.Node[j][1];
		double y3 = FuckingNet.Node[k][1];
		
		vector<vector<double>> D{
		{1,1,1},
		{x1,x2,x3},
		{y1,y2,y3}
		};

		double DetD = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

		vector<vector<double>> D_1{
		{x2 * y3 - x3 * y2, y2 - y3, x3 - x2},
		{x3 * y1 - x1 * y3, y3 - y1, x1 - x3},
		{x1 * y2 - x2 * y1, y1 - y2, x2 - x1}
		};
		D_1[0][1] = +0.;
		D_1[2][2] = +0.;
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				D_1[i][j] /= DetD;
			}
		}
		vector<vector<double>> G = BuildG(D_1, DetD);
		vector<vector<double>> M = BuildM(gamma , DetD);

		vector<vector<double>> A(3);
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				A[i].push_back(G[i][j]+M[i][j]);

			}
		}
		return A;
	}

	void BuildProfile()
	{
		vector<vector<int>> profile(FuckingNet.Node.size());


		for (size_t i = 0; i < FuckingNet.Elements.size(); i++)
		{
			for (size_t j = 1; j < 3; j++)
			{
				for (size_t k = 0; k < j; k++)
				{
					int current = FuckingNet.Elements[i][j];
					int node = FuckingNet.Elements[i][k];
					if (!count(profile[current].begin(), profile[current].end(),node))
					{
						if (profile[current].size() !=0 && profile[current][profile[current].size() - 1] > node)
						{
							for (int l = 0; l < profile[current].size(); l++)
							{
								if (node>profile[current][l])
								{
									profile[current].insert(profile[current].begin() + l, node);
									break;
								}
							}
						}
						else
						{
							profile[current].push_back(node);
						}
					}
				}
			}
		}
		
		ia.push_back(1);
		int count = 0;
		for (size_t i = 1; i < FuckingNet.Node.size(); i++)
		{
			ia.push_back(ia[i - 1] + count);
			count = 0;
			for (size_t j = 0; j < profile[i].size(); j++)
			{
				ja.push_back(profile[i][j]);
				count++;
			}
		}
		ia.push_back(ia[ia.size()-1] + count);
		al = vector<double>(ia[ia.size() - 1] - 1);
		au = vector<double>(ia[ia.size() - 1] - 1);
		di = vector<double>(FuckingNet.Node.size());
	}

	void ToGLobal(Matrix &A, vector<int> &el)
	{
		for (size_t i = 0; i < 3; i++)
		{
			di[el[i]] = di[el[i]] + A[i][i];
		}

		for (int i = 0; i < 3; i++)
		{
			int ibeg = ia[el[i]]-1;
			for (int j = 0; j < i; j++)
			{
				int iend = ia[el[i] + 1] - 1;
				while (ja[ibeg] != el[j])
				{
					int ind = (ibeg + iend) / 2 + 1;
					if (ja[ind]<=el[j])
					{
						ibeg = ind;
					}
					else
					{
						iend = ind;
					}
				}
				al[ibeg] += A[i][j];
				au[ibeg] += A[j][i];
				ibeg++;
			}
		}

	}

	vector<Matrix> BuildGlobal()
	{
		vector<Matrix> Locals;
		for (size_t i = 0; i < FuckingNet.Elements.size(); i++)
		{
			vector<int> element = FuckingNet.Elements[i];
			Locals.push_back(BuildLocal(element[0],element[1],element[2]));
			Locals[i][0][2] = 3;
			ToGLobal(Locals[i],element);
		}
		return Locals;
	}

private:

};


int main()
{
	Net testNet;
	testNet.Node.push_back({2,0});
	testNet.Node.push_back({ 2,1 });
	testNet.Node.push_back({ 3,1 });
	testNet.Node.push_back({ 2,14 });
	testNet.Node.push_back({ 7,4 });

	testNet.Elements.push_back({ 0,1,2 });
	testNet.Elements.push_back({ 2,3,4 });
	testNet.Elements.push_back({ 1,2,3 });




	//testNet.BuildNet(0, 2, 0, 2, 1, 1);

	Eq equity(testNet);
	equity.BuildProfile();
	vector<Matrix> Test = equity.BuildGlobal();

    std::cout << "Hello World!\n";
}

