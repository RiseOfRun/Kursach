#include <iostream>
#include <vector>

using namespace std;

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
		Node = vector<vector<double>>((nx+1)*(ny+1));
		Node.push_back(vector<double>{xmin, ymin});
		double y = ymin;

		for (int i = 0; i < ny; i++)
		{
			double x = xmin;
			for (int j = 0; j < nx; j++)
			{
				Node[j + 1].push_back(x+hx);
				Node[j + 1].push_back(y);
				Node[i + nx + 1].push_back(x);
				Node[i + nx + 1].push_back(y+hy);
				Elements.push_back({ j+i*(nx+1),j + 1+i*(nx+1), j + (nx+1)*(i+1)});
				x = x + hx;	
			}
			y = y + hy;
		}

		for (int i = ny; i>0; i--)
		{
			for (int j = nx; j >0; j--)
			{
				Elements.push_back({ j + i * (nx + 1),j - 1 + i * (nx + 1), j+i*(nx+1)-nx-1 });
			}
		}
	}
private:
};


class Eq
{
public:
	Net FuckingNet;
	double Lambda;
	Eq()
	{
		FuckingNet.BuildNet(0, 2, 0, 2, 2, 2);
	}

	vector<vector<double>> BuildG(vector<vector<double>>& D_1,  double DetD)
	{
		vector<vector<double>> G(3);
		double multix = Lambda * DetD / 2;
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				G[i].push_back(multix*(D_1[i][1]*D_1[j][1]- D_1[i][2] * D_1[j][2]));
			}
		}
		return G;
	}
	void BuildLocal(int i, int j, int k)
	{
		double x1 = FuckingNet.Node[i][0];
		double x2 = FuckingNet.Node[j][0];
		double x3 = FuckingNet.Node[k][0];
		double y1 = FuckingNet.Node[i][1];
		double y2 = FuckingNet.Node[j][1];
		double y3 = FuckingNet.Node[k][1];
		
		vector<vector<double>> D{
		vector<double>{1,1,1},
		vector<double> {x1,x2,x3},
		vector<double> {y1,y2,y3}
		};

		double DetD = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

		vector<vector<double>> D_1{
		vector<double> {x2 * y3 - x3 * y2, y2 - y3, x3 - x2},
		vector<double> {x3 * y1 - x1 * y3, y3 - y1, x1 - x3},
		vector<double> {x1 * y2 - x2 * y1, y1 - y2, x2 - x1}
		};

		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				D_1[i][j] /= DetD;
			}
		}
		vector<vector<double>> G = BuildG(D_1, DetD);
		
		double mesP = 1 / 2 * abs(DetD); //площадь треугольного кусочка области
		double L1, L2, L3;

	}
	~Eq();

private:

};


int main()
{
	Net testNet;
	testNet.BuildNet(0, 2, 0, 2, 2, 2);
    std::cout << "Hello World!\n";
}

