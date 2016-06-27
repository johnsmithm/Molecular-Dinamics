#include <math.h>
#include <map>
#include <vector>
#include<fstream>
#include <string>
#include <fstream>
#include <list>


using namespace std;



typedef vector<vector<double> > Vector3D;
typedef vector<double > Vector1D;

typedef list<int>  Cell;

typedef vector<Cell > TableCells;

#define x first
#define y second
#define pb push_back
#define mp make_pair
#define forn(i,n) for(size_t i=0;i<n;++i)


class ParameterReader{
	public:
		ParameterReader(string filename){
			readParameters(filename);
		}
		template<typename Type>
		inline void GetParameter(const std::string& key, Type &value) {
			if(IsDefined(key)){
				value = parameters[key];
			}
		}
	private:
		bool readParameters(const std::string& filename){
			//cerr<<filename<<"\n";
			ifstream in(filename);
			string s;
			while(in>>s){
					in>>parameters[s];				
			}
			return true;
		}
		inline bool IsDefined(const std::string& key) const{
			if(key == "name") return true;

			return (parameters.count(key)!=0);
		}

		map<string,string> parameters;
};



