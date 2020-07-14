#include <fstream> //too slow, use cstdio instead
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <vector>
#include <list>
#include <unordered_map>
#include <string>
#include <mutex>
#include <memory>
#include <omp.h>

using namespace std;
namespace fs = std::filesystem;

template <typename T>
struct Vec3 {
    T coord[3];
    int index;
    Vec3 (T x, T y, T z, int i) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
        index = i; 
    }

    Vec3 () = default;

    template <typename T2>
    float distanceSquare(const Vec3<T2> other ) const {
        float dx = static_cast<float>(coord[0]) - other.coord[0];
        float dy = static_cast<float>(coord[1]) - other.coord[1];
        float dz = static_cast<float>(coord[2]) - other.coord[2];
        return dx*dx + dy*dy + dz*dz;
    }
};

typedef Vec3<float> Point;
typedef Vec3<int> Index3d;

template<typename T>
ostream & operator<< (ostream & out, const Vec3<T> & p) {
    return out<< "{" << p.coord[0] << " " << p.coord[1] << " " << p.coord[2] << " " << p.index << "}";
}

struct OctNode {
    std::mutex objmutex;
    vector<Point> points;
    OctNode()  {
        points.reserve(50);
    }
};

class OctGrid {
public:
    OctGrid(const float lower_bound[], const float upper_bound[], const int parts_per_dimention = 50) {
        for (int i=0; i<3; i++) {
            float d = (upper_bound[i] - lower_bound[i]) * 0.02f;
            upperBound[i] = upper_bound[i] + d;
            lowerBound[i] = lower_bound[i] - d;
            parts_per_dim[i] = parts_per_dimention;
            dim[i] = (upperBound[i] - lowerBound[i])/parts_per_dim[i];
        }
        octs = make_unique<OctNode[]>(parts_per_dim[0] * parts_per_dim[1] * parts_per_dim[2]);
    }

    OctGrid(const float lower_bound[], const float upper_bound[], const int parts_per_dimention[]) {
        for (int i=0; i<3; i++) {
            float d = (upper_bound[i] - lower_bound[i]) * 0.02f;
            upperBound[i] = upper_bound[i] + d;
            lowerBound[i] = lower_bound[i] - d;
            parts_per_dim[i] = parts_per_dimention[i];
            dim[i] = (upperBound[i] - lowerBound[i])/parts_per_dim[i];
        }
        octs = make_unique<OctNode[]>(parts_per_dim[0] * parts_per_dim[1] * parts_per_dim[2]);
    }

    inline int getOctIndex(int i, int j, int k) {
        return i*parts_per_dim[1]*parts_per_dim[2] + j*parts_per_dim[2] +k;
    }

    inline int decomposeOctIndex(int index, int dim) {
        switch (dim)
        {
        case 0:
            return index / (parts_per_dim[1] * parts_per_dim[2]);
            break;
        case 1:
            return (index / parts_per_dim[2]) % parts_per_dim[1];
            break;
        case 2:
            return index % parts_per_dim[2];
            break;
        default:
            return -1; // error
        }
    }

    inline int getOctIndex(const float coord[]) {
        int dx = static_cast<int>((coord[0] - lowerBound[0])/dim[0]);
        int dy = static_cast<int>((coord[1] - lowerBound[1])/dim[1]);
        int dz = static_cast<int>((coord[2] - lowerBound[2])/dim[2]);
        return getOctIndex(dx, dy, dz);
    }

    void insert(const Point & p) {
        int targetOct = getOctIndex(p.coord);
        std::lock_guard lock(octs[targetOct].objmutex);
        octs[targetOct].points.push_back(p);
    }

    inline bool bounded(int i, int j, int k) {
        return i>=0 && j>=0 && k>=0 && i<parts_per_dim[0] && j<parts_per_dim[1] && k<parts_per_dim[2];
    }

    Point find(const float x, const float y, const float z) {
        return find(Point(x,y,z,0));
    }

    Point find(const Point & p) {

        //bound the point p, though the nearest point found may be not so accurate if p is outside the octgrid.
        //correct way is to project point p to the 6 facets of the grid and find the nearest octblock.
        Point fixedP = p;
        for (int i=0; i<3; i++) {
            if (fixedP.coord[i]>upperBound[i]) fixedP.coord[i] = upperBound[i]- 0.02 * (upperBound[i] - lowerBound[i]);
            if (fixedP.coord[i]<lowerBound[i]) fixedP.coord[i] = lowerBound[i]+ 0.02 * (upperBound[i] - lowerBound[i]);
        }

        list<int> searchqueue;
        unordered_map<int, bool> visited;
        int currentOct = getOctIndex(fixedP.coord);
        int curI = decomposeOctIndex(currentOct, 0);
        int curJ = decomposeOctIndex(currentOct, 1);
        int curK = decomposeOctIndex(currentOct, 2);
        for (int i=-1; i<=1; i++) {
            for (int j=-1; j<=1; j++) {
                for (int k=-1; k<=1; k++) {
                    if (!bounded(curI+i, curJ+j, curK+k)) continue;
                    int curOct = getOctIndex(curI+i, curJ+j, curK+k);
                    searchqueue.push_back(curOct);
                    visited[curOct] = true;
                }
            }
        }

        float minDist = 1e38f;
        Point minPoint(0,0,0,-1);
        int searchSize = 0;
        int searchSizeLowerBound = searchqueue.size();
        int searchSizeUpperBound = 2000000; // This value will only be achieved when none of the octnodes visited has points inside it. 2000000 should run within 0.5 sec. 
        while (!searchqueue.empty() &&
            (searchSize < searchSizeLowerBound || 
            (minPoint.index==-1 && searchSize<searchSizeUpperBound)
            )) 
        {
            searchSize++;
            int currentOct = searchqueue.front();
            searchqueue.pop_front();
            for (auto it: octs[currentOct].points) {
                float d = p.distanceSquare(it);
                if (d<minDist) {
                    minDist = d;
                    minPoint = it;
                }
            }

            int curI = decomposeOctIndex(currentOct, 0);
            int curJ = decomposeOctIndex(currentOct, 1);
            int curK = decomposeOctIndex(currentOct, 2);
            for (int i=-1; i<=1; i++) {
                for (int j=-1; j<=1; j++) {
                    for (int k=-1; k<=1; k++) {
                        if (!bounded(curI+i, curJ+j, curK+k)) continue;
                        int curOct = getOctIndex(curI+i, curJ+j, curK+k);
                        if (visited.count(curOct)!=0) continue; 
                        searchqueue.push_back(curOct);
                        visited[curOct] = true;
                    }
                }
            }
        }
        
        return minPoint;
    }

private:
    float lowerBound[3];
    float upperBound[3];
    float dim[3];
    int parts_per_dim[3];
    unique_ptr<OctNode[]> octs;
};

int main(int argc, char *argv[]) {
    fs::path p1 = argv[0];
    p1.remove_filename();
    p1 /= fs::path("data");
    
    std::string dataFilePath;
    for (const auto & p: fs::directory_iterator(p1)) {
        
        dataFilePath = p.path().string();
        break;
    }
    std::cout<<dataFilePath<<std::endl;
    ifstream dataFile(dataFilePath);
    
    string dummy;
    getline(dataFile, dummy);

    vector<Point> allpoints;
    allpoints.reserve(3000000);
    while (dataFile) {
        float x,y,z;
        int ind;
        dataFile>>ind>>x>>y>>z;
        allpoints.emplace_back(x,y,z,ind);
    }
    

    int num_threads = 1;
    
    #pragma omp parallel 
    {
        #pragma omp single 
        {
            num_threads = omp_get_num_threads();
        }
    }

    const int total_num = allpoints.size();
    const int partition_num = (total_num + num_threads - 1) / num_threads;
    vector<Point> upb(num_threads, Point(-1e-38f, -1e-38f, -1e-38f, 0));
    vector<Point> lob(num_threads, Point(+3e+38f, +3e+38f, +3e+38f, 0));
    Point upper(-1e-38f, -1e-38f, -1e-38f, 0);
    Point lower(+3e+38f, +3e+38f, +3e+38f, 0);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i<num_threads; i++) {
        for (int j = i * partition_num; j < total_num && j < (i+1) * partition_num; j++) {
            for (int k = 0; k<3; k++) {
                if (allpoints[j].coord[k]<lob[i].coord[k]) {
                    lob[i].coord[k] = allpoints[j].coord[k];
                }
                if (allpoints[j].coord[k]>upb[i].coord[k]) {
                    upb[i].coord[k] = allpoints[j].coord[k];
                }
            }
        }
        #pragma omp critical
        {
            for (int j=0; j<3; j++) {
                if (lob[i].coord[j]<lower.coord[j]) {
                    lower.coord[j] = lob[i].coord[j];
                }
                if (upb[i].coord[j]>upper.coord[j]) {
                    upper.coord[j] = upb[i].coord[j];
                }
            }
        }
    }

    cout<<"num_thread: "<<num_threads<<endl;
    cout<<"total num: "<<total_num<<endl;
    cout<<"upper bound: "<<upper<<endl;
    cout<<"lower bound: "<<lower<<endl;
    cout<<"first point:" <<allpoints[0]<<endl;

    OctGrid octgrid(lower.coord, upper.coord);

    cout<<"Building Data Structure"<<endl;
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i<total_num; i++) {
        octgrid.insert(allpoints[i]);
    }

    cout<<"Building complete."<<endl
        <<"Enter coordinate query in the form 'x y z' (do not include quotes), type exit to exit"<<endl;

    string stdinLine;
    std::getline(cin, stdinLine);
    while (stdinLine.substr(0,4)!="exit") {
        stringstream ss(stdinLine);
        Point t;
        ss >> t.coord[0] >> t.coord[1] >> t.coord[2];
        if (ss.fail()) {
            cout << "bad format" << endl;
        } else {
            cout << "nearest point" << octgrid.find(t) << endl;
        }
        std::getline(cin, stdinLine);
    }

    return 0;
}