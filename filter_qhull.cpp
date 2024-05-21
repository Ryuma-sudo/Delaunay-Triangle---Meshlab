#include "filter_qhull.h"
#include "qhull_tools.h"
#include <vcg/complex/algorithms/convex_hull.h>

using namespace std;
using namespace vcg;

QhullPlugin::QhullPlugin()
{
	typeList = {FP_DELAUNAY};

	for (ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QhullPlugin::~QhullPlugin()
{
}

QString QhullPlugin::pluginName() const
{
	return "FilterTest2";
}

QString QhullPlugin::filterName(ActionIDType f) const
{
	switch (f) {
	case FP_DELAUNAY: return QString("Delaunay Triangulation 2");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::pythonFilterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_DELAUNAY: return QString("Convert cloud point to delaunay triangulation");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_DELAUNAY: return QString("Convert cloud point to delaunay triangulation");
	default: assert(0);
	}
	return QString("Error: Unknown Filter");
}

QhullPlugin::FilterClass QhullPlugin::getClass(const QAction* a) const
{
	switch (ID(a)) {
	case FP_DELAUNAY: return FilterClass(FilterPlugin::Remeshing);
	default: assert(0);
	}
	return FilterClass(0);
}

RichParameterList QhullPlugin::initParameterList(const QAction* action, const MeshModel& m)
{
	RichParameterList parlst;
	switch (ID(action)) {
	case FP_DELAUNAY:
	default: break; // do not add any parameter for the other filters
	}
	return parlst;
}
// Triangulation start here

coordT* readpointsFromMesh(int* numpoints, int* dimension, MeshModel& m)
{
	coordT *points, *coords;

	coords = points = (coordT*) malloc((*numpoints) * (*dimension) * sizeof(coordT));

	int                    cnt = 0;
	CMeshO::VertexIterator vi;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
		if (!(*vi).IsD()) {
			for (int ii = 0; ii < *dimension; ++ii)
				*(coords++) = (*vi).P()[ii];
			++cnt;
		}
	assert(cnt == m.cm.vn);

	return (points);
}

// check if point P3 is on the left of (or lies on) the line P1 --> P2 or not

bool isOnLeft(int P1, int P2, int P3, coordT** List)
{
	return (
		(List[0][P2] - List[0][P1]) * (List[1][P3] - List[1][P1]) -
			(List[1][P2] - List[1][P1]) * (List[0][P3] - List[0][P1]) >= 0
		);
}

// locate the triangle enclosing the given point

int triLocate(int** Vertex, int** Edge, int P, coordT** List, int numTri)
{
	int triLoc = numTri;
	while (true) {
		int i = 0;
		for (i = 0; i < 3; i++) {
			int v1 = Vertex[i][triLoc];
			int v2 = Vertex[(i + 1) % 3][triLoc];
			if (!isOnLeft(v1, v2, P, List)) {
				triLoc = Edge[(i + 2) % 3][triLoc];
				break;
			}
		}
		if (i == 3)
			break;
	}
	return triLoc;
}

// check if point P is inside the circumcircle of triangle T

bool isInCircumcircle(coordT** List, int P, int P1, int P2, int P3)
{
	coordT x13 = List[0][P1] - List[0][P3];
	coordT x23 = List[0][P2] - List[0][P3];
	coordT x1P = List[0][P1] - List[0][P];
	coordT x2P = List[0][P2] - List[0][P];
	coordT y13 = List[1][P1] - List[1][P3];
	coordT y23 = List[1][P2] - List[1][P3];
	coordT y1P = List[1][P1] - List[1][P];
	coordT y2P = List[1][P2] - List[1][P];

	coordT cosA = x13 * x23 + y13 * y23;
	coordT cosB = x2P * x1P + y2P * y1P;

	if (cosA >= 0 && cosB >= 0) {
		return false;
	}
	else if (cosA < 0 && cosB < 0) {
		return true;
	}
	else {
		coordT sinA  = x13 * y23 - x23 * y13;
		coordT sinB  = x2P * y1P - x1P * y2P;
		coordT sinAB = sinA * cosB + cosA * sinB;
		if (sinAB < 0) {
			return true;
		}
		else
			return false;
	}

}

bool isInCircumcircle(int P, int P1, int P2, int P3, coordT** List)
{
	coordT centerX, centerY;
	coordT a1, b1, d1, a2, b2, d2;
	coordT radius, checkDistance;

	// Make the below calculation easier

	a1 = -2.0 * List[0][P1] + 2.0 * List[0][P2];
	b1 = -2.0 * List[1][P1] + 2.0 * List[1][P2];
	d1 = -List[0][P1] * List[0][P1] - List[1][P1] * List[1][P1] + List[0][P2] * List[0][P2] +
		 List[1][P2] * List[1][P2];
	a2 = -2.0 * List[0][P1] + 2.0 * List[0][P3];
	b2 = -2.0 * List[1][P1] + 2.0 * List[1][P3];
	d2 = -List[0][P1] * List[0][P1] - List[1][P1] * List[1][P1] + List[0][P3] * List[0][P3] +
		 List[1][P3] * List[1][P3];

	// Using determinant to calculate the coordinate x, and y of center point

	centerX = (d1 * b2 - d2 * b1) / (a1 * b2 - a2 * b1);
	centerY = (d2 * a1 - d1 * a2) / (a1 * b2 - a2 * b1);
	radius = sqrt(
        (centerX - List[0][P1]) * (centerX - List[0][P1]) +
        (centerY - List[1][P1]) * (centerY - List[1][P1]));

	// Distance between check point and center point

	checkDistance = sqrt(
		(centerX - List[0][P]) * (centerX - List[0][P]) +
		(centerY - List[1][P]) * (centerY - List[1][P]));
	return checkDistance <= radius;
}

//void flipEdge(int curTri, int adjTri, int L_index, int** Vertex, int** Edge)
//{
//	int P = Vertex[0][curTri];
//	// flip diagonal line by overlapping triangle i and adTri with 2 new triangle (net triangle + 0)
//
//	int edge1   = Edge[1][curTri];
//	int edge2   = Edge[2][curTri];
//	int adEdge1 = Edge[(L_index + 2) % 3][adjTri];
//	int adEdge2 = Edge[(L_index + 1) % 3][adjTri];
//	int V1      = Vertex[1][curTri];
//	int V2      = Vertex[2][curTri];
//	int L       = Vertex[L_index][adjTri]; // opposite with P by the sharing edge of
//										   // triangle curTri and adjTri
//
//	// triangle curTri
//
//	Vertex[0][curTri] = P;
//	Vertex[1][curTri] = V1;
//	Vertex[2][curTri] = L;
//	Edge[0][curTri]   = adEdge2;
//	Edge[1][curTri]   = adjTri;
//	Edge[2][curTri]   = edge2;
//
//	// triangle adjTri
//
//	Vertex[0][adjTri] = P;
//	Vertex[1][adjTri] = L;
//	Vertex[2][adjTri] = V2;
//	Edge[0][adjTri]   = adEdge1;
//	Edge[1][adjTri]   = edge1;
//	Edge[2][adjTri]   = curTri;
//}

std::map<std::string, QVariant> QhullPlugin::applyFilter(
	const QAction*           filter,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	vcg::CallBackPos* /* cb*/)
{
	qhT  qh_qh = {};
	qhT* qh    = &qh_qh;

	switch (ID(filter)) {
	case FP_DELAUNAY: {
		MeshModel& m         = *md.mm();
		MeshModel& nm        = *md.addNewMesh("", "Delaunay Triangulation");
		int        dim       = 3;
		int        numpoints = m.cm.vn;
		coordT*    points;                                // contain the coordinate of all the points
		points = readpointsFromMesh(&numpoints, &dim, m);
		
		int** Vertex = (int**) malloc(3 * sizeof(int*)); // contain the index of vertex
		int** Edge   = (int**) malloc(3 * sizeof(int*)); // contain index of adjacent triangle
		for (int i = 0; i < 3; i++) {
			Vertex[i] = (int*) malloc((2 * numpoints + 1) * sizeof(int));
			Edge[i]   = (int*) malloc((2 * numpoints + 1) * sizeof(int));
		}

		// create super triangle

		coordT minX, minY, maxX, maxY;
		minX = minY = 1e9;
		maxX = maxY = -1e9;
		for (int i = 0; i < numpoints; i++)
			for (int j = 0; j < dim; j++) {
				if (j == 0) { // dimension x
					minX = min(minX, points[i * dim + j]);
					maxX = max(maxX, points[i * dim + j]);
				}
				else if (j == 1) { // dimension y
					minY = min(minY, points[i * dim + j]);
					maxY = max(maxY, points[i * dim + j]);
				}
			}
		coordT dX = (maxX - minX) * 1;
		coordT dY = (maxY - minY) * 1;

		// coordinate of the super triangle vertecies

		coordT superX[3] = {minX - dX, maxX + dX * 3, minX - dX};
		coordT superY[3] = {minY - dY * 3, maxY + dY * 3, maxY + dY * 3};

		// create a list of all point

		coordT** List = (coordT**) malloc(2 * sizeof(coordT*));
		List[0]       = (coordT*)  malloc((numpoints + 3) * sizeof(coordT));
		List[1]       = (coordT*)  malloc((numpoints + 3) * sizeof(coordT));
		for (int i = 0; i < numpoints; i++) {
			List[0][i] = points[i * dim];
			List[1][i] = points[i * dim + 1];
		}

		// add super triangle vertecies to the list

		for (int i = numpoints; i < numpoints + 3; i++) {
			List[0][i] = superX[i - numpoints];
			List[1][i] = superY[i - numpoints];
		}

		//for (int i = 0; i < numpoints + 3; i++) {
		//	log("%f, %f, 0", List[0][i], List[1][i]);
		//}
	
		// adding super triangle to Vertex and Edge

		Vertex[0][1] = numpoints;
		Vertex[1][1] = numpoints + 1;
		Vertex[2][1] = numpoints + 2;
		Edge[0][1] = 0;
		Edge[1][1] = 0;
		Edge[2][1] = 0;

		// loop over each point
		int numTri = 1;
		for (int ii = 0; ii < 10; ii++) {
			int P = ii; // current point
			int T = triLocate(Vertex, Edge, P, List, numTri); // locate the enclosing triangle

			// overlap triangle T with one of 3 new triangle

			int A  = Edge[0][T];
			int B  = Edge[1][T];
			int C  = Edge[2][T];
			int V1 = Vertex[0][T];
			int V2 = Vertex[1][T];
			int V3 = Vertex[2][T];

			Vertex[0][T] = P;
			Vertex[1][T] = V2;
			Vertex[2][T] = V3;
			Edge[0][T]   = A;
			Edge[1][T]   = numTri + 1;
			Edge[2][T]   = numTri + 2;

			// create 2 new triangle (total net triangle + 2)
			
			// triangle numTri + 1

			numTri++;
			Vertex[0][numTri] = P;
			Vertex[1][numTri] = V3;
			Vertex[2][numTri] = V1;
			Edge[0][numTri]   = B;
			Edge[1][numTri]   = numTri + 1;
			Edge[2][numTri]   = T;

			// triangle numTri + 2

			numTri++;
			Vertex[0][numTri] = P;
			Vertex[1][numTri] = V1;
			Vertex[2][numTri] = V2;
			Edge[0][numTri]   = C;
			Edge[1][numTri]   = T;
			Edge[2][numTri]   = numTri - 1;

			log("succesfully create 2 new triangle");

			/* creating queue and put all above triangle into queue*/

			vector<int> queue;
			queue.push_back(T);
			queue.push_back(numTri - 1);
			queue.push_back(numTri);
			int size = queue.size();

			for (int cnt = 0; cnt < size; cnt++) {
				log("%d", cnt);
				int curTri = queue[cnt]; // current triangle
				int adjTri = Edge[0][curTri]; // adjacent triangle through point P of current triangle
				if (adjTri == 0) { // skip if this triangle has no adjacent triangle (oppsite to point P)
					continue;
				}

				// find index of point L opposite with P by the sharing edge of triangle curTri and adjTri

				int L_index = 0;
				for (L_index = 0; L_index < 3; L_index++) {
					if (Vertex[L_index][adjTri] != Vertex[1][curTri] && Vertex[L_index][adjTri] != Vertex[2][curTri]) {
						break;
					}
				}
				if (L_index >= 3)
					continue;
				//flipEdge(curTri, adjTri, L_index, Vertex, Edge);
				
				if (isInCircumcircle(P, Vertex[1][curTri], Vertex[2][curTri], Vertex[L_index][adjTri], List))
				{
					log("point p is on the circle");
					
					// flip diagonal line by overlapping triangle i and adTri with 2 new triangle (net triangle + 0)

					int edge1   = Edge[1][curTri];
					int edge2   = Edge[2][curTri];
					int adEdge1 = Edge[(L_index + 2) % 3][adjTri];
					int adEdge2 = Edge[(L_index + 1) % 3][adjTri];
					int V1      = Vertex[1][curTri];
					int V2      = Vertex[2][curTri];
					int L       = Vertex[L_index][adjTri]; // opposite with P by the sharing edge of
														  // triangle curTri and adjTri

					// triangle curTri

					Vertex[0][curTri] = P;
					Vertex[1][curTri] = V1;
					Vertex[2][curTri] = L;
					Edge[0][curTri]   = adEdge2;
					Edge[1][curTri]   = adjTri;
					Edge[2][curTri]   = edge2;

					// triangle adjTri

					Vertex[0][adjTri] = P;
					Vertex[1][adjTri] = L;
					Vertex[2][adjTri] = V2;
					Edge[0][adjTri]   = adEdge1;
					Edge[1][adjTri]   = edge1;
					Edge[2][adjTri]   = curTri;

					// update edge for other adjacent triangle

					for (int i = 0; i < 3; i++) {
						if (Edge[i][edge1] == curTri) {
							Edge[i][edge1] = adjTri;
							break;
						}
					}
					for (int i = 0; i < 3; i++) {
						if (Edge[i][adEdge2] == adjTri) {
							Edge[i][adEdge2] = curTri;
							break;
						}
					}

					// add these 2 triangle to queue

					log("flip successed");
					queue.push_back(curTri);
					queue.push_back(adjTri);
					size += 2;
				}
			}
			queue.clear();
		}
		// adding face to mesh

		for (int i = 1; i < numTri; i++) {

			// skip face containing super triangle vertices

			if (Vertex[0][i] >= numpoints || Vertex[1][i] >= numpoints ||
				Vertex[2][i] >= numpoints) {
				continue;
			}

			Point3d p0 = {List[0][Vertex[0][i]], List[1][Vertex[0][i]], 0};
			Point3d p1 = {List[0][Vertex[1][i]], List[1][Vertex[1][i]], 0};
			Point3d p2 = {List[0][Vertex[2][i]], List[1][Vertex[2][i]], 0};
			tri::Allocator<CMeshO>::AddFace(nm.cm, p0, p1, p2);
		}

		// clear pointer

		for (int i = 0; i < 3; i++) {
			free(Vertex[i]);
			free(Edge[i]);
		}
		free(Vertex);
		free(Edge);
		free(List[0]);
		free(List[1]);
		free(List);

	} break;
	default: wrongActionCalled(filter);
	}

	return std::map<std::string, QVariant>();
}
MESHLAB_PLUGIN_NAME_EXPORTER(QhullPlugin)
