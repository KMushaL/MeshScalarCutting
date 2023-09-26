#include <iostream>
#include "ag2_exact_traits.h"

std::vector<Eigen::Vector3d> croppedAP;
std::vector<Eigen::Vector3d> sites3d;
std::vector<Eigen::Vector3d> clipSegments;

std::vector<Point_2> threeTurnTwo(const std::vector<Eigen::Vector3d>& face, double& x, Eigen::MatrixXd& rotate) {
	std::vector<Point_2> tri;
	Eigen::Vector3d v1 = face[0], v2 = face[1], v3 = face[2];

	Eigen::Vector3d xN = (v2 - v1).cross(v3 - v1);
	xN.normalize();
	Eigen::Vector3d yN = (v2 - v1);
	yN.normalize();
	Eigen::Vector3d zN = xN.cross(yN);
	zN.normalize();
	rotate << xN, yN, zN;

	v1 = rotate.transpose() * v1;
	v2 = rotate.transpose() * v2;
	v3 = rotate.transpose() * v3;

	tri.emplace_back(Point_2(v1[1], v1[2]));
	tri.emplace_back(Point_2(v2[1], v2[2]));
	tri.emplace_back(Point_2(v3[1], v3[2]));
	x = v1[0];
	return tri;
}

Eigen::Vector3d twoTurnThree(Point_2 vert, double x, Eigen::MatrixXd rotate) {
	Eigen::Vector3d vert3D = rotate * Eigen::Vector3d(x, vert.x(), vert.y());
	std::cout << vert3D.transpose()<<std::endl;
	return vert3D;
}

void VoriDoCut(std::vector<Eigen::Vector3d>& face) {
	double x = 0.0;
	Eigen::MatrixXd rotate;
	std::vector<Point_2> tri = threeTurnTwo(face, x, rotate);

	Polygon_2 tri_boundary(tri.begin(), tri.end());

	for (auto i : tri)
		sites3d.emplace_back(twoTurnThree(i, x, rotate));

	Apollonius_graph ag;
	std::vector<Point_2> sites = tri;
	for (int i = 0; i < sites.size(); ++i)
	{
		const auto sitePoint = sites[i];
		Apollonius_graph::Site_2 site(sitePoint, 0);
		ag.insert(site);
	}

	Iso_rectangle_2 bbox(-2000, -2000, 2000, 2000);

	Cropped_voronoi_from_apollonius vor(bbox);

	// iterate to extract Voronoi diagram edges around each vertex
	/*Apollonius_graph::Finite_vertices_iterator vit;*/
	BoundaryApolloniusGraph bag(tri_boundary);
	for (auto vit = ag.finite_vertices_begin(); vit != ag.finite_vertices_end(); ++vit) {
		std::cout << "Vertex " << vit->site().point() << std::endl;
		Apollonius_graph::Edge_circulator ec = ag.incident_edges(vit), done(ec);
		if (ec != 0) {
			do {
				ag.draw_dual_edge(*ec, vor);
			} while (++ec != done);
		}
		//print the cropped Voronoi diagram edges as segments
		/*std::copy(vor.m_cropped_vd.begin(), vor.m_cropped_vd.end(),
			std::ostream_iterator<Segment_2>(std::cout, "\n"));
		std::cout << "=========\n";*/
		for (const auto& vor_seg : vor.m_cropped_vd)
		{
			croppedAP.emplace_back(twoTurnThree(vor_seg.source(),x,rotate));
			croppedAP.emplace_back(twoTurnThree(vor_seg.target(), x, rotate));
			//ag_dual_segs.emplace_back(vor_seg);
			bag << vor_seg;
		}
		vor.reset();
	}

	for (const auto& clip_seg : bag.clipSegments)
	{
		clipSegments.emplace_back(twoTurnThree(clip_seg.vertex(0), x, rotate));

		clipSegments.emplace_back(twoTurnThree(clip_seg.vertex(1), x, rotate));
	}

	std::cout << "###################\n";
	//extract the entire cropped Voronoi diagram
	//ag.draw_dual(vor);
}


int main()
{
	std::string inPath("data/mesh.obj");
	Eigen::MatrixXd vertMat;
	Eigen::MatrixXi faceMat;

	igl::read_triangle_mesh(inPath, vertMat, faceMat);

	std::vector<Eigen::Vector3d> verts;
	std::vector<Eigen::Vector3i>faces;
	for (int i = 0; i < vertMat.rows(); ++i) {
		Eigen::Vector3d vertPos;
		vertPos = vertMat.row(i);
		verts.emplace_back(vertPos);
	}

	for (int i = 0; i < faceMat.rows(); ++i) {
		auto face = faceMat.row(i);
		faces.emplace_back(face);
	}

	for (int i = 0;i<faces.size();++i)
	{
		std::vector<Eigen::Vector3d> face;
		for(int j = 0;j<3;++j)
			face.emplace_back(verts[faces[i][j]]);
		VoriDoCut(face);
	}

	const std::string sites_out_file = "C:\\temp\\demo1\\data\\sites.obj";
	const std::string out_file = "C:\\temp\\demo1\\data\\out_vor.obj";
	const std::string clip_out_file = "C:\\temp\\demo1\\data\\clip_out_vor.obj";

	std::ofstream sites_out(sites_out_file);
	//std::cout << sites3d.size() << std::endl;
	for (auto i : sites3d)
		sites_out << "v " << i.transpose() << " \n";
	sites_out.close();

	std::ofstream out(out_file);
	int edgeIdx = 1;
	for (int i = 0; i < croppedAP.size(); ++i) {
		out << "v " << croppedAP[i].transpose() << " \nv " << croppedAP[++i].transpose() << std::endl;
		//std::cout << "v1: " << vor_seg.vertex(0) << " 0 v2: " << vor_seg.vertex(1) << " 0" << std::endl;
		out << "l " << edgeIdx << " " << edgeIdx + 1 << std::endl;
		edgeIdx += 2;
	}
	out.close();

	std::ofstream clip_out(clip_out_file);
	edgeIdx = 1;
	for (int i = 0; i < clipSegments.size(); ++i) {
		out << "v " << clipSegments[i].transpose() << " \nv " << clipSegments[++i].transpose() << std::endl;
		//std::cout << "v1: " << vor_seg.vertex(0) << " 0 v2: " << vor_seg.vertex(1) << " 0" << std::endl;
		out << "l " << edgeIdx << " " << edgeIdx + 1 << std::endl;
		edgeIdx += 2;
	}
	clip_out.close();

	return 0;
}