#include <iostream>
#include "utils/String.hpp"
#include "core/test/CoreTest.hpp"

using namespace mscut;

// Follow is your old code, I hope you can compare with my code instead of discarding them.
struct YourCode {
	void VoriDoCut(std::vector<Eigen::Vector3d>& face) {
		//double x = 0.0;
		//Eigen::MatrixXd rotate;
		//std::vector<Point_2> tri/* = getTriLocalCoord(face, x, rotate)*/;

		//Polygon_2 tri_boundary(tri.begin(), tri.end());

		//for (auto i : tri)
		//	//sites3d.emplace_back(getTriGlobalCoord(i, x, rotate));

		//	Apollonius_graph ag;
		//std::vector<Point_2> sites = tri;
		//for (int i = 0; i < sites.size(); ++i)
		//{
		//	const auto sitePoint = sites[i];
		//	Apollonius_graph::Site_2 site(sitePoint, 0);
		//	//ag.insert(site);
		//}

		//Iso_rectangle_2 bbox(-2000, -2000, 2000, 2000);

		//CroppedVoronoiFromApollonius vor(bbox);

		// iterate to extract Voronoi diagram edges around each vertex
		/*Apollonius_graph::Finite_vertices_iterator vit;*/
		//BoundaryApolloniusGraph bag(tri_boundary);
		//for (auto vit = ag.finite_vertices_begin(); vit != ag.finite_vertices_end(); ++vit) {
		//	std::cout << "Vertex " << vit->site().point() << std::endl;
		//	Apollonius_graph::Edge_circulator ec = ag.incident_edges(vit), done(ec);
		//	if (ec != 0) {
		//		do {
		//			ag.draw_dual_edge(*ec, vor);
		//		} while (++ec != done);
		//	}
		//	//print the cropped Voronoi diagram edges as segments
		//	/*std::copy(vor.m_cropped_vd.begin(), vor.m_cropped_vd.end(),
		//		std::ostream_iterator<Segment_2>(std::cout, "\n"));
		//	std::cout << "=========\n";*/
		//	for (const auto& vor_seg : vor.m_cropped_vd)
		//	{
		//		croppedAP.emplace_back(getTriGlobalCoord(vor_seg.source(), x, rotate));
		//		croppedAP.emplace_back(getTriGlobalCoord(vor_seg.target(), x, rotate));
		//		//ag_dual_segs.emplace_back(vor_seg);
		//		bag << vor_seg;
		//	}
		//	vor.reset();
		//}

		/*for (const auto& clip_seg : bag.clipSegments)
		{
			clipSegments.emplace_back(getTriGlobalCoord(clip_seg.vertex(0), x, rotate));

			clipSegments.emplace_back(getTriGlobalCoord(clip_seg.vertex(1), x, rotate));
		}*/

		std::cout << "###################\n";
		//extract the entire cropped Voronoi diagram
		//ag.draw_dual(vor);
	}

	void yourCodeInMain()
	{
		//for (int i = 0; i < faceMat.rows(); ++i) {
		//	auto face = faceMat.row(i);
		//	faces.emplace_back(face);
		//}

		//for (int i = 0; i < faces.size(); ++i)
		//{
		//	std::vector<Eigen::Vector3d> face;
		//	for (int j = 0; j < 3; ++j)
		//		face.emplace_back(verts[faces[i][j]]);
		//	VoriDoCut(face);
		//}

		//const std::string sites_out_file = R"(.\data\sites.obj)";
		//const std::string out_file = R"(.\data\\out_vor.obj)";
		//const std::string clip_out_file = R"(.\data\clip_out_vor.obj)";

		//std::ofstream sites_out(sites_out_file);
		////std::cout << sites3d.size() << std::endl;
		//for (auto i : sites3d)
		//	sites_out << "v " << i.transpose() << " \n";
		//sites_out.close();

		//std::ofstream out(out_file);
		//int edgeIdx = 1;
		//for (int i = 0; i < croppedAP.size(); ++i) {
		//	out << "v " << croppedAP[i].transpose() << " \nv " << croppedAP[++i].transpose() << std::endl;
		//	//std::cout << "v1: " << vor_seg.vertex(0) << " 0 v2: " << vor_seg.vertex(1) << " 0" << std::endl;
		//	out << "l " << edgeIdx << " " << edgeIdx + 1 << std::endl;
		//	edgeIdx += 2;
		//}
		//out.close();

		//std::ofstream clip_out(clip_out_file);
		//edgeIdx = 1;
		//for (int i = 0; i < clipSegments.size(); ++i) {
		//	out << "v " << clipSegments[i].transpose() << " \nv " << clipSegments[++i].transpose() << std::endl;
		//	//std::cout << "v1: " << vor_seg.vertex(0) << " 0 v2: " << vor_seg.vertex(1) << " 0" << std::endl;
		//	out << "l " << edgeIdx << " " << edgeIdx + 1 << std::endl;
		//	edgeIdx += 2;
		//}
		//clip_out.close();
	}
};

using ScalarFunc = core::MSCuttingModel::ScalarFunc;

Scalar val(const Eigen::Vector3d& p)
{
	return p.x() * p.x() + p.y() * p.y() + p.z() * p.z() - 1;
}

Eigen::Vector3d grad(const Eigen::Vector3d& p)
{
	return 2 * p;
}

int main()
{
	const std::string mesh_file(MODEL_DIR DELIMITER R"(singleTri.obj)");
	const static int numSamplePoints = 0;

	ScalarFunc scalarFunc = { val, grad };
	core::MSCuttingModel mscModel(mesh_file, scalarFunc, numSamplePoints);
	const std::string ad_vis_file = str_util::concatFilePath(VIS_DIR, mscModel.modelName, (std::string)"apollonius_diagram.obj");
	mscModel.launch(ad_vis_file);

	/*unit_test::testSamplingPoint(mesh_file, 2);

	unit_test::testLocalGlobalTransoform(mesh_file, numSamplePoints);

	unit_test::testComputeAD(mesh_file, numSamplePoints);*/

	return 0;
}