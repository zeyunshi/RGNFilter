#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>

#include "permutohedral.h"
#include "vtk.h"
#include "obj.h"


    void computeFaceCenters(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, Eigen::MatrixXd& centers) {
        centers.resize(3, tris.cols());
        for (int f = 0; f < tris.cols(); ++f)
            centers.col(f) = (pts.col(tris(0, f)) + pts.col(tris(1, f)) + pts.col(tris(2, f))) / 3.0;
    }

    void computeFaceNormals(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, Eigen::MatrixXd& nls, Eigen::VectorXd& areas) {
        nls.resize(3, tris.cols());
        areas.resize(tris.cols());
        for (int i = 0; i < tris.cols(); ++i) {
            Eigen::Vector3d v1 = pts.col(tris(1, i)) - pts.col(tris(0, i));
            Eigen::Vector3d v2 = pts.col(tris(2, i)) - pts.col(tris(0, i));
            Eigen::Vector3d normal = v1.cross(v2);
            const double len = normal.norm();
            areas[i] = 0.5 * len;
            if (len < 1e-8)
                nls.col(i) = Eigen::VectorXd::Zero(3, 1);
            else
                nls.col(i) = normal / len;
        }
    }

    void computePointNormals(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, 
        const Eigen::MatrixXd& nls, Eigen::MatrixXd& pt_nls) {
        pt_nls.resize(3, pts.cols());
        pt_nls.setZero();
        for (int i = 0; i < tris.cols(); ++i) {
            for (int t = 0; t < 3; ++t) {
                pt_nls.col(tris(t, i)) += nls.col(i);
            }
        }
        for (int i = 0; i < pts.cols(); ++i)
            pt_nls.col(i).normalize();
    }

    double computeAverageEdgeLength(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris) {
        double sum_len = 0;
        int edge_num = 0;
        for (int i = 0; i < tris.cols(); ++i) {
            for (int j = 0; j < 3; j++) {
                Eigen::Vector3d v1 = pts.col(tris(j, i));
                Eigen::Vector3d v2 = pts.col(tris((j + 1) % 3, i));
                sum_len += (v1 - v2).norm();
                ++edge_num;
            }
        }
        return (edge_num > 0 ? sum_len / edge_num : 0);
    }

    int rollingGuidanceNormalFiltering(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, 
        const Eigen::MatrixXd &centers, const Eigen::MatrixXd& nls, 
        const Eigen::VectorXd& areas, Eigen::MatrixXd& result_nls,
        const double sigma_r, const double sigma_s, const int iter_num) {
        const double avg_edge_len = computeAverageEdgeLength(pts, tris);
        const double EPSILON = 1e-20;
        
        Eigen::MatrixXd pos(6, tris.cols());
        pos.block(0, 0, 3, tris.cols()) = centers / (sigma_s * avg_edge_len);

        Eigen::MatrixXd vals = nls;
        for (int i = 0; i < tris.cols(); ++i)
            vals.col(i) *= areas[i];
        
        for (int i = 0; i < iter_num; ++i) {
            pos.block(3, 0, 3, pos.cols()) = result_nls / sigma_r;

            PermutohedralLattice lattice(pos.rows(), vals.rows(), pos.cols());
            std::vector<float> p(6), v(3);
            for (int i = 0; i < tris.cols(); ++i) {
                for (int t = 0; t < 3; ++t) {
                    p[t] = static_cast<float>(pos(t, i));
                    v[t] = static_cast<float>(vals(t, i));
                }
                for (int t = 3; t < 6; ++t)
                    p[t] = static_cast<float>(pos(t, i));
                lattice.splat(&p[0], &v[0]);
            }

            lattice.blur();

            lattice.beginSlice();
            for (int i = 0; i < tris.cols(); ++i) {
                lattice.slice(&v[0]);
                for (int t = 0; t < 3; ++t)
                    result_nls(t, i) = v[t];
            }

            for (int j = 0; j < result_nls.cols(); j++)
            {
                double d = result_nls.col(j).norm();
                if (d < EPSILON) 
                    result_nls.col(j) = nls.col(j);
                else 
                    result_nls.col(j) /= d;
            }
        }

        return 0;
    }

    void computeRotation(const Eigen::Vector3d& vec_a, const Eigen::Vector3d& vec_b, Eigen::Matrix3d& rot) {
        double tmp = vec_a.norm() * vec_b.norm();
        if (tmp < 1e-20) {
            rot.setIdentity();
            return;
        }
        Eigen::Vector3d axis = vec_a.cross(vec_b);
        const double len = axis.norm();
        if (len < 1e-20)
            axis = Eigen::Vector3d(0, 0, 1);
        else
            axis /= len;
        double cos_alpha = vec_a.dot(vec_b) / tmp;
        if (cos_alpha < -1) 
            cos_alpha = -1;
        if (cos_alpha > 1) 
            cos_alpha = 1;
        const double alpha = acos(cos_alpha);
        rot = Eigen::AngleAxisd(acos(cos_alpha), axis).toRotationMatrix();
    }

    void getTriangleVertexGradient(const Eigen::VectorXd& a, const Eigen::VectorXd& b, Eigen::VectorXd& gradient) {
        double len2 = 0;
        Eigen::VectorXd c = a - b;

        double dotres = a.dot(c);
        double lenc2 = c.squaredNorm();
        double ratio = dotres / (lenc2);

        gradient = (b - a) * ratio + a;

        len2 = gradient.squaredNorm();
        gradient /= len2;
    }

    double computeTriangleArea(const Eigen::VectorXd& source, const Eigen::VectorXd& vleft, const Eigen::VectorXd& vright) {
        Eigen::Vector3d v1 = vleft - source;
        Eigen::Vector3d v2 = vright - source;
        Eigen::Vector3d nls = v1.cross(v2);
        return 0.5 * nls.norm();
    }

    void computeTriangleDivergence(const Eigen::VectorXd& source, const Eigen::VectorXd& vleft, const Eigen::VectorXd& vright,
        double &area, Eigen::VectorXd& div) {
        Eigen::VectorXd s_l = source - vleft;
        Eigen::VectorXd s_r = source - vright;
        Eigen::VectorXd l_r = vleft - vright;
        Eigen::VectorXd l_s = -s_l;
        Eigen::VectorXd r_s = -s_r;
        Eigen::VectorXd r_l = -l_r;

        //����iT
        Eigen::VectorXd hs;  getTriangleVertexGradient(s_l, s_r, hs);
        Eigen::VectorXd hl;  getTriangleVertexGradient(l_r, l_s, hl);
        Eigen::VectorXd hr;  getTriangleVertexGradient(r_s, r_l, hr);

        //gradient field
        Eigen::VectorXd wx = hl * (l_s[0]) + hr * (r_s[0]);
        Eigen::VectorXd wy = hl * (l_s[1]) + hr * (r_s[1]);
        Eigen::VectorXd wz = hl * (l_s[2]) + hr * (r_s[2]);

        //S��
        area = computeTriangleArea(source, vleft, vright);

        //divergence
        div.resize(3);
        div[0] = wx.dot(hs) * area;
        div[1] = wy.dot(hs) * area;
        div[2] = wz.dot(hs) * area;
    }

    void computeDivergence(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, const Eigen::MatrixXd &centers,
        const Eigen::MatrixXd& origin_nls, const Eigen::MatrixXd& filtered_nls,  
        Eigen::MatrixXd& divs, std::vector<bool>& fixed) {
        divs.resize(pts.cols(), 3);
        divs.setZero();
        Eigen::Matrix3d rot, rot_pts;
        rot.setIdentity();
        Eigen::VectorXd div(3);
        std::vector<double> pt_areas(pts.cols(), 0);
        double area;
        for (int f = 0; f < tris.cols(); ++f) {
            computeRotation(origin_nls.col(f), filtered_nls.col(f), rot);
            for (int t = 0; t < 3; ++t)
                rot_pts.col(t) = rot * (pts.col(tris(t, f)) - centers.col(f)) + centers.col(f);
            for (int t = 0; t < 3; ++t) {
                const int t1 = (t + 1) % 3;
                computeTriangleDivergence(rot * rot_pts.col(t), rot * rot_pts.col(t1),
                    rot * rot_pts.col((t + 2) % 3), area, div);
                for (int j = 0; j < 3; ++j)
                    divs(tris(t1, f), j) -= div[j];
                pt_areas[tris(t1, f)] += area;
            }
        }

        for (int i = 0; i < pts.cols(); ++i) {
            for (int t = 0; t < 3; ++t)
                divs(i, t) /= pt_areas[i];
        }

        // fix first point 
        for (int i = 0; i < 1; ++i) {
            fixed[i] = true;
            for (int t = 0; t < 3; ++t)
                divs(i, t) = pts(t, i);
        }

        // fix first triangle 
        //for (int f = 0; f < tris.cols() / 5; ++f) {
        //    for (int t = 0; t < 3; ++t) {
        //        fixed[tris(t, f)] = true;
        //        for (int d = 0; d < 3; ++d) {
        //            divs(tris(t, f), d) = pts(d, tris(t, f));
        //        }
        //    }
        //}

        //for (int i = 0; i < pts.cols(); ++i) {
        //    bool too_large = false;
        //    for (int t = 0; t < 3; ++t)
        //        if (fabs(divs(i, t)) > 0.02)
        //            too_large = true;
        //    if (too_large) {
        //        fixed[i] = true;
        //        for (int t = 0; t < 3; ++t)
        //            divs(i, t) = pts(t, i);
        //    }
        //}

        // fix some special points 
        //std::vector<int> fixed_pts = { 264,871,206,456,536,278 };
        //for (int i = 0; i < fixed_pts.size(); ++i) {
        //    fixed[fixed_pts[i]] = true;
        //    for (int d = 0; d < 3; ++d)
        //        divs(fixed_pts[i], d) = pts(d, fixed_pts[i]);
        //}
    }

    void rotationDisconnected(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, const Eigen::MatrixXd &centers,
        const Eigen::MatrixXd& origin_nls, const Eigen::MatrixXd& filtered_nls) {
        std::ofstream ofs("divergence.obj");

        Eigen::Matrix3d rot;
        rot.setIdentity();
        for (int f = 0; f < tris.cols(); ++f) {
            computeRotation(origin_nls.col(f), filtered_nls.col(f), rot);
            for (int t = 0; t < 3; ++t) {
                Eigen::VectorXd pt = rot * pts.col(tris(t, f)) - (rot * centers.col(f) - centers.col(f));
                ofs << "v " << pt[0] << ' ' << pt[1] << ' ' << pt[2] << '\n';
            }
        }
        for (int f = 0; f < tris.cols(); ++f) {
            ofs << "f " << 3 * f + 1 << ' ' << 3 * f + 2 << ' ' << 3 * f + +3 << '\n';
        }
        ofs.close();
    }

    void computeUniformLaplacian3(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, Eigen::SparseMatrix<double>& L) {
        std::vector<Eigen::Triplet<double>> triple;
        triple.reserve(pts.cols() * 7 * 3);

        std::vector<std::set<int> > adj_map(pts.cols());
        for (int f = 0; f < tris.cols(); ++f) {
            for (int t = 0; t < 3; ++t) {
                adj_map[tris(t, f)].insert(tris((t + 1) % 3, f));
                adj_map[tris((t + 1) % 3, f)].insert(tris(t, f));
            }
        }

        for (int i = 0; i < pts.cols(); ++i) {
            const std::set<int>& neighbors = adj_map[i];
            const int vi3 = 3 * i;
            double r = -1.0 / neighbors.size();
            for (auto itr = neighbors.begin(); itr != neighbors.end(); ++itr) {
                const int vj3 = 3 * *itr;
                for (int t = 0; t < 3; ++t)
                    triple.emplace_back(vi3 + t, vj3 + t, r);
            }
            for (int t = 0; t < 3; ++t)
                triple.emplace_back(vi3 + t, vi3 + t, 1);
        }

        L.resize(pts.cols() * 3, pts.cols() * 3);
        L.setFromTriplets(triple.begin(), triple.end());
    }

    void computeUniformLaplacian(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, 
        const std::vector<bool> &fixed, 
        Eigen::SparseMatrix<double>& L) {
        std::vector<Eigen::Triplet<double>> triple;
        triple.reserve(pts.cols() * 7);

        std::vector<std::set<int> > adj_map(pts.cols());
        for (int f = 0; f < tris.cols(); ++f) {
            for (int t = 0; t < 3; ++t) {
                adj_map[tris(t, f)].insert(tris((t + 1) % 3, f));
                adj_map[tris((t + 1) % 3, f)].insert(tris(t, f));
            }
        }

        for (int i = 0; i < pts.cols(); ++i) {
            if (fixed[i]) {
                triple.emplace_back(i, i, 1);
                continue;
            }
            const std::set<int>& neighbors = adj_map[i];
            double r = -1.0 / neighbors.size();
            for (auto itr = neighbors.begin(); itr != neighbors.end(); ++itr) {
                triple.emplace_back(i, *itr, r);
            }
            triple.emplace_back(i, i, 1);
        }

        L.resize(pts.cols(), pts.cols());
        L.setFromTriplets(triple.begin(), triple.end());
    }

    double sinValue(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
    {
        double lab = (b - a).norm();
        double lac = (c - a).norm();
        return ((b - a).cross(c - a)).norm() / (lab * lac);
    }

    double cosValue(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
    {
        double lab = (b - a).norm();
        double lac = (c - a).norm();
        double lbc = (b - c).norm();
        double lab2 = lab * lab;
        double lac2 = lac * lac;
        double lbc2 = lbc * lbc;
        return (lab2 + lac2 - lbc2) / (2.0 * lab * lac);
    }

    double maxValue(const double value1, const double value2)
    {
        if (value1 > value2)
            return value1;
        else
            return value2;
    }

    double computeCotWeight(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c) {
        double cosx = cosValue(a, b, c);
        double sinx = maxValue(sinValue(a, b, c), 1e-8);
        double cotx = cosx / sinx;
        return maxValue(cotx, 1e-8);
    }

    void computeCotCofficients(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris,
        std::vector<double>& voronoi_areas, std::vector<std::map<int, double> >& cot_map) {

    }

    void computeCotLaplacian(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris,
        const std::vector<bool>& fixed,
        Eigen::SparseMatrix<double>& L) {
        std::vector<Eigen::Triplet<double>> triple;
        triple.reserve(pts.cols() * 7);

        std::vector<double> voronoi_areas(pts.cols(), 0);
        std::vector<std::map<int, double> > cot_map(pts.cols());
        for (int f = 0; f < tris.cols(); ++f) {
            for (int t = 0; t < 3; ++t) {
                const int v0 = tris(t, f), v1 = tris((t + 1) % 3, f), v2 = tris((t + 2) % 3, f);
                double weight = computeCotWeight(pts.col(v2), pts.col(v0), pts.col(v1));
                cot_map[v0][v1] += 0.5*weight;
                cot_map[v1][v0] += 0.5*weight;
                voronoi_areas[v0] += (pts.col(v0) - pts.col(v2)).squaredNorm() * weight * 0.125;
                voronoi_areas[v1] += (pts.col(v1) - pts.col(v2)).squaredNorm() * weight * 0.125;
            }
        }

        for (int i = 0; i < pts.cols(); ++i) {
            if (fixed[i]) {
                triple.emplace_back(i, i, 1);
                continue;
            }
            const std::map<int, double>& neighbors = cot_map[i];
            double sum = 0, w = 0;
            for (auto itr = neighbors.begin(); itr != neighbors.end(); ++itr) {
                w = itr->second / voronoi_areas[i];
                triple.emplace_back(i, itr->first, -w);
                sum += w;
            }
            triple.emplace_back(i, i, sum);
        }

        L.resize(pts.cols(), pts.cols());
        L.setFromTriplets(triple.begin(), triple.end());
    }

    void generateIdentityMatrix(const int dim, Eigen::SparseMatrix<double>& I) {
        I.resize(dim, dim);
        std::vector<Eigen::Triplet<double>> triple;
        triple.reserve(dim);
        for (int i = 0; i < dim; i++)
            triple.emplace_back(i, i, 1);
        I.setFromTriplets(triple.begin(), triple.end());
        I.makeCompressed();
    }

    void computeNormalConstraint(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, const Eigen::MatrixXd &nls,
        const int dim, Eigen::SparseMatrix<double>& C, Eigen::VectorXd& b) {
        std::vector<Eigen::Triplet<double>> triple;
        triple.reserve(pts.cols());
        b.resize(pts.cols());
        b.setZero();
        for (int i = 0; i < pts.cols(); ++i) {
            const Eigen::VectorXd& n = nls.col(i);
            const Eigen::VectorXd& p = pts.col(i);
            triple.emplace_back(i, i, n[dim]);
            b[i] += p[dim] * n[dim];
        }
        C.resize(pts.cols(), pts.cols());
        C.setFromTriplets(triple.begin(), triple.end());
    }

    int vertexUpdatingLaplace(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, 
        const Eigen::MatrixXd &centers, const Eigen::MatrixXd& origin_nls,
        const Eigen::MatrixXd& filtered_nls, Eigen::MatrixXd& result_pts) {
        // generate poisson structure 
        Eigen::MatrixXd divs;
        std::vector<bool> fixed(pts.cols(), false);
        computeDivergence(pts, tris, centers, origin_nls, filtered_nls, divs, fixed);
        Eigen::SparseMatrix<double> L;
        computeCotLaplacian(pts, tris, fixed, L);

        if (1) {
            std::ofstream ofs("divergencex.vtk");
            tri2vtk(ofs, &pts(0, 0), pts.cols(), &tris(0, 0), tris.cols());
            Eigen::VectorXd xd = divs.col(0);
            point_data(ofs, &xd[0], pts.cols(), "xd");
        }
        if (1) {
            std::ofstream ofs("laplacex.vtk");
            tri2vtk(ofs, &pts(0, 0), pts.cols(), &tris(0, 0), tris.cols());
            Eigen::VectorXd xl = L * (pts.row(0).transpose());
            point_data(ofs, &xl[0], pts.cols(), "xl");
        }
        if (1) {
            std::ofstream ofs("diff.vtk");
            tri2vtk(ofs, &pts(0, 0), pts.cols(), &tris(0, 0), tris.cols());
            Eigen::VectorXd xl = L * (pts.row(0).transpose()) - divs.col(0);
            for (int i = 0; i < xl.size(); ++i)
                xl[i] = fabs(xl[i]);
            point_data(ofs, &xl[0], pts.cols(), "xdiff");
        }

        // solve poisson equation 
        { 
            // with fixed points constraints 
            Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
            L.makeCompressed();
            solver.compute(L);
            result_pts = pts;
            for (int t = 0; t < 3; ++t) {
                Eigen::VectorXd b = divs.col(t);
                result_pts.row(t) = solver.solve(b).transpose();
            }
        }

        return 0;
    }

    void computeFaceGradient(const Eigen::Vector3d &v0, const Eigen::Vector3d &v1, const Eigen::VectorXd &v2, Eigen::Matrix3d& grad) {
        Eigen::Vector3d e0 = v1 - v0; 
        Eigen::Vector3d e1 = v2 - v1;
        Eigen::Vector3d e2 = v0 - v2;

        Eigen::Vector3d n = e0.cross(e1);
        double len2 = n.squaredNorm();
        if (len2 < 1e-20) // for degenerated case
        {
            for (int i = 0; i < 3; i++)
                grad.col(i) = Eigen::Vector3d::Zero();
        }
        else
        {
            grad.col(0) = n.cross(e1) / len2;
            grad.col(1) = n.cross(e2) / len2;
            grad.col(2) = n.cross(e0) / len2;
        }
    }

    void computeGradientMatrix(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, Eigen::SparseMatrix<double>& G) {
        std::vector<Eigen::Triplet<double>> triple;
        triple.reserve(9 * tris.cols());
        Eigen::Matrix3d grad;
        for (int i = 0; i < tris.cols(); ++i) {
            computeFaceGradient(pts.col(tris(0, i)), pts.col(tris(1, i)), pts.col(tris(2, i)), grad);
            for (int k = 0; k < 3; ++k) {
                for (int j = 0; j < 3; ++j) {
                    triple.emplace_back(3 * i + j, tris(k, i), grad(j, k));
                }
            }
        }

        G.resize(3 * tris.cols(), pts.cols());
        G.reserve(Eigen::VectorXi::Constant(3 * tris.cols(), 3));
        G.setFromTriplets(triple.begin(), triple.end());
    }

    void computeMassMatrix(const Eigen::VectorXd& areas, Eigen::SparseMatrix<double>& M) {
        std::vector<Eigen::Triplet<double>> triple;
        triple.reserve(3 * areas.size());
        for (int i = 0; i < areas.size(); ++i) {
            for (int t = 0; t < 3; ++t)
                triple.emplace_back(3*i+t, 3*i+t, areas[i]);
        }
        M.resize(3 * areas.size(), 3 * areas.size());
        M.setFromTriplets(triple.begin(), triple.end());
        const double avg_area = areas.sum() / areas.size();
        M /= avg_area;
    }

    void computeDivMatrix(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris,
        const Eigen::MatrixXd& origin_nls, const Eigen::MatrixXd& filtered_nls, Eigen::MatrixXd& Dv) {
        Dv.resize(3 * tris.cols(), 3);

        Eigen::Matrix3d rot, tri_pts, rot_pts;
        Eigen::Matrix3d grad;
        rot.setIdentity();
        for (int f = 0; f < tris.cols(); ++f) {
            computeRotation(origin_nls.col(f), filtered_nls.col(f), rot);
            for (int t = 0; t < 3; ++t) {
                tri_pts.col(t) = pts.col(tris(t, f));
            }
            rot_pts = rot * tri_pts;
            computeFaceGradient(rot_pts.col(0), rot_pts.col(1), rot_pts.col(2), grad);
            Dv.block(3 * f, 0, 3, 3) = grad * rot_pts.transpose();
        }
    }

    int vertexUpdatingGradient(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris,
        const Eigen::VectorXd &areas, const Eigen::MatrixXd& origin_nls,
        const Eigen::MatrixXd& filtered_nls, Eigen::MatrixXd& result_pts) {
        Eigen::SparseMatrix<double> G;
        computeGradientMatrix(pts, tris, G);
        Eigen::SparseMatrix<double> Gt = G.transpose();
        
        Eigen::SparseMatrix<double> M;
        computeMassMatrix(areas, M);

        Eigen::SparseMatrix<double> LV = Gt * M * G;


        Eigen::SparseMatrix<double> I;
        generateIdentityMatrix(pts.cols(), I);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        double wl = 3000;
        Eigen::SparseMatrix<double> A = wl * LV + I;
        solver.analyzePattern(A);
        solver.factorize(A);

        Eigen::MatrixXd Dv;
        computeDivMatrix(pts, tris, origin_nls, filtered_nls, Dv);
        Eigen::MatrixXd deltaV = Gt * M * Dv;

        Eigen::MatrixXd b = pts.transpose() + wl * deltaV;
        Eigen::MatrixXd x = solver.solve(b);
        
        result_pts = x.transpose();
        return 0;
    }

    int tri2vtk_with_normal(const char *file_name, 
        const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, const Eigen::MatrixXd& nls) {
        std::ofstream ofs(file_name);
        if (ofs.fail())
            return __LINE__;
        tri2vtk(ofs, &pts(0, 0), pts.cols(), &tris(0, 0), tris.cols());
        ofs << "CELL_DATA " << tris.cols() << "\n";
        ofs << "NORMALS  " << "normals" << " float\n";
        for (int i = 0; i < tris.cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                if (2 != j)
                    ofs << nls(j, i) << ' ';
                else
                    ofs << nls(j, i);
            }
            ofs << "\n";
        }
        return 0;
    }

    int rgdSmoothing(const Eigen::MatrixXd& pts, const Eigen::MatrixXi& tris, Eigen::MatrixXd& result_pts) {
        Eigen::MatrixXd nls;
        Eigen::VectorXd areas;
        computeFaceNormals(pts, tris, nls, areas);
        //tri2vtk_with_normal("origin_normal.vtk", pts, tris, nls);

        Eigen::MatrixXd centers;
        computeFaceCenters(pts, tris, centers);

        Eigen::MatrixXd  filtered_nls(3, tris.cols());
        filtered_nls.setZero();
        const double sigma_r = 0.2, sigma_s = 5.0;
        const int iter_num = 5;
        rollingGuidanceNormalFiltering(pts, tris, centers, nls, areas, filtered_nls, sigma_r, sigma_s, iter_num);
        
        //tri2vtk_with_normal("filtered_normal.vtk", pts, tris, filtered_nls);

        //rotationDisconnected(pts, tris, centers, nls, filtered_nls);

        int ret = vertexUpdatingGradient(pts, tris, areas, nls, filtered_nls, result_pts);
        //tri2vtk_with_normal("rgdsmoothing.vtk", result_pts, tris, filtered_nls);
        return ret;
    }

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "#[Usage] RGNFilter input.obj output.obj" << std::endl;
        return __LINE__;
    }
    std::ifstream ifs(argv[1]);
    if (ifs.fail()) {
        std::cerr << "#[error] can't open : " << argv[1] << std::endl;
        return __LINE__;
    }
    Eigen::MatrixXd pts;
    Eigen::MatrixXi tris;
    if (0 != obj2tri(ifs, pts, tris)) {
        std::cerr << "#[error] obj format error!" << std::endl;
        return __LINE__;
    }
    
    Eigen::MatrixXd result_pts;
    clock_t start, end;
    start = clock();
    rgdSmoothing(pts, tris, result_pts);
    end = clock();
    std::cout << "Filter time : " << (double)(end - start) / 1000.0 << " sec" << std::endl;

    std::ofstream ofs(argv[2]);
    if (ofs.fail()) {
        std::cerr << "#[error] can't open : " << argv[2] << std::endl;
        return __LINE__;
    }
    tri2obj(ofs, &result_pts(0, 0), result_pts.cols(), &tris(0, 0), tris.cols());

    std::cerr << "#[info] completed..." << std::endl;

    return 0;
}