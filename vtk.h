#ifndef COMMON_VTK_H
#define COMMON_VTK_H

template <typename OS, typename POINTMESH>
void point2vtk(OS &os, const POINTMESH &pmt)
{
    os << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";

	os << "POINTS " << pmt.pcoord_.size() << " float\n";
	for (size_t i = 0; i < pmt.pcoord_.size(); ++i)
		os << pmt.pcoord_[i][0] << " " << pmt.pcoord_[i][1] << " " << pmt.pcoord_[i][2] << "\n";

	os << "CELLS " << pmt.pcoord_.size() << " " << pmt.pcoord_.size() * 2 << "\n";
	for (size_t i = 0; i < pmt.pcoord_.size(); ++i)
        os << 1 << " " << i << "\n";

	os << "CELL_TYPES " << pmt.pcoord_.size() << "\n";
	for (size_t i = 0; i < pmt.pcoord_.size(); ++i)
        os << 1 << "\n";
}

template <typename OS, typename POINTMESH>
void line2vtk(OS &os, const POINTMESH &pmt)
{

    os << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";

    os << "POINTS " << pmt.pcoord_.size() << " float\n";
    for (size_t i = 0; i < pmt.pcoord_.size(); ++i)
        os << pmt.pcoord_[i][0] << " " << pmt.pcoord_[i][1] << " " << pmt.pcoord_[i][2] << "\n";

    os << "CELLS " << pmt.pcoord_.size()/2 << " " << pmt.pcoord_.size()/2 * 3 << "\n";
    for (size_t i = 0; i < pmt.pcoord_.size()/2; ++i)
        os << 2 << " " << 2*i+0 << " " << 2*i+1 << "\n";

    os << "CELL_TYPES " << pmt.pcoord_.size() / 2 << "\n";
    for (size_t i = 0; i < pmt.pcoord_.size() / 2; ++i)
        os << 3 << "\n";
}

template <typename OS, typename FLOAT, typename INT>
void point2vtk(OS &os,
	const FLOAT *node, size_t node_num,
	const INT *points, size_t points_num)
{
	os << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";

	os << "POINTS " << node_num << " float\n";
	for (size_t i = 0; i < node_num; ++i)
		os << node[i * 3 + 0] << " " << node[i * 3 + 1] << " " << node[i * 3 + 2] << "\n";

	os << "CELLS " << points_num << " " << points_num * 2 << "\n";
	for (size_t i = 0; i < points_num; ++i)
		os << 1 << " " << points[i] << "\n";

	os << "CELL_TYPES " << points_num << "\n";
	for (size_t i = 0; i < points_num; ++i)
		os << 1 << "\n";
}

template <typename OS, typename FLOAT, typename INT>
void line2vtk(
	OS &os,
	const FLOAT *node, size_t node_num,
	const INT *line, size_t line_num)
{
	os << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";

	os << "POINTS " << node_num << " float\n";
	for (size_t i = 0; i < node_num; ++i)
		os << node[i * 3 + 0] << " " << node[i * 3 + 1] << " " << node[i * 3 + 2] << "\n";

	os << "CELLS " << line_num << " " << line_num * 3 << "\n";
	for (size_t i = 0; i < line_num; ++i)
		os << 2 << " " << line[i * 2 + 0] << " " << line[i * 2 + 1] << "\n";

	os << "CELL_TYPES " << line_num << "\n";
	for (size_t i = 0; i < line_num; ++i)
		os << 3 << "\n";
}

/*
template <typename OS, typename TRIMESH>
void tri2vtk(OS &os, TRIMESH & tmp)
{
	os << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";

	os << "POINTS " << tmp.verts().size() << " float\n";

	std::map<type_define::CVI, size_t> cvi2id;
	int i = 0;
	for (auto it = tmp.verts().begin(); it != tmp.verts().end(); ++it){
		os << it->pos_[0] << " " << it->pos_[1] << " " << it->pos_[2] << "\n";
		cvi2id[it] = i++;
	}

	os << "CELLS " << tmp.faces().size() << " " << tmp.faces().size() * 4 << "\n";
	for (auto it = tmp.faces().begin(); it != tmp.faces().end(); ++it){
		os << 3;

		type_define::CEI ei = it->edge();
		for (size_t i = 0; i < 3; ++i){
			type_define::CVI vi = ei->vert();
			auto v_it = cvi2id.find(vi);
			if (v_it == cvi2id.end()){
				throw std::logic_error("can not find vi");
			}
			os << " " << v_it->second;
			ei = ei->next();
		}

		os << "\n";
	}
	os << "CELL_TYPES " << tmp.faces().size() << "\n";
	for (size_t i = 0; i < tmp.faces().size(); ++i)
		os << 5 << "\n";
}
*/
template <typename OS, typename FLOAT, typename INT>
void tri2vtk(
	OS &os,
	const FLOAT *node, size_t node_num,
	const INT *tri, size_t tri_num)
{
	os << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";

	os << "POINTS " << node_num << " float\n";
	for (size_t i = 0; i < node_num; ++i)
		os << node[i * 3 + 0] << " " << node[i * 3 + 1] << " " << node[i * 3 + 2] << "\n";

	os << "CELLS " << tri_num << " " << tri_num * 4 << "\n";
	for (size_t i = 0; i < tri_num; ++i)
		os << 3 << "  " << tri[i * 3 + 0] << " " << tri[i * 3 + 1] << " " << tri[i * 3 + 2] << "\n";
	os << "CELL_TYPES " << tri_num << "\n";
	for (size_t i = 0; i < tri_num; ++i)
		os << 5 << "\n";
}

template <typename OS, typename Iterator, typename INT>
void vtk_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table")
{
    os << "SCALARS " << value_name << " float\nLOOKUP_TABLE " << table_name << "\n";
    for(int i = 0; i < size; ++i, ++first)
        os << *first << "\n";
}

template <typename OS, typename Iterator, typename INT>
void vtk_data_rgba(OS &os, Iterator first, INT size, const char *value_name,
                   const char *table_name = "my_table")
{
    os << "COLOR_SCALARS " << value_name << " 4\n";//\nLOOKUP_TABLE " << table_name << "\n";
    for(size_t i = 0; i < size; ++i)
    {
        for(size_t j = 0; j < 4; ++j,++first)
        {
            if(j != 3)
                os << *first << " ";
            else
                os << *first;
        }
        os << "\n";
    }
}

template <typename OS, typename Iterator, typename INT>
void point_data(OS &os, Iterator first, INT size, char *value_name, char *table_name = "my_table")
{
	os << "POINT_DATA " << size << "\n";
	vtk_data(os, first, size, value_name, table_name);
}

template <typename OS, typename Iterator, typename INT>
void cell_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table")
{
    os << "CELL_DATA " << size << "\n";
    vtk_data(os, first, size, value_name, table_name);
}

template <typename OS, typename Iterator, typename INT>
void cell_data_rgba(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table")
{
    os << "CELL_DATA " << size << "\n";
    vtk_data_rgba(os, first, size, value_name, table_name);
}

template <typename OS, typename Iterator, typename INT>
void point_data_rgba(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table")
{
    os << "POINT_DATA " << size << "\n";
    vtk_data_rgba(os, first, size, value_name, table_name);
}

template <typename OS, typename Iterator, typename INT>
void cell_data_rgba_and_scalar(OS &os, Iterator rgba_first, Iterator scalar_first, INT size,
                               const char *rgba_value_name, const char *scalar_value_name,
                               const char *table_name = "my_table")
{
    os << "CELL_DATA " << size << "\n";
    vtk_data_rgba(os, rgba_first, size, rgba_value_name, table_name);
    vtk_data(os, scalar_first, size, scalar_value_name, table_name);

}


#endif // COMMON_VTK_H
