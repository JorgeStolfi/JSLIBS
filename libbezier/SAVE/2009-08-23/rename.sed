#! /bin/sed -f
# Last edited on 2005-02-12 07:52:48 by stolfi
# 
s/dg_bezier/bz_basic/g
s/dg_bezpatch/bz_patch/g
s/test_dg_bezier/test_bz/g
s/dg_BezIndex/bz_index_t/g
s/dg_interval/interval/g
s/dg_Interval/interval_t/g
s/bzn_/bz_/g
s/bz_DDim/bz_ddim_t/g
s/bz_RDim/bz_rdim_t/g
s/bz_AxisIndex/bz_axis_index_t/g
s/bz_Axis/bz_axis_t/g
s/bz_CPIndex/bz_cpindex_t/g
s/bz_Patch/bz_patch_t/g
s/bz_Degrees/bz_degrees_t/g
s/bz_Degree/bz_degree_t/g
s/bz_patch_t_/bz_patch_/g
s/bz_get_control_point/bz_patch_control_point/g
s/bz_eval/bz_patch_eval/g
s/dg_SignedDir/bz_signed_dir_t/g
s/dg_FaceIndex/bz_face_index_t/g

s/dg_point_box_map/box_point_map/g
s/dg_point_box_unmap/box_point_unmap/g
s/dg_axis_belongs/box_axis_belongs/g
s/dg_axis_include/box_axis_include/g
s/dg_axis_exclude/box_axis_exclude/g
s/dg_axes_count/box_axes_count/g
s/dg_axes_complement/box_axes_complement/g
s/dg_axis/box_axis/g
s/dg_find_axis/box_find_axis/g
s/dg_face_dimension/box_face_dimension/g
s/dg_face_position/box_face_position/g
s/dg_face_signature/box_face_signature/g
s/dg_face_rel_box/box_face_rel_box/g
s/dg_face_print/box_face_print/g


