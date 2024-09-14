
/*---------------------------------------------------------------------------
	ZaRan	-	A Totally Automatic CFD Software
	Copyright (C) ,Since 2020
-------------------------------------------------------------------------------
License
	This file is part of ZaRan.

!	@file	geometry_types.h
!	@brief	CGAL link typedefs.
!	@author	Liu Guangying.
!   @date  2024.05.08
!   @location ShenZhen
\*---------------------------------------------------------------------------*/
#ifndef CGAL_TYPES_H // 防止重复包含
#define CGAL_TYPES_H

//-------------------------------------------------------------------------------Alpha_shape部分
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
//-------------------------------------------------------------------------------cdt部分
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

//-------------------------------------------------------------------------------Alpha_shape部分
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>         Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>               Delaunay;
typedef CGAL::Alpha_shape_2<Delaunay>                        Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator         Alpha_shape_vertices_iterator;
typedef Alpha_shape_2::Point_iterator                        Point_iterator;
//-----------------------------------------------------------------------------cdt部分
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb_cdt;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb_cdt;
typedef CGAL::Triangulation_data_structure_2<Vb_cdt, Fb_cdt> Tds_cdt;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds_cdt> CDT;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Finite_faces_iterator Finite_faces_iterator;



#endif // CGAL_TYPES_H